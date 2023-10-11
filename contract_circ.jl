function gen_tn_sc3d_circ(dem::DEM{3}, max_bond_dim::Int; svd_eps::Float64=1e-16, do_gauging::Bool=false, canonicalness_target::Float64=1e-4)
    n_detectors = size(dem.H, 1)
   
    # the X and Y positions in stim take steps of 2
    pos_checks = [convert(Tuple{Int,Int,Int}, round.((0.5*pos[1], 0.5*pos[2], pos[3]))) for pos in dem.pos_checks]
    xmin,xmax = extrema([pos[1] for pos in pos_checks])
    ymin,ymax = extrema([pos[2] for pos in pos_checks])
    zmin,zmax = extrema([pos[3] for pos in pos_checks])
    println(dem.n_bits, " ", n_detectors)
    println("$xmin:$xmax $ymin:$ymax $zmin:$zmax")

    # build empty 3D tensor network
    tn3d = generate_empty_cubictn(xmin, xmax, ymin, ymax, zmin, zmax)
    for j = 1:n_detectors
        idx = Index(2)
        pos = pos_checks[j]
        @assert length(open_indices(tn3d,pos)) == 1
        T = ITensor([1., 0.], idx, open_indices(tn3d,pos)[1])
        contract1!(tn3d, pos, T)
    end
    open_index = Dict{Pos3D,Index}()
    for x=xmin:xmax for y=ymin:ymax for z=zmin:zmax
        @assert length(open_indices(tn3d,(x,y,z))) == 1
        open_index[x,y,z] = open_indices(tn3d,(x,y,z))[1]
    end end end

    expval = 0
    # add faults 
    for i = 1:dem.n_bits
        println("$i / ",dem.n_bits)
        affected_dets = [j for j=1:n_detectors if dem.H[j,i]==1]
        n_dets = length(affected_dets)

        if (dem.L[i]==0)
            val0,val1 = (1. - dem.p[i], dem.p[i])
        else
            # Hadamard trick
            val0,val1 = (1. - dem.p[i], -dem.p[i])
        end

        if n_dets == 1 # single-site hyperedge
            idx_middle = Index(2)
            idx_out = Index(2)
            pos = pos_checks[affected_dets[1]]
            T = ITensor([val0,val1], idx_middle)
            T *= ITensor(CHECK3, open_index[pos], idx_middle, idx_out)
            contract1!(tn3d, pos, T)
            open_index[pos] = idx_out
        else
            @assert n_dets in [2,3,4]
            # keep track of "connected" components inside affected detectors
            components = Vector(1:n_dets)
            # store the resulting bonds
            bonds = Set{Set{Pos3D}}() # set of bonds, where each bond is a set of 2 sites
            # check out which detectors are neighbors
            for i = 1:n_dets
                posi = pos_checks[affected_dets[i]]
                for j = i+1:n_dets
                    posj = pos_checks[affected_dets[j]]
                    if components[i]!=components[j] && sum(abs.(posi .- posj)) < 1.5
                        # i and j are in the same connected component
                        replace!(components, components[j]=>components[i])
                        push!(bonds, Set([posi,posj]))
                    end
                end
            end

            function bonds_between(pos1::Pos3D, pos2::Pos3D)
                dx = sign(pos2[1]-pos1[1])
                dy = sign(pos2[2]-pos1[2])
                dz = sign(pos2[3]-pos1[3])

                path_bonds = Vector{Set{Pos3D}}()
                curr = pos1
                while curr != pos2
                    if curr[1]!=pos2[1]
                        next = curr .+ (dx,0,0)
                    elseif curr[2]!=pos2[2]
                        next = curr .+ (0,dy,0)
                    else
                        next = curr .+ (0,0,dz)
                    end
                    push!(path_bonds, Set([curr,next]))
                    curr = next
                    if length(bonds) > 100
                        error("stuck in loop")
                    end
                end
                return path_bonds
            end
            function closest_distance(comp1::Int, comp2::Int)
                @assert comp1 != comp2
                mindist = 99999999
                mindist_bonds = Vector{Set{Pos3D}}()
                for i = 1:n_dets if components[i] == comp1
                    posi = pos_checks[affected_dets[i]]
                    for j = 1:n_dets if components[j] == comp2
                        posj = pos_checks[affected_dets[j]]
                        d = sum(abs.(posi .- posj))
                        if d < mindist
                            mindist = d
                            mindist_bonds = bonds_between(posi, posj)
                        end
                    end end
                end end
                return mindist, mindist_bonds
            end
            # while there are more than one connected component, find the smallest distances
            # between the component, then merge the two closest ones
            while length(unique(components)) > 1
                unique_comps = unique(components)
                combinations = [(c1,c2) for c1 in unique_comps for c2 in unique_comps if c2!=c1]
                dist = [closest_distance(c1,c2) for (c1,c2) in combinations]
                idx = argmin([val[1] for val in dist])
                c1,c2 = combinations[idx]
                replace!(components, c2=>c1)
                bonds = union(bonds, Set(dist[idx][2]))
            end

            # figure out which sites are used for "bridging"
            main_sites = [pos_checks[j] for j in affected_dets]
            bridging_sites = setdiff(union(bonds...), main_sites)

            # open up sites
            temp_open_index = Dict{Pos3D,Index}()
            for pos in main_sites 
                idx_new = Index(2)
                idx_middle = Index(2)
                idx_temp = Index(2)
                T = ITensor(CHECK3, open_index[pos], idx_middle, idx_new)
                arr = (pos == main_sites[1]) ? [val0 0. ; 0. val1] : [1. 0. ; 0. 1.]
                T *= ITensor(arr, idx_middle, idx_temp)
                contract1!(tn3d, pos, T)
                open_index[pos] = idx_new
                temp_open_index[pos] = idx_temp
            end
            for pos in bridging_sites 
                idx_temp = Index(2)
                T = ITensor([1., 1.], idx_temp)
                contract1!(tn3d, pos, T)
                temp_open_index[pos] = idx_temp
            end

            # apply 2-qubit gates
            for b in bonds
                bvec = collect(b)
                @assert length(bvec) == 2
                pos1,pos2 = (bvec[1], bvec[2])
                idx1 = Index(2)
                idx2 = Index(2)
                idxmid = Index(2)
                T =  ITensor(DELTA3, temp_open_index[pos1], idxmid, idx1)
                T *= ITensor(DELTA3, temp_open_index[pos2], idxmid, idx2)
                expval += contract2_SU!(tn3d, pos1, pos2, T, Vector{Index}([idx1]), max_bond_dim, svd_eps)
                expval += renormalize_tensor!(tn3d, pos1)
                expval += renormalize_tensor!(tn3d, pos2)
                expval += renormalize_bond_matrix!(tn3d, pos1, pos2)
                temp_open_index[pos1] = idx1
                temp_open_index[pos2] = idx2
            end

            # close sites
            for pos in main_sites 
                T = ITensor([1.,1.], temp_open_index[pos])
                contract1!(tn3d, pos, T)
            end
            for pos in bridging_sites
                T = ITensor([1.,1.], temp_open_index[pos])
                contract1!(tn3d, pos, T)
            end

            # gauging step
            if do_gauging
                while canonicalness(tn3d) > canonicalness_target 
                    for x=xmin:xmax for y=ymin:ymax for z=zmin:zmax
                        for pos2 in neighbors(tn3d,(x,y,z)) if pos2>(x,y,z)
                            if dim(commonind( tn3d[(x,y,z)], bond_matrix(tn3d, (x,y,z),pos2) )) > 1
                            expval += identity_SU!(tn3d, (x,y,z), pos2, max_bond_dim, svd_eps)
                            expval += renormalize_tensor!(tn3d, (x,y,z))
                            expval += renormalize_tensor!(tn3d, pos2)
                            expval += renormalize_bond_matrix!(tn3d, (x,y,z), pos2)
                            end
                        end end
                    end end end
                end
            end
        end
    end

    return tn3d, expval, open_index
end

function contract_cubictn(tn::ContractionNetwork, max_bond_dim::Int, sc_bond_dim::Int ; svd_eps::Float64=1e-16, do_gauging::Bool=false, canonicalness_threshold::Float64=1e-4, split_dim::Int=2, ignore_bond_tensors::Bool=false)
    xmin,xmax = extrema([k[1] for (k,v) in tn.tensors])
    ymin,ymax = extrema([k[2] for (k,v) in tn.tensors])
    zmin,zmax = extrema([k[3] for (k,v) in tn.tensors])

    # split up the bond tensors and put them onto sites
    if !ignore_bond_tensors
        for x=xmin:xmax for y=ymin:ymax for z=zmin:zmax
            pos1 = (x,y,z)
            for pos2 in neighbors(tn, pos1) if pos2 > pos1
                idx1 = commonind(tn[pos1], bond_matrix(tn, pos1, pos2))
                idx2 = commonind(tn[pos2], bond_matrix(tn, pos1, pos2))
                idxnew = Index(dim(idx1))
                arr = sqrt.(Matrix(bond_matrix(tn,pos1,pos2), idx1, idx2))
                tn[pos1] *= ITensor(arr, idx1, idxnew)
                tn[pos2] *= ITensor(arr, idx2, idxnew)
            end end
        end end end
    end

    # set up PEPS 
    peps = generate_empty_PEPS(xmin, xmax, ymin, ymax)
    open_edge = Dict{Tuple{Int,Int},Index}()
    for x=xmin:xmax for y=ymin:ymax
        open_edge[x,y] = open_indices(peps, (x,y))[1]
        tn[(x,y,zmin)] *= ITensor([1.], open_edge[x,y])
    end end

    # contract
    expval = 0
    for z=zmin:zmax
        # contract two site bonds on same level
        for x=xmin:xmax for y=ymin:ymax
            for (x2,y2) in neighbors(peps, (x,y)) if (x2,y2) > (x,y)
                downidx1 = commonind(peps[(x,y)], tn[(x,y,z)])
                downidx2 = commonind(peps[(x2,y2)], tn[(x2,y2,z)])
                link = commonind(tn[(x,y,z)], tn[(x2,y2,z)])

                U1,S1,V1 = svd(tn[(x,y,z)], [idx for idx in [downidx1, link] if !isnothing(idx)], maxdim=split_dim)
                U2,S2,V2 = svd(tn[(x2,y2,z)], [idx for idx in [downidx2, link] if !isnothing(idx)], maxdim=split_dim)
                #U1,S1,V1 = svd(tn[(x,y,z)], [idx for idx in [downidx1, link] if !isnothing(idx)], maxdim=max_bond_dim)
                #U2,S2,V2 = svd(tn[(x2,y2,z)], [idx for idx in [downidx2, link] if !isnothing(idx)], maxdim=max_bond_dim)
                upidx1 = commonind(U1,S1)
                upidx2 = commonind(U2,S2)
                tn[(x,y,z)] = S1*V1
                tn[(x2,y2,z)] = S2*V2

                expval += renormalize_tensor!(tn, (x,y,z))
                expval += renormalize_tensor!(tn, (x2,y2,z))
                h = Int(floor(log2(LinearAlgebra.norm( U1 ))))
                U1 /= exp2(h)
                expval += h
                h = Int(floor(log2(LinearAlgebra.norm( U2 ))))
                U2 /= exp2(h)
                expval += h
                expval += contract2_SU!(peps, (x,y), (x2,y2), U1*U2, Vector{Index}([upidx1]), max_bond_dim, svd_eps)
                expval += renormalize_tensor!(peps, (x,y))
                expval += renormalize_tensor!(peps, (x2,y2))
                expval += renormalize_bond_matrix!(peps, (x,y), (x2,y2))
                open_edge[x,y] = upidx1
                open_edge[x2,y2] = upidx2 
            end end
        end end

        # contract the final remaining R into the peps
        for x=xmin:xmax for y=ymin:ymax
            peps[(x,y)] *= tn[(x,y,z)]
            expval += renormalize_tensor!(peps, (x,y))
            if z < zmax
                @assert length(open_indices(peps, (x,y))) == 1
                open_edge[x,y] = open_indices(peps, (x,y))[1]
            end
        end end

        # gauge
        if z < zmax && do_gauging
            n_steps = 0
            while canonicalness(peps) > canonicalness_threshold && n_steps < 100
                for x=xmin:xmax for y=ymin:ymax
                    for (x2,y2) in neighbors(peps, (x,y)) if (x2,y2) > (x,y)
                        expval += identity_SU!(peps, (x,y), (x2,y2), max_bond_dim, svd_eps)
                        expval += renormalize_tensor!(peps, (x,y))
                        expval += renormalize_tensor!(peps, (x2,y2))
                        expval += renormalize_bond_matrix!(peps, (x,y), (x2,y2))
                    end end
                end end
            end
        end
    end

    # contract final PEPS with sweep contractor
    res = contract_PEPS(peps, sc_bond_dim, svd_eps)
    return (res[1], res[2]+expval)
end
