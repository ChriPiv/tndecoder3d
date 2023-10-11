const Pos3D = Tuple{Int,Int,Int}

struct CubicTNTensorDescription
    type::TensorType
    val::Tuple{Float64,Float64}
end
struct CubicTN
    xmin::Int
    xmax::Int
    ymin::Int
    ymax::Int
    zmin::Int
    zmax::Int

    sites::Dict{Pos3D,CubicTNTensorDescription}
    one_site_bonds::Dict{Pos3D,CubicTNTensorDescription}
    two_site_bonds::Dict{Set{Pos3D},CubicTNTensorDescription}
end

function gen_tn_sc3d_phenomenological(sc3d::CSSCode{3}, d::Int, err_rate::Float64, syndrome::Vector{UInt8})
    # Z-type stabilizers are the star ones (i.e. weight 6)
    # therefore we want to do the DEM picture of the Z sector
    num_Z_stabs = size(sc3d.Hz, 1)

    # find the dimensions of our qubic lattice
    pos_stab_z = [convert(Tuple{Int,Int,Int}, round.(0.5.*pos)) for pos in sc3d.pos_stab_z]
    pos_qubits = [0.5.*pos for pos in sc3d.pos_qubits]
    xmin, xmax = extrema([pos[1] for pos in pos_stab_z])
    ymin, ymax = extrema([pos[2] for pos in pos_stab_z])
    zmin, zmax = extrema([pos[3] for pos in pos_stab_z])

    # iterate through all Z stabs
    sites = Dict{Pos3D,CubicTNTensorDescription}()
    for i = 1:num_Z_stabs
        val = (syndrome[i]==1) ? (0.,1.) : (1.,0.)
        sites[pos_stab_z[i]] = CubicTNTensorDescription(TT_CHECK, val)
    end

    # iterate through all qubits
    one_site_bonds0 = Dict{Pos3D,CubicTNTensorDescription}()
    one_site_bonds1 = Dict{Pos3D,CubicTNTensorDescription}()
    two_site_bonds0 = Dict{Set{Pos3D},CubicTNTensorDescription}()
    two_site_bonds1 = Dict{Set{Pos3D},CubicTNTensorDescription}()
    for i = 1:sc3d.n_qubits
        val0 = (1-err_rate, err_rate)
        val1 = (sc3d.Lz[i]==0) ? (1-err_rate, err_rate) : (1-err_rate, -err_rate)
        pos = pos_qubits[i]
        @assert sum(isinteger.(pos)) == 2
        nbors = [convert(Tuple{Int,Int,Int}, floor.(pos)), convert(Tuple{Int,Int,Int}, ceil.(pos))]
        nbors = [p for p in nbors if xmin<=p[1]<=xmax && ymin<=p[2]<=ymax && zmin<=p[3]<=zmax]
        if length(nbors) == 1
            @assert !(nbors[1] in keys(one_site_bonds0))
            one_site_bonds0[nbors[1]] = CubicTNTensorDescription(TT_DELTA, val0)
            one_site_bonds1[nbors[1]] = CubicTNTensorDescription(TT_DELTA, val1)
        else
            key = Set([nbors[1], nbors[2]])
            @assert !(key in keys(two_site_bonds0))
            two_site_bonds0[key] = CubicTNTensorDescription(TT_DELTA, val0)
            two_site_bonds1[key] = CubicTNTensorDescription(TT_DELTA, val1)
        end
    end

    ctn0 = CubicTN(xmin, xmax, ymin, ymax, zmin, zmax, sites, one_site_bonds0, two_site_bonds0)
    ctn1 = CubicTN(xmin, xmax, ymin, ymax, zmin, zmax, sites, one_site_bonds1, two_site_bonds1)
    return ctn0, ctn1
end


function gen_tn_sc3d_dual(sc3d::CSSCode{3}, d::Int, err_rate::Float64, representative0::Vector{UInt8}, representative1::Vector{UInt8})
    # Z-type stabilizers are the star ones (i.e. weight 6)
    # therefore we want to do the gauge picture of the X sector
    num_Z_stabs = size(sc3d.Hz, 1)

    # find the dimensions of our qubic lattice
    pos_stab_z = [convert(Tuple{Int,Int,Int}, round.(0.5.*pos)) for pos in sc3d.pos_stab_z]
    pos_qubits = [0.5.*pos for pos in sc3d.pos_qubits]
    xmin, xmax = extrema([pos[1] for pos in pos_stab_z])
    ymin, ymax = extrema([pos[2] for pos in pos_stab_z])
    zmin, zmax = extrema([pos[3] for pos in pos_stab_z])

    # iterate through all Z stabs
    sites = Dict{Pos3D,CubicTNTensorDescription}()
    for i = 1:num_Z_stabs
        sites[pos_stab_z[i]] = CubicTNTensorDescription(TT_DELTA, (1.,1.))
    end

    # iterate through all qubits
    one_site_bonds0 = Dict{Pos3D,CubicTNTensorDescription}()
    one_site_bonds1 = Dict{Pos3D,CubicTNTensorDescription}()
    two_site_bonds0 = Dict{Set{Pos3D},CubicTNTensorDescription}()
    two_site_bonds1 = Dict{Set{Pos3D},CubicTNTensorDescription}()
    for i = 1:sc3d.n_qubits
        val0 = (representative0[i]==0) ? (1-err_rate, err_rate) : (err_rate, 1-err_rate)
        val1 = (representative1[i]==0) ? (1-err_rate, err_rate) : (err_rate, 1-err_rate)
        pos = pos_qubits[i]
        @assert sum(isinteger.(pos)) == 2
        nbors = [convert(Tuple{Int,Int,Int}, floor.(pos)), convert(Tuple{Int,Int,Int}, ceil.(pos))]
        nbors = [p for p in nbors if xmin<=p[1]<=xmax && ymin<=p[2]<=ymax && zmin<=p[3]<=zmax]
        if length(nbors) == 1
            @assert !(nbors[1] in keys(one_site_bonds0))
            one_site_bonds0[nbors[1]] = CubicTNTensorDescription(TT_CHECK, val0)
            one_site_bonds1[nbors[1]] = CubicTNTensorDescription(TT_CHECK, val1)
        else
            key = Set([nbors[1], nbors[2]])
            @assert !(key in keys(two_site_bonds0))
            two_site_bonds0[key] = CubicTNTensorDescription(TT_CHECK, val0)
            two_site_bonds1[key] = CubicTNTensorDescription(TT_CHECK, val1)
        end
    end

    ctn0 = CubicTN(xmin, xmax, ymin, ymax, zmin, zmax, sites, one_site_bonds0, two_site_bonds0)
    ctn1 = CubicTN(xmin, xmax, ymin, ymax, zmin, zmax, sites, one_site_bonds1, two_site_bonds1)
    return ctn0, ctn1
end

function contract_exact(ctn::CubicTN)
    xmin = ctn.xmin
    xmax = ctn.xmax
    ymin = ctn.ymin
    ymax = ctn.ymax
    zmin = ctn.zmin
    zmax = ctn.zmax
    neighbors(x,y,z) = [(x2,y2,z2) for (x2,y2,z2) in [(x+1,y,z),(x-1,y,z),(x,y+1,z),(x,y-1,z),(x,y,z+1),(x,y,z-1)] if (xmin<=x2<=xmax && ymin<=y2<=ymax && zmin<=z2<=zmax)]

    # generate indices
    indices = Dict{Set{Pos3D},Index}()
    for x=xmin:xmax for y=ymin:ymax for z=zmin:zmax
        for pos2 in neighbors(x,y,z) if pos2 > (x,y,z)
            indices[Set([(x,y,z), pos2])] = Index(2)
        end end
    end end end

    # build 3D TN
    tn3d = Dict{Pos3D,ITensor}()
    for x=xmin:xmax for y=ymin:ymax for z=zmin:zmax
        # create site
        t = ctn.sites[x,y,z]
        i = [indices[Set([(x,y,z),pos2])] for pos2 in neighbors(x,y,z)]
        if (x,y,z) in keys(ctn.one_site_bonds)
            osb_idx = Index(2)
            push!(i, osb_idx)
        end
        n_nbors = length(i)
        arr = (t.type == TT_DELTA) ? delta_arr(n_nbors, t.val[1], t.val[2]) : check_arr(n_nbors, t.val[1], t.val[2])
        T = ITensor(arr, i...)

        # contract one-site bond
        if (x,y,z) in keys(ctn.one_site_bonds)
            tosb = ctn.one_site_bonds[x,y,z]
            T *= ITensor([tosb.val[1], tosb.val[2]], osb_idx)
        end

        # contract two-site bonds
        i = Index(2)
        for (x2,y2,z2) in neighbors(x,y,z) if (x2,y2,z2) > (x,y,z)
            tbond = ctn.two_site_bonds[Set([(x,y,z), (x2,y2,z2)])]
            arrbond = (tbond.type == TT_DELTA) ? delta_arr(2, tbond.val[1], tbond.val[2]) : check_arr(2, tbond.val[1], tbond.val[2])
            T *= ITensor(arrbond, i, indices[Set([(x,y,z),(x2,y2,z2)])])
            T *= delta(i, indices[Set([(x,y,z),(x2,y2,z2)])])
        end end

        tn3d[(x,y,z)] = T
    end end end

    # contract 3D TN
    T = ITensor([1.])
    expval = 0
    for x=xmin:xmax for y=ymin:ymax for z=zmin:zmax
        T *= tn3d[(x,y,z)]
        h = Int(floor(log2(LinearAlgebra.norm( T ))))
        T /= exp2(h)
        expval += h
    end end end

    return (storage(T)[1], expval)
end


function contract(ctn::CubicTN, max_bond_dim::Int, sc_bond_dim::Int, svd_eps::Float64=1e-16, do_gauging::Bool=false, canonicalness_threshold::Float64=1e-4 )
    xmin = ctn.xmin
    xmax = ctn.xmax
    ymin = ctn.ymin
    ymax = ctn.ymax
    zmin = ctn.zmin
    zmax = ctn.zmax

    # set up PEPS 
    peps = generate_empty_PEPS(xmin, xmax, ymin, ymax)
    open_edge = Dict{Tuple{Int,Int},Index}()
    for x=xmin:xmax for y=ymin:ymax
        open_edge[x,y] = open_indices(peps, (x,y))[1]
    end end

    # contract
    expval = 0
    for z=zmin:zmax
        # contract value tensor on each site
        for x=xmin:xmax for y=ymin:ymax
            t = ctn.sites[x,y,z]
            if z == zmin
                arr = [t.val[1], t.val[2]]
            else
                if t.type == TT_DELTA
                    arr = [t.val[1] 0 ; 0 t.val[2]]
                else
                    arr = [t.val[1] t.val[2] ; t.val[2] t.val[1]]
                end
            end
            idx = Index(2)
            contract1!(peps, (x,y), ITensor(arr, open_edge[x,y], idx))
            open_edge[x,y] = idx
        end end

        # contract one site bonds
        for x=xmin:xmax for y=ymin:ymax
            if (x,y,z) in keys(ctn.one_site_bonds)
                t = ctn.one_site_bonds[x,y,z]
                idx_out = Index(2)
                idx_bond = Index(2)
                arr = (ctn.sites[x,y,z].type == TT_DELTA) ? DELTA3 : CHECK3
                T = ITensor(arr, open_edge[x,y], idx_out, idx_bond)
                T *= ITensor([t.val[1], t.val[2]], idx_bond)
                contract1!(peps, (x,y), T)
                open_edge[x,y] = idx_out
            end
        end end

        # contract two site bonds on same level
        for x=xmin:xmax for y=ymin:ymax
            for (x2,y2) in neighbors(peps, (x,y)) if (x2,y2) > (x,y)
                t1 = ctn.sites[x,y,z]
                t2 = ctn.sites[x2,y2,z]
                tbond = ctn.two_site_bonds[Set([(x,y,z), (x2,y2,z)])]

                if tbond.type == TT_DELTA
                    arrbond = [tbond.val[1] 0 ; 0 tbond.val[2]]
                else
                    arrbond = [tbond.val[1] tbond.val[2] ; tbond.val[2] tbond.val[1]]
                end
                arr1 = (t1.type==TT_DELTA) ? DELTA3 : CHECK3
                arr2 = (t2.type==TT_DELTA) ? DELTA3 : CHECK3

                idx1_out = Index(2)
                idx2_out = Index(2)
                idx1_bond = Index(2)
                idx2_bond = Index(2)
                T  = ITensor(arr1, open_edge[x,y], idx1_out, idx1_bond)
                T *= ITensor(arrbond, idx1_bond, idx2_bond)
                T *= ITensor(arr2, open_edge[x2,y2], idx2_out, idx2_bond)
                expval += contract2_SU!(peps, (x,y), (x2,y2), T, Vector{Index}([idx1_out]), max_bond_dim, svd_eps)
                expval += renormalize_tensor!(peps, (x,y))
                expval += renormalize_tensor!(peps, (x2,y2))
                expval += renormalize_bond_matrix!(peps, (x,y), (x2,y2))
                open_edge[x,y] = idx1_out
                open_edge[x2,y2] = idx2_out

            #for i = 1:0#while canonicalness(peps) > canonicalness_threshold && n_steps < 100
            #    for x=xmin:xmax for y=ymin:ymax
            #        for (x2,y2) in neighbors(peps, (x,y)) if (x2,y2) > (x,y)
            #            expval += identity_SU!(peps, (x,y), (x2,y2), max_bond_dim, svd_eps)
            #            expval += renormalize_tensor!(peps, (x,y))
            #            expval += renormalize_tensor!(peps, (x2,y2))
            #            expval += renormalize_bond_matrix!(peps, (x,y), (x2,y2))
            #        end end
            #    end end
            #end
            end end
        end end

        # contract subsequent two-site bond
        for x=xmin:xmax for y=ymin:ymax
            if z == zmax
                t = ctn.sites[x,y,z]
                arr = (t.type == TT_DELTA) ? [1. 1.] : [1. 0.]
                T = ITensor(arr, open_edge[x,y])
                contract1!(peps, (x,y), T)
            else
                t = ctn.two_site_bonds[Set([(x,y,z), (x,y,z+1)])]
                if t.type == TT_DELTA
                    arr = [t.val[1] 0 ; 0 t.val[2]]
                else
                    arr = [t.val[1] t.val[2] ; t.val[2] t.val[1]]
                end
                idx = Index(2)
                T = ITensor(arr, open_edge[x,y], idx)
                contract1!(peps, (x,y), T)
                open_edge[x,y] = idx
            end
        end end
       
    end

    # contract final PEPS with sweep contractor
    res = contract_PEPS(peps, sc_bond_dim, svd_eps)
    return (res[1], res[2]+expval)
end


function contract_PEPS(peps::ContractionNetwork{P}, sc_bond_dim::Int, svd_eps::Float64=1e-16) where P
    xmin, xmax = extrema([pos[1] for (pos,v) in peps.tensors])
    ymin, ymax = extrema([pos[2] for (pos,v) in peps.tensors])

    # contract bond matrices 
    for x=xmin:xmax for y=ymin:ymax
        for (x2,y2) in neighbors(peps, (x,y)) if (x2,y2) > (x,y)
            peps.tensors[x,y] *= bond_matrix(peps, (x,y), (x2,y2))
        end end
    end end

    # we sweep from ymin to ymax
    expval = 0
    mps = MPS([peps.tensors[x,ymin] for x=xmin:xmax])
    for y = ymin+1:ymax
        # contract
        for x = xmin:xmax
            i = x - xmin + 1
            mps.data[i] *= peps.tensors[x,y]
        end
        # normalize
        for i = 1:length(mps.data)
            h = Int(floor(log2(LinearAlgebra.norm( mps.data[i] ))))
            mps.data[i] /= exp2(h)
            expval += h
        end
        # truncate
        truncate!(mps, maxdim=sc_bond_dim, cutoff=svd_eps)
        # normalize
        for i = 1:length(mps.data)
            h = Int(floor(log2(LinearAlgebra.norm( mps.data[i] ))))
            mps.data[i] /= exp2(h)
            expval += h
        end
    end

    # get output
    T = ITensor([1.])
    for T2 in mps.data
        T *= T2
    end
    val = storage(T)[1]
    h = Int(floor(log2(LinearAlgebra.norm( val ))))
    val /= exp2(h)
    expval += h

    return val, expval
end

