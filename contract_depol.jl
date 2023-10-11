function gen_tn_general(sc3d::CSSCode{3}, d::Int, PTensor::Matrix{Float64}, syndrome_x::Vector{UInt8}, syndrome_z::Vector{UInt8})
    pos_qubits = [convert(Pos3D, pos) for pos in sc3d.pos_qubits]
    pos_stab_x = [convert(Pos3D, pos) for pos in sc3d.pos_stab_x]
    pos_stab_z = [convert(Pos3D, pos) for pos in sc3d.pos_stab_z]
    xmin,xmax = extrema([pos[1] for pos in [pos_qubits ; pos_stab_x ; pos_stab_z]])
    ymin,ymax = extrema([pos[2] for pos in [pos_qubits ; pos_stab_x ; pos_stab_z]])
    zmin,zmax = extrema([pos[3] for pos in [pos_qubits ; pos_stab_x ; pos_stab_z]])
    # build empty 3D tensor network
    tnI = generate_empty_cubictn(xmin, xmax, ymin, ymax, zmin, zmax)
    tnX = generate_empty_cubictn(xmin, xmax, ymin, ymax, zmin, zmax)
    tnY = generate_empty_cubictn(xmin, xmax, ymin, ymax, zmin, zmax)
    tnZ = generate_empty_cubictn(xmin, xmax, ymin, ymax, zmin, zmax)
    for x=xmin:xmax, y=ymin:ymax, z=zmin:zmax
        # remove bond tensor connections
        # (it's kinda hacky that I use the ContractionNetwork object here)
        tnI[(x,y,z)] = ITensor([1.])
        tnX[(x,y,z)] = ITensor([1.])
        tnY[(x,y,z)] = ITensor([1.])
        tnZ[(x,y,z)] = ITensor([1.])
    end

    # build a map of indices
    indices = Dict{Set{Pos3D},Index}()
    for i=1:sc3d.n_qubits, j=1:size(sc3d.Hx,1)
        if sc3d.Hx[j,i] == 1
            @assert pos_qubits[i] in neighbors(tnX, pos_stab_x[j])
            indices[Set([pos_qubits[i], pos_stab_x[j]])] = Index(2)
        end
    end
    for i=1:sc3d.n_qubits, j=1:size(sc3d.Hz,1)
        if sc3d.Hz[j,i] == 1
            @assert pos_qubits[i] in neighbors(tnX, pos_stab_z[j])
            indices[Set([pos_qubits[i], pos_stab_z[j]])] = Index(2)
        end
    end
    
    # add check tensors
    for j=1:size(sc3d.Hx,1)
        pos = pos_stab_x[j]
        inds = Vector{Index}()
        for pos2 in neighbors(tnX, pos)
            if Set([pos,pos2]) in keys(indices)
                push!(inds, indices[Set([pos,pos2])])
            end
        end
        if syndrome_x[j] == 0
            arr = check_arr(length(inds), 1., 0.)
        else
            arr = check_arr(length(inds), 0., 1.)
        end
        tnI[pos] = ITensor(arr, inds...)
        tnX[pos] = ITensor(arr, inds...)
        tnY[pos] = ITensor(arr, inds...)
        tnZ[pos] = ITensor(arr, inds...)
    end
    for j=1:size(sc3d.Hz,1)
        pos = pos_stab_z[j]
        inds = Vector{Index}()
        for pos2 in neighbors(tnX, pos)
            if Set([pos,pos2]) in keys(indices)
                push!(inds, indices[Set([pos,pos2])])
            end
        end
        if syndrome_z[j] == 0
            arr = check_arr(length(inds), 1., 0.)
        else
            arr = check_arr(length(inds), 0., 1.)
        end
        tnI[pos] = ITensor(arr, inds...)
        tnX[pos] = ITensor(arr, inds...)
        tnY[pos] = ITensor(arr, inds...)
        tnZ[pos] = ITensor(arr, inds...)
    end

    # add depol tensors
    for i=1:sc3d.n_qubits
        pos = pos_qubits[i]

        # find indices of connected X and Z stabilizers
        inds_X = Vector{Index}()
        inds_Z = Vector{Index}()
        for pos2 in neighbors(tnX, pos)
            if Set([pos,pos2]) in keys(indices)
                if pos2 in pos_stab_x
                    push!(inds_X, indices[Set([pos,pos2])])
                else
                    push!(inds_Z, indices[Set([pos,pos2])])
                end
            end
        end

        # internal indices
        idx_X = Index(2)
        idx_Z = Index(2)
        
        T = ITensor(PTensor, idx_X, idx_Z)
        T *= ITensor(delta_arr(length(inds_Z)+1, 1., 1.), inds_Z..., idx_Z)
        T *= ITensor(delta_arr(length(inds_X)+1, 1., 1.), inds_X..., idx_X)
        tnI[pos] = T 

        T = ITensor(PTensor, idx_X, idx_Z)
        val = (sc3d.Lx[i]==1) ? (1., -1.) : (1., 1.)
        T *= ITensor(delta_arr(length(inds_X)+1, val[1], val[2]), inds_X..., idx_X)
        T *= ITensor(delta_arr(length(inds_Z)+1, 1., 1.), inds_Z..., idx_Z)
        tnX[pos] = T 

        T = ITensor(PTensor, idx_X, idx_Z)
        val = (sc3d.Lx[i]==1) ? (1., -1.) : (1., 1.)
        T *= ITensor(delta_arr(length(inds_X)+1, val[1], val[2]), inds_X..., idx_X)
        val = (sc3d.Lz[i]==1) ? (1., -1.) : (1., 1.)
        T *= ITensor(delta_arr(length(inds_Z)+1, val[1], val[2]), inds_Z..., idx_Z)
        tnY[pos] = T 
         
        T = ITensor(PTensor, idx_X, idx_Z)
        T *= ITensor(delta_arr(length(inds_X)+1, 1., 1.), inds_X..., idx_X)
        val = (sc3d.Lz[i]==1) ? (1., -1.) : (1., 1.)
        T *= ITensor(delta_arr(length(inds_Z)+1, val[1], val[2]), inds_Z..., idx_Z)
        tnZ[pos] = T 
    end

    # add dummy indices to make sure it's truly cubic
    for x=xmin:xmax, y=ymin:ymax, z=zmin:zmax
        for pos2 in neighbors(tnX, (x,y,z))
            if length(commoninds(tnX[(x,y,z)], tnX[pos2])) == 0
                idx = Index(1)
                tnI[(x,y,z)] *= ITensor([1.], idx)
                tnI[pos2] *= ITensor([1.], idx)
                tnX[(x,y,z)] *= ITensor([1.], idx)
                tnX[pos2] *= ITensor([1.], idx)
                tnY[(x,y,z)] *= ITensor([1.], idx)
                tnY[pos2] *= ITensor([1.], idx)
                tnZ[(x,y,z)] *= ITensor([1.], idx)
                tnZ[pos2] *= ITensor([1.], idx)
            end
        end
    end

    return tnI, tnX, tnY, tnZ 
end

function gen_tn_depol(sc3d::CSSCode{3}, d::Int, err_rate::Float64, syndrome_x::Vector{UInt8}, syndrome_z::Vector{UInt8})
    return gen_tn_general(sc3d, d, [(1. - err_rate) err_rate/3. ; err_rate/3. err_rate/3.], syndrome_x, syndrome_z)
end
