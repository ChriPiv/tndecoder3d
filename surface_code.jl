function stab_3DSC_Z(d::Int)
    ind=Dict{Tuple{Int,Int,Int},Int}()
    for z=1:2d-1, y=0:2d-2, x=0:2d-2
        if isodd(x)+isodd(y)+isodd(z)==1
            ind[x,y,z]=length(ind)+1
        end
    end
    n=d*(3d^2-4d+2)
    ns=d*(3d^2-4d+2)
    is=1
    S=falses(n,ns)
    stab_pos = Vector{Tuple{Float64,Float64,Float64}}()
    for z=1:2d-1, y=0:2d-2, x=0:2d-2
        if isodd(x)+isodd(y)+isodd(z)==2
            if iseven(x)
                S[ind[x,y-1,z],is]=true
                S[ind[x,y+1,z],is]=true
                (z>1) && (S[ind[x,y,z-1],is]=true)
                (z<2d-1) && (S[ind[x,y,z+1],is]=true)
                push!(stab_pos, (x,y,z))
            elseif iseven(y)
                if y==0
                    S[ind[x-1,y,z],is]=true
                    S[ind[x+1,y,z],is]=true
                    (z>1) && (S[ind[x,y,z-1],is]=true)
                    (z<2d-1) && (S[ind[x,y,z+1],is]=true)
                    push!(stab_pos, (x,y,z))
                else
                    is-=1
                end
            else
                S[ind[x-1,y,z],is]=true
                S[ind[x+1,y,z],is]=true
                S[ind[x,y-1,z],is]=true
                S[ind[x,y+1,z],is]=true
                push!(stab_pos, (x,y,z))
            end
            is+=1
        end
    end

    qubit_pos = Vector{Tuple{Float64,Float64,Float64}}(undef, n)
    for (k,v) in ind
        qubit_pos[v] = convert(Tuple{Float64,Float64,Float64}, k)
    end

    return S[:,1:is-1], stab_pos, qubit_pos
end

function log_3DSC_Z(d::Int)::BitVector
    L=falses(d*(3d^2-4d+2))
    L[1:d*(3d-2):end].=true
    return L
end

function stab_3DSC_X(d::Int)
    ind=Dict{Tuple{Int,Int,Int},Int}()
    for z=1:2d-1, y=0:2d-2, x=0:2d-2
        if isodd(x)+isodd(y)+isodd(z)==1
            ind[x,y,z]=length(ind)+1
        end
    end
    n=d*(3d^2-4d+2)
    ns=(d-1)*d^2
    is=1
    S=falses(n,ns)
    stab_pos = Vector{Tuple{Float64,Float64,Float64}}(undef, ns)
    for z=1:2d-1, y=0:2d-2, x=0:2d-2
        if iseven(x)&&iseven(y)&&iseven(z)
            (x>0) && (S[ind[x-1,y,z],is]=true)
            (x<2d-2) && (S[ind[x+1,y,z],is]=true)
            (y>0) && (S[ind[x,y-1,z],is]=true)
            (y<2d-2) && (S[ind[x,y+1,z],is]=true)
            S[ind[x,y,z-1],is]=true
            S[ind[x,y,z+1],is]=true
            stab_pos[is] = (x,y,z)
            is+=1
        end
    end

    qubit_pos = Vector{Tuple{Float64,Float64,Float64}}(undef, n)
    for (k,v) in ind
        qubit_pos[v] = convert(Tuple{Float64,Float64,Float64}, k)
    end

    return S, stab_pos, qubit_pos
end

function log_3DSC_X(d::Int)::BitVector
    L=falses(d*(3d^2-4d+2))
    # L[d^2+1:d*(3d-2)].=true
    L[1:d^2].=true
    return L
end


function unrotated_3D_surface_code(d::Int)
    # we switch X and Z to match panqec's convention
    Hx, pos_stab_x, pos_qubits = stab_3DSC_Z(d)
    Hx = convert(Matrix{UInt8}, transpose(Hx))
    Lx = convert(Vector{UInt8}, log_3DSC_Z(d))
    Hz, pos_stab_z, pos_qubits = stab_3DSC_X(d)
    Hz = convert(Matrix{UInt8}, transpose(Hz))
    Lz = convert(Vector{UInt8}, log_3DSC_X(d))
    n_qubits = size(Hz, 2)
    return CSSCode{3}(n_qubits, Hx, Hz, Lx, Lz, pos_qubits, pos_stab_x, pos_stab_z)
end

