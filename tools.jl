DELTA3 = zeros(Float64, 2, 2, 2)
DELTA3[begin] = 1.
DELTA3[end] = 1.
CHECK3 = zeros(Float64, 2, 2, 2)
CHECK3[1,1,1] = 1.
CHECK3[1,2,2] = 1.
CHECK3[2,2,1] = 1.
CHECK3[2,1,2] = 1.

function delta_arr(num_legs::Int64, val1::Float64, val2::Float64)
    t = zeros(2*ones(Int,num_legs)...)
    t[begin] = val1
    t[end] = val2
    return t
end
function check_arr(num_legs::Int64, val1::Float64, val2::Float64)
    t = zeros(2*ones(Int,num_legs)...)
    for i = 1:2^num_legs
        t[i] = (count_ones(i-1)%2 == 0) ? val1 : val2
    end 
    return t
end

function get_res(res0::Tuple{Float64,Int}, res1::Tuple{Float64,Int})
    if res0[1] == 0.
        return UInt8(1)
    elseif res1[1] == 0.
        return UInt8(0)
    end

    if res0[2] == res1[2]
        return UInt8(abs(res1[1]) > abs(res0[1]))
    else
        return UInt8(res1[2] > res0[2])
    end
end


function permute_TN(tn::BitFlipTN, perm::Vector{Int})
    @assert length(tn) == length(perm)
    tn2 = deepcopy(tn)
    iperm = invperm(perm)
    permute!(tn2, perm)
    for i=1:length(tn2)
        tn2[i].adj = [iperm[j] for j in tn2[i].adj]
    end
    return tn2
end


function rank_mod2(A::Matrix{UInt8})
    m = size(A,1)
	n = size(A,2)
    # WLOG assume we have more cols than rows 
    tranposed = false
    if n < m 
        transposed = true
        A = transpose(A)
        n, m = m, n
    end
	
    basis_size = 0 
    M = zeros(UInt8, (m, m)) # reduced A in lower triangular form
    free_indices= zeros(Int32, m) # indices of non-fixed bits that form basis
    for i = 1:n
        v = copy(A[:,i])
        for j = 1:m
            if v[j] == 1
                if free_indices[j] != 0
                    v = (v .+ M[:,j]) .% 2
                else
                    # add to M
                    M[:,j] = v
                    free_indices[j] = i
                    basis_size += 1
                    break
                end
            end
        end
        if basis_size == m
            break
        end
    end

    return basis_size
end

function linearly_independent_rows(A::Matrix{UInt8})
    num_rows = size(A)[1]
    B = zeros(UInt8, size(A)...)

    idx = 1
    retval = Vector{Int}()
    for i = 1:num_rows
        B[idx,:] = A[i,:]
        if rank_mod2(B[1:idx,:]) == idx
            # row was linearly independent
            idx += 1
            push!(retval, i)
        end
    end
    return retval
end

function _solve_mod2_lineq!(A::Matrix{UInt8}, b::Vector{UInt8})
    # we assume that A is full-rank here.
    # A and b get modified in-place, solution is stored in b

    # first we bring A to upper triangular form
    n = size(A)[1]
    for i = 1:n # iterate through columns
        if A[i,i] == 0
            # need to swap row with another one
            found_swapping_row = false
            for j = i+1:n
                if A[j,i] == 1
                    A[i,:], A[j,:] = A[j,:], A[i,:]
                    b[i], b[j] = b[j], b[i]
                    found_swapping_row = true
                    break
                end
            end
            if not found_swapping_row
                error("Did not find row to swap with")
            end
        end
        # make sure all entries below (i,i) are set to zero
        for j = i+1:n
            if A[j,i] == 1
                A[j,:] = (A[j,:] + A[i,:] ) .% 2
                b[j] = (b[j] + b[i] ) .% 2
            end
        end
    end 

    # now make A diagonal
    for i = 1:n # iterate through columns
        for j = 1:i-1 # iterate through rows above i
            if A[j,i] == 1
                A[j,:] = (A[j,:] + A[i,:] ) .% 2
                b[j] = (b[j] + b[i] ) .% 2
            end
        end
    end 
end



function find_representatives(H::Matrix{UInt8}, l::Vector{UInt8}, syndromes::Vector{UInt8})
	H_stacked = [H ; transpose(l)]
	m = size(H_stacked,1)
	n = size(H_stacked,2)

    # goal: find linearly first k linaerly independent columns of H_stacked
    basis_size = 0 
    M = zeros(UInt8, (m, m)) # reduced H in lower triangular form
    free_indices= zeros(Int32, m) # indices of non-fixed bits that form basis
    for i = 1:n
        v = copy(H_stacked[:,i])
        for j = 1:m
            if v[j] == 1
                if free_indices[j] != 0
                    v = (v .+ M[:,j]) .% 2
                else
                    # add to M
                    M[:,j] = v
                    free_indices[j] = i
                    basis_size += 1
                    break
                end
            end
        end
        if basis_size == m
            break
        end
    end
    if basis_size != m
        error("Couldn't find full linearly independent basis set")
    end

	H_reduced = H_stacked[:,free_indices]
	y0 = vcat(syndromes, Vector{UInt8}([0]))
	y1 = vcat(syndromes, Vector{UInt8}([1]))
	_solve_mod2_lineq!(copy(H_reduced[:,1:m]), y0)
	_solve_mod2_lineq!(copy(H_reduced[:,1:m]), y1)
	x0 = zeros(UInt8, n)
	x1 = zeros(UInt8, n)
	x0[free_indices] = y0
	x1[free_indices] = y1
	return x0, x1

end

function gen_mle_table(H::Matrix{UInt8}, l::Vector{UInt8}, p::Vector{Float64})
	m = size(H,1)
	n = size(H,2)
    
	table_size = 2^(m + 1)
    mle_table = zeros(table_size)
    temp_buffer = zeros(table_size)

    mle_table[1] = 1 
    for i = 0:n-1
        # iterate through each fault
        idx = sum(H[:,i+1] .* (2 .^ (m-1:-1:0))) # error events triggered by fault
        idx = 2*idx + l[i+1] # add logical information to index
        # update rule
        for j = 0:table_size-1
            temp_buffer[j+1] = mle_table[j+1] * (1-p[i+1]) + mle_table[ (j⊻idx) + 1] * p[i+1]
        end
        for j = 0:table_size-1
            mle_table[j+1] = temp_buffer[j+1]
        end
    end 
    return mle_table
end
function gen_mle_table(H::Matrix{UInt8}, l::Vector{UInt8}, p::Float64)
    return gen_mle_table(H, l, p*ones(Float64, length(l)))
end

function mle_table_lookup(H::Matrix{UInt8}, mle_table::Vector{Float64}, mmt_vector::Vector{UInt8})
    m= size(H,1)
	n = size(H,2)

	data_as_num = sum( mmt_vector .* 2 .^ (m-1:-1:0) )
    idx0 = data_as_num * 2 
    idx1 = data_as_num * 2 + 1 
    p0 = mle_table[idx0+1]
    p1 = mle_table[idx1+1]

    return p0, p1
end

