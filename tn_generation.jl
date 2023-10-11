struct CSSCode{N}
    n_qubits::Int
    Hx::Matrix{UInt8}
    Hz::Matrix{UInt8}
    Lx::Vector{UInt8}
    Lz::Vector{UInt8}
    pos_qubits::Vector{NTuple{N,Float64}}
    pos_stab_x::Vector{NTuple{N,Float64}}
    pos_stab_z::Vector{NTuple{N,Float64}}
end

struct StabilizerCode{N}
    n_qubits::Int
    stabs::Matrix{UInt8} # rows are of size 2n to store x and z part
    Lx::Vector{UInt8} # vector of size 2n to store x and z part
    Lz::Vector{UInt8} # vector of size 2n to store x and z part
    pos_qubits::Vector{NTuple{N,Float64}}
    pos_stabs::Vector{NTuple{N,Float64}}
end

struct DEM{N}
    n_bits::Int
    H::Matrix{UInt8}
    L::Vector{UInt8}
    p::Vector{Float64}
    pos_checks::Vector{NTuple{N,Float64}}
end

@enum TensorType TT_DELTA=1 TT_CHECK=2 DEPOL=3
mutable struct BitFlipTensor{N}
    adj::Vector{Int}
    type::TensorType
    val::Tuple{Float64,Float64}
    pos::NTuple{N, Float64}
end
const BitFlipTN{N} = Vector{BitFlipTensor{N}} where N

mutable struct DepolTensor{N}
    adj::Vector{Int}
    type::TensorType
    adj_is_x::Vector{Bool} # true=X and false=Z. Only relevant for depol tensors.
    val::Tuple{Float64,Float64} # (1-p) and p/3 for depol tensor
    pos::NTuple{N, Float64}
end
const DepolTN{N} = Vector{DepolTensor{N}} where N



function get_arr(tensor::BitFlipTensor)::Array{Float64}
    if tensor.type == STT_DELTA
        num_legs = length(tensor.adj)
        t = zeros(2*ones(Int, num_legs)...)
        t[1] = tensor.val[1]
        t[end] = tensor.val[2]
        return t
    elseif tensor.type == STT_CHECK
        num_legs = length(tensor.adj)
        t = zeros(2*ones(Int, num_legs)...)
        for i = 0:(2^num_legs)-1
            t[i+1] = (count_ones(i)%2 == 0) ? tensor.val[1] : tensor.val[2]
        end
        return t
    else
        error("Invalid tensor type")
    end
end



function css_to_stab(css::CSSCode{N})::StabilizerCode{N} where N
    return StabilizerCode{N}(
        css.n_qubits,
        [css.Hx ; zeros(UInt8, size(css.Hx,1), size(css.Hz,2)) ;; zeros(UInt8, size(css.Hz,1), size(css.Hx,2)) ; css.Hz],
        [css.Lx; zeros(UInt8, css.n_qubits)],
        [zeros(UInt8, css.n_qubits); css.Lz],
        css.pos_qubits,
        [css.pos_stab_x; css.pos_stab_z]
    )
end


function dem_from_stim(stim_model::PyObject)
    # possible instructions in a detector error model are: detector, logical_observable, error, shift_detectors, repeat
    # the flattened() method removes the last two
    model_string = stim_model.flattened().__str__()

    # count number of errors and detectors
    num_error = 0
    num_detector = 0
    for line in split(model_string, "\n")
        if split(line, "(")[1] == "error"
            num_error += 1
        elseif split(line, "(")[1] == "detector"
            num_detector += 1
        end
    end

    # create dem objects
    H = zeros(UInt8, num_detector, num_error)
    L = zeros(UInt8, num_error)
    p = zeros(Float64, num_error)
    pos_checks = Vector{NTuple{3,Float64}}(undef, num_detector)

    # populate objects
    err_idx = 0
    for line in split(model_string, "\n")
        instr_name = split(line, "(")[1]
        arguments = split( split(line, "(")[2], ")")[1]
        arguments = [strip(s) for s in split(arguments, ",")]
        targets = strip(split(line, ")")[2])
        targets = [strip(s) for s in split(targets, " ")]

        if instr_name == "error"
            err_idx += 1
            p[err_idx] = parse(Float64, arguments[1])
            for t in targets
                if t == "L0"
                    L[err_idx] = 1
                elseif t[1] == 'D'
                    det_idx = parse(Int, t[2:end]) + 1
                    H[det_idx, err_idx] = 1
                else
                    error("Invalid error target: ", t)
                end
            end
        elseif instr_name == "detector"
            posx = parse(Float64, arguments[1])
            posy = parse(Float64, arguments[2])
            posz = parse(Float64, arguments[3])
            idx = parse(Int, targets[1][2:end]) + 1
            pos_checks[idx] = (posx,posy,posz)
        else
            error("Unsupported instruction: ", instr_name)
        end
    end
    
    return DEM{3}(num_error,H, L, p, pos_checks) 
end


