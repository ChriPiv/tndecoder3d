using MKL
using ArgParse
using JSON
using PyCall
@pyimport pymatching
include("../tndecoder3d.jl")
using .TNDecoder3D

s = ArgParseSettings()
@add_arg_table s begin
    "--distance", "-d"
        arg_type = Int
        help = "code distance"
        nargs = '+'
        required = true
    "--error_rate", "-p"
        arg_type = Int 
        help = "physical error rate in unis of 1e-4"
        nargs = '+'
        required = true
    "--out_file", "-o"
        arg_type = String
        help = "file to which output data is written to"
        required = true
     "--sector"
        arg_type = String
        help = "sector to consider, may be either 'X' or 'Z'. These correspond to the 23% and 3% sector respectively."
        default = "Z"
    "--output_write_interval"
        arg_type = Int
        help = "after how many samples should the output file be updated"
        default = 10000
    "--num_samples", "-n"
        arg_type = Int
        help = "number of Monte Carlo samples"
        default = 1_000_000
end
args = parse_args(s)
@assert args["sector"] in ["X", "Z"]

data = Dict{Tuple{Int,Int}, Tuple{Int,Int}}() # (distance,err_rate,bond_dim) -> (num_succ,num_fail)
for d in args["distance"] for p_int in args["error_rate"]
    data[d,p_int] = (0,0)
end end

for d in args["distance"]
    sc3d = unrotated_3D_surface_code(d)
    if args["sector"] == "Z"
        matching_decoder = pymatching.Matching(sc3d.Hz)
    else
        matching_decoder = pymatching.Matching(sc3d.Hx)
    end
    for p_int in args["error_rate"]
        p = p_int / 10000.
        for i = 1:args["num_samples"]
            err = convert(Vector{UInt8}, rand(sc3d.n_qubits) .< p)
            if args["sector"] == "Z"
                syndrome = convert(Vector{UInt8}, (sc3d.Hz * err) .% 2)
                err_L = sum(err .& sc3d.Lz) % 2
                prediction = matching_decoder.decode(syndrome)
                prediction_L = sum(prediction .& sc3d.Lz) % 2
            else
                syndrome = convert(Vector{UInt8}, (sc3d.Hx * err) .% 2)
                err_L = sum(err .& sc3d.Lx) % 2
                prediction = matching_decoder.decode(syndrome)
                prediction_L = sum(prediction .& sc3d.Lx) % 2
            end
            if prediction_L == err_L
                data[d,p_int] = data[d,p_int] .+ (1,0)
            else
                data[d,p_int] = data[d,p_int] .+ (0,1)
            end 

            if i % args["output_write_interval"] == 0
                json_str = JSON.json(Dict(
                    "args" => args,
                    "data" => data
                ))
                open(args["out_file"], "w") do f
                    write(f, json_str)
                end
            end

        end
    end
end
