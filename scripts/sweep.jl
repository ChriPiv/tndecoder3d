using MKL
using ArgParse
using JSON
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
    "--bond_dimension", "-b"
        arg_type = Int
        help = "3D bond dimensions"
        nargs = '+'
        required = true
    "--out_file", "-o"
        arg_type = String
        help = "file to which output data is written to"
        required = true
    "--mps_bond_dimension"
        arg_type = Int
        help = "2D bond dimensions"
        default = 8
     "--sector"
        arg_type = String
        help = "sector to consider, may be either 'X' or 'Z'. These correspond to the 23% and 3% sector respectively."
        default = "Z"
    "--output_write_interval"
        arg_type = Int
        help = "after how many samples should the output file be updated"
        default = 20
    "--num_samples", "-n"
        arg_type = Int
        help = "number of Monte Carlo samples"
        default = 1_000_000
    "--svd_eps"
        arg_type = Int
        help = "cuttof in power of 10 passed to ITensors's svd. Default is 16 (corresponding to 1e-16)"
        default = 16
end
args = parse_args(s)
@assert args["sector"] in ["X", "Z"]

data = Dict{Tuple{Int,Int,Int}, Tuple{Int,Int}}() # (distance,err_rate,bond_dim) -> (num_succ,num_fail)
for d in args["distance"] for p_int in args["error_rate"] for bd in args["bond_dimension"]
    data[d,p_int,bd] = (0,0)
end end end

svd_eps = 10. ^ (-args["svd_eps"]) 

for d in args["distance"]
    sc3d = unrotated_3D_surface_code(d)
    for p_int in args["error_rate"]
        p = p_int / 10000.
        for i = 1:args["num_samples"]
            err = convert(Vector{UInt8}, rand(sc3d.n_qubits) .< p)
            for bd in args["bond_dimension"]
                if args["sector"] == "Z"
                    syndrome = convert(Vector{UInt8}, (sc3d.Hz * err) .% 2)
                    err_L = sum(err .& sc3d.Lz) % 2
                    tn0,tn1 = gen_tn_sc3d_phenomenological(sc3d, d, p, syndrome)
                    GC.gc()
                    res = contract(tn1, bd, args["mps_bond_dimension"], svd_eps)
                    GC.gc()
                    dec_L = (res[1]>0.) ? 0 : 1
                else
                    syndrome = convert(Vector{UInt8}, (sc3d.Hx * err) .% 2)
                    repr0, repr1 = find_representatives(sc3d.Hx, sc3d.Lx, syndrome)
                    err_L = sum(err .& sc3d.Lx) % 2
                    tn0,tn1 = gen_tn_sc3d_dual(sc3d, d, p, repr0, repr1)
                    GC.gc()
                    res0 = contract(tn0, bd, args["mps_bond_dimension"], svd_eps)
                    GC.gc()
                    res1 = contract(tn1, bd, args["mps_bond_dimension"], svd_eps)
                    GC.gc()
                    dec_L = get_res(res0,res1)
                end

                if dec_L == err_L
                    data[d,p_int,bd] = data[d,p_int,bd] .+ (1,0)
                else
                    data[d,p_int,bd] = data[d,p_int,bd] .+ (0,1)
                end 
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
