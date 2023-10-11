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

data = Dict{Tuple{Int,Int,Int}, Tuple{Int,Int}}() # (distance,err_rate,bond_dim) -> (num_succ,num_fail)
for d in args["distance"] for p_int in args["error_rate"] for bd in args["bond_dimension"]
    data[d,p_int,bd] = (0,0)
end end end

svd_eps = 10. ^ (-args["svd_eps"]) 

function tuple_add(tlist)
    max_exp = maximum([t[2] for t in tlist])
    base = sum([ t[1] * 2. ^ (t[2]-max_exp) for t in tlist ])
    if base == 0.
        return (0., -10000000)
    else
        h = Int(floor(log2(abs(base))))
        return (base/exp2(h) , max_exp+h)
    end
end

for d in args["distance"]
    sc3d = unrotated_3D_surface_code(d)
    for p_int in args["error_rate"]
        p = p_int / 10000.
        for i = 1:args["num_samples"]
            err = convert(Vector{UInt8}, rand(sc3d.n_qubits) .< p)
            for bd in args["bond_dimension"]
                # to sample depol noise, sample a random number in [0,1]
                # range (0,p/3) : X error
                # range (p/3,2*p/3) : Y error
                # range (2p/3,p) : Z error
                # range (p,1) : no error
                err = rand(sc3d.n_qubits)
                err_x = err .< (p * 2/3)
                err_z = (err .> (p * 1/3)) .& (err .< p)
                syndrome_x = convert(Vector{UInt8}, (sc3d.Hx * err_z) .% 2)
                syndrome_z = convert(Vector{UInt8}, (sc3d.Hz * err_x) .% 2)
                err_Lx = sum(err_z .& sc3d.Lx) % 2 
                err_Lz = sum(err_x .& sc3d.Lz) % 2

                if err_Lx==0 && err_Lz==0 idx_ref = 1
                elseif err_Lx==1 && err_Lz==0 idx_ref = 2
                elseif err_Lx==1 && err_Lz==1 idx_ref = 3
                elseif err_Lx==0 && err_Lz==1 idx_ref = 4
                else error("???")
                end

                tnI, tnX, tnY, tnZ = gen_tn_depol(sc3d, d, p, syndrome_x, syndrome_z)

                resI = (0., -10000000)
                resX = (0., -10000000)
                resY = (0., -10000000)
                resZ = (0., -10000000)
                try
                    resI = contract_cubictn(deepcopy(tnI), bd, args["mps_bond_dimension"] ; svd_eps=svd_eps, split_dim=4, ignore_bond_tensors=true)
                catch e end
                 try
                    resX = contract_cubictn(deepcopy(tnX), bd, args["mps_bond_dimension"] ; svd_eps=svd_eps, split_dim=4, ignore_bond_tensors=true)
                catch e end
                try
                    resY = contract_cubictn(deepcopy(tnY), bd, args["mps_bond_dimension"] ; svd_eps=svd_eps, split_dim=4, ignore_bond_tensors=true)
                catch e end
                try
                    resZ = contract_cubictn(deepcopy(tnZ), bd, args["mps_bond_dimension"] ; svd_eps=svd_eps, split_dim=4, ignore_bond_tensors=true)
                catch e end

                A00 = tuple_add([ ( resI[1],resI[2]), ( resY[1],resY[2]), ( resZ[1],resZ[2]), ( resX[1],resX[2]) ])
                A01 = tuple_add([ ( resI[1],resI[2]), (-resY[1],resY[2]), ( resZ[1],resZ[2]), (-resX[1],resX[2]) ])
                A10 = tuple_add([ ( resI[1],resI[2]), ( resY[1],resY[2]), (-resZ[1],resZ[2]), (-resX[1],resX[2]) ])
                A11 = tuple_add([ ( resI[1],resI[2]), (-resY[1],resY[2]), (-resZ[1],resZ[2]), ( resX[1],resX[2]) ])
                idx_dec = argmax([ reverse(t) for t in [A00,A01,A10,A11] ])

                if idx_ref == idx_dec
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


