using JLD2
using MKL
using ArgParse
using JSON
include("../tndecoder3d.jl")
using .TNDecoder3D
using PyCall
using ITensors
@pyimport stim

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
"--truncation_dimension", "-t"
    arg_type = Int
    help = "bond dimension of pre-compressed network"
    default = 4 
"--split_dimension", "-s"
    arg_type = Int
    help = "bond dimension used during SVD splitting process"
    default = 8 
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

svd_eps = 10. ^ (-args["svd_eps"]) 

data = Dict{Tuple{Int,Int,Int}, Tuple{Int,Int}}() # (distance,err_rate,bond_dim) -> (num_succ,num_fail)
for d in args["distance"] for p_int in args["error_rate"] for bd in args["bond_dimension"]
data[d,p_int,bd] = (0,0)
end end end

for d in args["distance"]
    for p_int in args["error_rate"]
        p = p_int / 10000.
        circuit = stim.Circuit.generated("surface_code:rotated_memory_x", distance=d, rounds=d, after_clifford_depolarization=p, before_measure_flip_probability=p, after_reset_flip_probability=p)
        model = circuit.detector_error_model(decompose_errors=false)
        dem = dem_from_stim(model)
        sampler = model.compile_sampler()

        #tn3d,expval,open_index = gen_tn_sc3d_circ(dem, 6)
        tn3d,expval,open_index = load_object("/cpivetea/circtn"*string(d)*".jld2")

        pos_checks = [convert(Tuple{Int,Int,Int}, round.((0.5*pos[1], 0.5*pos[2], pos[3]))) for pos in dem.pos_checks]
        xmin,xmax = extrema([pos[1] for pos in pos_checks])
        ymin,ymax = extrema([pos[2] for pos in pos_checks])
        zmin,zmax = extrema([pos[3] for pos in pos_checks])

        final_truncation_dim = args["truncation_dimension"]
        if final_truncation_dim != 0
            for x=xmin:xmax for y=ymin:ymax for z=zmin:zmax
                for pos2 in TNDecoder3D.neighbors(tn3d,(x,y,z)) if pos2>(x,y,z)
                    expval += TNDecoder3D.identity_SU!(tn3d, (x,y,z), pos2, final_truncation_dim, svd_eps)
                    expval += TNDecoder3D.renormalize_tensor!(tn3d, (x,y,z))
                    expval += TNDecoder3D.renormalize_tensor!(tn3d, pos2)
                    expval += TNDecoder3D.renormalize_bond_matrix!(tn3d, (x,y,z), pos2)
                end end 
            end end end 
        end 
        tn3d_base = tn3d


        for i = 1:args["num_samples"]
            tn3d = deepcopy(tn3d_base)
            samples, logicals, _ = sampler.sample(1)
            detectors = convert(Vector{UInt8}, samples[1,:])
            logical_err = convert(UInt8, logicals[1])

            for i = 1:length(detectors)
                pos = pos_checks[i]
                if detectors[i] == 0
                    tn3d[pos] *= ITensor([1., 0.], open_index[pos])
                else
                    tn3d[pos] *= ITensor([0., 1.], open_index[pos])
                end
            end 
            for x=xmin:xmax for y=ymin:ymax for z=zmin:zmax
                if length(open_indices(tn3d, (x,y,z))) == 1
                    tn3d[(x,y,z)] *= ITensor([1.], open_indices(tn3d,(x,y,z))[1])
                end
            end end end

            for bd in args["bond_dimension"]
                res = contract_cubictn(deepcopy(tn3d), bd, args["mps_bond_dimension"]; svd_eps=svd_eps , split_dim=args["split_dimension"])
                dec_L = (res[1]>0.) ? 0 : 1

                if dec_L == logical_err
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
