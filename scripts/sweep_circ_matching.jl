using JLD2
using MKL
using ArgParse
using JSON
include("../tndecoder3d.jl")
using .TNDecoder3D
using PyCall
using ITensors
@pyimport stim
@pyimport pymatching 

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
"--output_write_interval"
    arg_type = Int
    help = "after how many samples should the output file be updated"
    default = 20
"--num_samples", "-n"
    arg_type = Int
    help = "number of Monte Carlo samples"
    default = 1_000_000
end
args = parse_args(s)


data = Dict{Tuple{Int,Int}, Tuple{Int,Int}}() # (distance,err_rate) -> (num_succ,num_fail)
for d in args["distance"] for p_int in args["error_rate"]
data[d,p_int] = (0,0)
end end 

for d in args["distance"]
    for p_int in args["error_rate"]
        p = p_int / 10000.
        circuit = stim.Circuit.generated("surface_code:rotated_memory_x", distance=d, rounds=d, after_clifford_depolarization=p, before_measure_flip_probability=p, after_reset_flip_probability=p)
        true_model = circuit.detector_error_model(decompose_errors=false)
        approximate_model = circuit.detector_error_model(decompose_errors=true)
        sampler = true_model.compile_sampler()
        matching = pymatching.Matching.from_detector_error_model(approximate_model)

        for i = 1:args["num_samples"]
            samples, logicals, _ = sampler.sample(1)
            detectors = convert(Vector{Int8}, samples[1,:])
            logical_err = convert(Int8, logicals[1])

            predicted_observables = matching.decode(detectors)

            if logical_err == predicted_observables[1]
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
