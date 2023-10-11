include("../tndecoder3d.jl")
using .TNDecoder3D
using PyCall
using JLD2
@pyimport stim

function gen_circ_tn(d::Int, p::Float64, bd::Int)
    circuit = stim.Circuit.generated("surface_code:rotated_memory_x", distance=d, rounds=d, after_clifford_depolarization=p, before_measure_flip_probability=p, after_reset_flip_probability=p)
    model = circuit.detector_error_model(decompose_errors=false)
    dem = dem_from_stim(model)

    tn3d,expval,open_index = gen_tn_sc3d_circ(dem, bd; do_gauging=false, canonicalness_target=1e-4)

    save_object("/scratch/cpivetea/circtn"*string(d)*".jld2", (tn3d,expval,open_index))
    
end


gen_circ_tn(3, 0.01, 16)
