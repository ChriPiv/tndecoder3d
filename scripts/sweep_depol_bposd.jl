using MKL
using ArgParse
using JSON
using PyCall
include("../tndecoder3d.jl")
using .TNDecoder3D
@pyimport ldpc
@pyimport panqec.codes as pqc

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

data = Dict{Tuple{Int,Int}, Tuple{Int,Int}}() # (distance,err_rate,bond_dim) -> (num_succ,num_fail)
for d in args["distance"] for p_int in args["error_rate"]
    data[d,p_int] = (0,0)
end end


function build_dem(sc3d::PyObject, p::Float64)
    Hx = convert(Matrix{Int8}, sc3d.Hx.todense())
    Hz = convert(Matrix{Int8}, sc3d.Hz.todense())
    n_stab_x = size(Hx,1)
    n_stab_z = size(Hz,1)
    H = zeros(Int8, n_stab_x+n_stab_z, 3*sc3d.n)

    # depol noise [1-p, p/3, p/3, p/3] corresponds to independent X,Y,Z noise
    # with strength q:
    q = 0.5 - sqrt(9. - 12*p) / 6.
    p_vec = q * ones(Float64, 3*sc3d.n)

    for i = 0:sc3d.n-1
        H[:,3*i+1] = [zeros(Int8, n_stab_x) ; Hz[:,i+1]]        # X error
        H[:,3*i+2] = [Hx[:,i+1] ; Hz[:,i+1]]                    # Y error
        H[:,3*i+3] = [Hx[:,i+1] ; zeros(Int8, n_stab_z)]        # Z error
    end

    return H, p_vec
end

for d in args["distance"]
    sc3d = pqc.surface_3d.Planar3DCode(d,d,d)
    for p_int in args["error_rate"]
        p = p_int / 10000.
        
        Hx = convert(Matrix{Int8}, sc3d.Hx.todense())
        Hz = convert(Matrix{Int8}, sc3d.Hz.todense())
        Lx = convert(Vector{Int8}, sc3d.logicals_x[1,:])[1:sc3d.n]
        Lz = convert(Vector{Int8}, sc3d.logicals_z[1,:])[sc3d.n+1:end]

        if d == 3
            order = 19
        else
            order = 60
        end
        decX = ldpc.bposd_decoder(
            Hx,
            channel_probs=ones(sc3d.n)*2*p/3,
            max_iter=sc3d.n,
            bp_method="ms",
            ms_scaling_factor=0, #min sum scaling factor. If set to zero the variable scaling factor method is used
            osd_method="osd_cs", #the OSD method. Choose from:  1) "osd_e", "osd_cs", "osd0"
            osd_order=order #the osd search depth
        )
        decZ = ldpc.bposd_decoder(
            Hz,
            channel_probs=ones(sc3d.n)*2*p/3,
            max_iter=sc3d.n,
            bp_method="ms",
            ms_scaling_factor=0, #min sum scaling factor. If set to zero the variable scaling factor method is used
            osd_method="osd_cs", #the OSD method. Choose from:  1) "osd_e", "osd_cs", "osd0"
            osd_order=order #the osd search depth
        )
        
        for i = 1:args["num_samples"]
            # to sample depol noise, sample a random number in [0,1]
            # range (0,p/3) : X error
            # range (p/3,2*p/3) : Y error
            # range (2p/3,p) : Z error
            # range (p,1) : no error
            err = rand(sc3d.n)
            err_x = err .< (p * 2/3)
            err_z = (err .> (p * 1/3)) .& (err .< p)
            syndrome_x = (Hx * err_z) .% 2
            syndrome_z = (Hz * err_x) .% 2
            err_Lx = sum(err_z .& Lx) % 2 
            err_Lz = sum(err_x .& Lz) % 2

            #error_estimate_x = decZ.decode(syndrome_z)
            #new_probs = zeros(sc3d.n)
            #for i = 1:sc3d.n
            #    new_probs[i] = (error_estimate_x[i]==1) ? 0.5 : (p/3) / (1. - 2*p/3) 
            #end
            #decX.update_channel_probs(new_probs)
            #error_estimate_z = decX.decode(syndrome_x)

            error_estimate_z = decX.decode(syndrome_x)
            new_probs = zeros(sc3d.n)
            for i = 1:sc3d.n
                new_probs[i] = (error_estimate_z[i]==1) ? 0.5 : (p/3) / (1. - 2*p/3) 
            end
            decZ.update_channel_probs(new_probs)
            error_estimate_x = decZ.decode(syndrome_z)

            dec_LX = sum( error_estimate_z .& Lx) % 2
            dec_LZ = sum( error_estimate_x .& Lz) % 2

            if dec_LX == err_Lx && dec_LZ == err_Lz
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
