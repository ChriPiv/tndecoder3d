The scripts in this directory evaluate the performance of various decoders by running them on randomly generated error patterns:

* `sweep.jl`: 3D TN decoder on point or loop sector (specifiec by the `--sector` argument) of 3D surface code
* `sweep_matching.jl`: Matching decoder on point sector of 3D surface code
* `sweep_bposd.py`: BP+OSD on loop sector of 3D surface code
* `sweep_depol.jl`: 3D TN decoder on 3D surface code with depolarizing noise
* `sweep_depol_bposd.jl`: BP+OSD decoder on 3D surface code with depolarizing noise
* `sweep_circ.jl`: 3D TN decoder on rotated surface code with circuit-level noise
* `sweep_circ_matching.jl`: Matching decoder decoder on rotated surface code with circuit-level noise

The exact parameters accepted by the scripts vary slightly on the task performed.


To run the circuit-level TN decoder in `sweep_circ.jl`, you need to first generate the pre-compressed tensor network which is stored on the disk.
You can do this using the function in `precompress_circuit_tn.jl`.
Check out the two files to see where the pre-compressed network is stored, you will most likely have to change the file paths.



