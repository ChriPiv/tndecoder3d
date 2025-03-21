This repository contains source code used to generate the data for figure 1 of the paper _Tensor Network Decoding Beyond 2D_ by Christophe Piveteau, Christopher T. Chubb and Joseph M. Renes ([https://arxiv.org/abs/2310.10722](https://arxiv.org/abs/2310.10722)).

The code is encapsulated in a Julia package which can be loaded by including the file `tndecoder3d.jl`.
You will need a working Julia environment (version 1.8 was used in the paper) with the packages `ITensors` and `PyCall` installed.
If you're using a newer version of ITensors (>=0.7), you will also require the `ITensorMPS` package.
Following Python packages are used through PyCall and must be pre-installed:

* stim
* pymatching
* ldpc
* panqec

To see code examples demonstrating how to use the package (specifically for the point sector, loop sector and depolarizing noise on the unrotated 3D surface code as well as circuit-level noise on the rotated 3D surface code), see the scripts in the directory `scripts`.

