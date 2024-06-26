<p align="center">
    <img src="UVaTRMM_logo.png" alt="UvaTRMM logo">
</p>

# UVaTRMM: Distributed TRMM for commodity clusters

UVaTRMM is an open source C implementation for distributed TRMM (Triangular Matrix-Matrix Multiplication).
It leverages MPI to enable distributed memory executions.

Remark: By commodity clusters we mean small systems equipped with, at most, 40 nodes.

## Compiling UVaTRMM

```bash
$ cd UVaTRMM
$ make all
```

Dependencies: 

* MKL (or equivalent, for TRMM and GEMM kernels)
* MPI (remember to properly link the chosen library)

** Remember to create the directory "exec" for the binaries generated during compilation
 
## Executing UVaTRMM

```bash
$ mpirun -np <NumberOfMPIProcesses> exec/uvatrmm.x <MatrixDimension> <reg/balanced> <boxes/trapezoid>
```

## For comparison against ScaLAPACK

We provide ScaLAPACK GEMM and TRMM sample codes ready to be compiled and compared with our implementation.

Compile them with:
```bash
$ cd ScaLAPACK
$ make
```

## For comparison against SLATE

Please download the SLATE latest version from its repository (https://github.com/icl-utk-edu/slate)

# How to cite

If you write a scientific paper describing research that makes substantive use of UVaTRMM, we would appreciate that you cite the following paper:

* Performance improvement of the triangular matrix product in commodity clusters

```BibTeX
	@article{Carratala2023_UVaTRMM,
	author = "Santamaria-Valenzuela, Inmaculada and Carratal{\'a}-S{\'a}ez, Roc{\'i}o and Torres, Yuri and Llanos, Diego R. and Gonzalez-Escribano, Arturo",
	title="Performance improvement of the triangular matrix product in commodity clusters",
	journal="The Journal of Supercomputing",
	year="2024",
	issn="1573-0484",
	doi="10.1007/s11227-024-06097-7",
	}   
```
Paper link: https://rdcu.be/dIFn2
