# sparse-ir-fortran - A Fortran interface for sparse-ir

This library provides a Fortran interface for `sparse-ir`.
At runtime, the installation of `sparse-ir` is not required.

First, please install `sparse-ir` with `xprec`.

## Creating a data file
1. Generate a data file:

```bash
> python3 dump.py 1e+4 1e-10 ir_nlambda4_ndigit10.dat
```

This generate the data file `ir_nlambda4_ndigit10.dat` containing sparse sampling points and transformation matrices for $\Lambda=10^4$ and $\epsilon = 10^{-10}$ ($\epsilon$ is a cut-off value for singular values).

2. Build object files and link them to your program:

```bash
> gfortran  -c -o sparse_ir.o sparse_ir.f90
> gfortran  -c -o sparse_ir_io.o sparse_ir_io.f90
```

You do not need `sparse_ir_preset.f90` if your program reads a data file at runtime.

## Embedding data in a Fortran source file
Data can be embedded in a source file.
This allows you to avoid loading data files at runtime.

The following command generates a source file containg data for a matrix of
$\Lambda=10^{-\mathrm{nlambda}}$ ($\Lambda=1,2,3,4$) and $\epsilon=10^{-\mathrm{ndigit}}$ (ndigit = 10).

```bash
> python3 mk_preset.py --nlambda 1 2 3 4 --ndigit 10 > sparse_ir_preset.f90
```

The generated file is consisitent with Fortran95 and can be compiled as follows:

```bash
> gfortran  -c -o sparse_ir.o sparse_ir.f90
> gfortran  -c -o sparse_ir_preset.o sparse_ir_preset.f90
```

## Using  from your Fortran program
Under construction
