# sparse-ir-fortran - A Fortran interface for sparse-ir

This library provides a Fortran interface for `sparse-ir`.
At runtime, the installation of `sparse-ir` is not required.

Before you use it, please install `sparse-ir` with `xprec`. The python modules of `sparse-ir` will be used to generate datasets of IR-basis functions in advance. You can output the datasets in human-readable format or embed them into your Fortran source files.

## Creating a data file
### 1. Generate a data file:

```bash
> python3 dump.py 1e+4 1e-10 ir_nlambda4_ndigit10.dat
```

This generates the data file `ir_nlambda4_ndigit10.dat` containing sparse sampling points and transformation matrices for Λ=10<sup>4</sup> and ε<sup>-10</sup>  (ε is a cut-off value for singular values).

### 2. Build object files and link them to your program:

```bash
> gfortran  -c -o sparse_ir.o sparse_ir.f90
> gfortran  -c -o sparse_ir_io.o sparse_ir_io.f90
```

You do not need `sparse_ir_preset.f90` if your program reads a data file at runtime.

## Embedding data in a Fortran source file
Data can be embedded in a source file.
This allows you to avoid loading data files at runtime.

### 1. Generate a fortran source file:

The following command generates a source file containg data for a matrix of
Λ=10<sup>nlambda</sup> (nlambda = 1, 2, 3, 4) and ε=10<sup>-ndigit</sup> (ndigit = 10).

```bash
> python3 mk_preset.py --nlambda 1 2 3 4 --ndigit 10 > sparse_ir_preset.f90
```

### 2. Build object files and link them to your program:

The generated file `sparse_ir_preset.f90` is consisitent with Fortran95 and can be compiled as follows:

```bash
> gfortran  -c -o sparse_ir.o sparse_ir.f90
> gfortran  -c -o sparse_ir_preset.o sparse_ir_preset.f90
```

You do not need `sparse_ir_io.f90` if your program uses embedded data.

## How to use sparse-ir modules in your Fortran program
If your program reads data from the file, you should declare the use of the sparse-ir modules to initiate IR basis objects as follows:

```fortran
program main
  use sparse_ir
  use sparse_ir_io
```

Or, if you want to use the embed data of IR basis objects, you should do this instead:

```fortran
program main
  use sparse_ir
  use sparse_ir_preset
```

In either case, the derived type of "`IR`" should be declared.

```fortran
  implicit none
  type(IR) :: ir_obj ! You can declare with any name you wish.
```

Hereafter it is assumed that you will set the parameters as Λ=10<sup>4</sup>, β=10<sup>3</sup> and ε<sup>-10</sup>.
You can store the IR basis objects into the derived type of "`IR`" as follows:

Using `sparse_ir_io`:

```fortran
  double precision :: beta ! inverse temperature
  ...
  beta = 1.0d3
  open(99, file="ir_nlambda4_ndigit10.dat", status='old') ! Any unit number is OK.
  ir_obj = read_ir(99, beta)
```

Using `sparse_ir_preset`:

```fortran
  double precision :: beta ! inverse temperature
  ...
  beta = 1.0d3
  ir_obj = mk_ir_preset(4, 10, beta)
```

Here you are ready to use the IR basis objects and call the IR basis subroutines for a given value of `beta` in your program. If you want to use the IR basis sets with different values of `beta`, you can reset `beta` with using the subroutine `set_beta`:

```fortran
  beta = 2.0d3
  call set_beta(ir_obj, beta)
```

This subroutine updates the objects depending on `beta`.

## Available objects and subroutines

A detailed explanation of the IR basis objects and subroutines is shown [HERE](https://spm-lab.github.io/sparse-ir-fortran/).
