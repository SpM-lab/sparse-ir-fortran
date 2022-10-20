sparse-ir-fortran - A Fortran interface for sparse-ir
=====================================================

This library provides a Fortran interface for [sparse-ir].  

[sparse-ir]: https://github.com/SpM-lab/sparse-ir


Installation
------------
This package relies on the python package [sparse-ir] to generate the basis
functions:

    pip install sparse-ir[xprec]

**Note**, however, once the basis function are generated, this package is
standalone.  You can simply check out the repository and copy the relevant
source files into your code directory:

    git clone https://github.com/SpM-lab/sparse-ir-fortran


Documentation and tutorial
--------------------------
Check out our [comprehensive tutorial], where self-contained
notebooks for several many-body methods - GF(2), GW, Eliashberg equations,
Lichtenstein formula, FLEX, ... - are presented.  There is also a Fortran
sample code.

Refer to the [API documentation] for more details on how to work
with the Fortran library.

This library is currently somewhat restricted.  There is also a fully
fledged  [Python library] and a [Julia library] available for the IR basis
and sparse sampling.

[comprehensive tutorial]: https://spm-lab.github.io/sparse-ir-tutorial
[API documentation]: https://spm-lab.github.io/sparse-ir-fortran
[Python library]: https://github.com/SpM-lab/sparse-ir
[Julia library]: https://github.com/SpM-lab/SparseIR.jl


Creating a the basis functions
------------------------------
In order to work with sparse-ir-fortran, you first have to construct the
basis functions.  For this, you have two options: you can either create
data files, which you load at runtime, or you can embed the data into
a Fortran module itself.

### Creating data files

 1. Generate a data file:

        python dump.py 1e+4 1e-10 ir_nlambda4_ndigit10.dat

    This generates the data file `ir_nlambda4_ndigit10.dat` containing sparse
    sampling points and transformation matrices for Λ=10<sup>4</sup> and
    ε=10<sup>-10</sup>  (ε is a cut-off value used when creating IR-basis
    objects).

 2. Build object files and link them to your program:

        gfortran  -c -o sparse_ir.o sparse_ir.f90
        gfortran  -c -o sparse_ir_io.o sparse_ir_io.f90

    You do not need `sparse_ir_preset.f90` if your program reads a data file at
    runtime.

### Embedding data in aFortran source file

 1. Generate a fortran source file:  The following command generates a source
    file containg data for a matrix of Λ=10<sup>nlambda</sup>
    (nlambda = 1, 2, 3, 4) and ε=10<sup>-ndigit</sup> (ndigit = 10).

        python mk_preset.py --nlambda 1 2 3 4 --ndigit 10 > sparse_ir_preset.f90

 2. Build object files and link them to your program:  The generated file 
    `sparse_ir_preset.f90` is consistent with Fortran95 and can be compiled 
    as follows:

        gfortran  -c -o sparse_ir.o sparse_ir.f90
        gfortran  -c -o sparse_ir_preset.o sparse_ir_preset.f90

    You do not need `sparse_ir_io.f90` if your program uses embedded data.


How to use sparse-ir modules in your Fortran program
----------------------------------------------------
If your program reads data from the file, you should declare the use of the
sparse-ir modules to initiate IR basis objects as follows:

```fortran
program main
  use sparse_ir
  use sparse_ir_io
```

Or, if you want to use the embed data of IR basis objects, you should do this
instead:

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

Hereafter it is assumed that you will set the parameters as Λ=10<sup>4</sup>, 
β=10<sup>3</sup> and ε<sup>-10</sup>.  You can store the IR basis objects into
the derived type of "`IR`" as follows:

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

Here you are ready to use the IR basis objects and call the IR basis subroutines
for a given value of `beta` in your program. If you want to use the IR basis
sets with different values of `beta`, you can reset `beta` with using the
subroutine `set_beta`:

```fortran
  beta = 2.0d3
  call set_beta(ir_obj, beta)
```

This subroutine updates the objects depending on `beta`.

License and citation
--------------------
This software is released under the MIT License.  See LICENSE for details.

If you find the intermediate representation, sparse sampling, or this software
useful in your research, please consider citing the following papers:

 - Hiroshi Shinaoka et al., [Phys. Rev. B 96, 035147]  (2017)
 - Jia Li et al., [Phys. Rev. B 101, 035144] (2020)
 - Markus Wallerberger et al., [arXiv 2206.11762] (2022)

If you are discussing sparse sampling in your research specifically, please
also consider citing an independently discovered, closely related approach, the
MINIMAX isometry method (Merzuk Kaltak and Georg Kresse,
[Phys. Rev. B 101, 205145], 2020).

[Phys. Rev. B 96, 035147]: https://doi.org/10.1103/PhysRevB.96.035147
[Phys. Rev. B 101, 035144]: https://doi.org/10.1103/PhysRevB.101.035144
[arXiv 2206.11762]: https://doi.org/10.48550/arXiv.2206.11762
[Phys. Rev. B 101, 205145]: https://doi.org/10.1103/PhysRevB.101.205145
