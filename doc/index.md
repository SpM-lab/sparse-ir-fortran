# sparse-ir-fortran - A Fortran interface for sparse-ir

This library provides a Fortran interface for `sparse-ir`.
At runtime, the installation of `sparse-ir` is not required.

Before you use it, please install `sparse-ir` with `xprec`. The python modules of `sparse-ir` will be used to generate datasets of IR-basis functions in advance. You can output the datasets in human-readable format or embed them into your Fortran source files.

## Creating a data file
### 1. Generate a data file:

```bash
> python3 dump.py 1e+4 1e-10 ir_nlambda4_ndigit10.dat
```

This generates the data file `ir_nlambda4_ndigit10.dat` containing sparse sampling points and transformation matrices for $\Lambda=10^4$ and $\epsilon = 10^{-10}$ ($\epsilon$ is a cut-off value for singular values).

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
$\Lambda=10^{\mathrm{nlambda}}$ (nlambda = 1, 2, 3, 4) and $\epsilon=10^{-\mathrm{ndigit}}$ (ndigit = 10).

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

Hereafter it is assumed that you will set the parameters as $\Lambda = 10^4$, $\beta = 10^3$, and $\epsilon = 10^{-10}$.
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
## Available objects

### `DOUBLE PRECISION:: IR%beta`
It returns the input value $\beta$ of the functions `read_ir` or `mk_ir_preset`.

### `DOUBLE PRECISION:: IR%lambda`
It returns the input value $\Lambda$ of the functions `read_ir` or `mk_ir_preset`. $\Lambda$ determines the sparseness of the IR basis. In the "`sparse-ir-fortran`" interface, this value determines which dataset to extract among ones with different values of $\Lambda$.

### `DOUBLE PRECISION:: IR%wmax`
It returns the value of $\Lambda/\beta$, where $\beta$ is `IR%beta`.

### `DOUBLE PRECISION:: IR%eps`
It returns the value of $\epsilon=10^{-\mathrm{ndigit}}$, where $\mathrm{ndigit}$ is one of input variables of the functions `read_ir` or `mk_ir_preset`.
The cut-off value $\epsilon$ is used to determine how small a singular value should be considered meaningful when the kernels are SVDecomposed for generating the data of IR basis. In the "`sparse-ir-fortran`" interface, this value determines which dataset to extract among ones with different cut-offs.

### `DOUBLE PRECISION:: IR%s(IR%size)`
It returns the singular values of SVD for generating the data of IR basis.

### `INTEGER:: IR%size`
It returns the size of `IR%s`.

### `DOUBLE PRECISION:: IR%tau(IR%ntau)`
It returns the values of $\tau$ of sparse sampling points of imaginary time.

### `INTEGER:: IR%ntau`
It returns the size of `IR%tau`.

### `INTEGER:: IR%freq_f(IR%nfreq_f)`
It returns the odd integers corresponding to sampling Matsubara frequencies for fermionic functions.

### `INTEGER:: IR%nfreq_f`
It is the number of sampling Matsubara frequencies for fermionic functions.

### `INTEGER:: IR%freq_b(IR%nfreq_b)`
It returns the even integers corresponding to sampling Matsubara frequencies for bosonic functions.

### `INTEGER:: IR%nfreq_b`
It is the number of sampling Matsubara frequencies for bosonic functions.

### `TYPE(DecomposedMatrix):: IR%uhat_f`
It refers to the derived type of `DecomposedMatrix` which contains `IR%uhat_f%a`, `IR%uhat_f%inv_s`, `IR%uhat_f%ut`, and `IR%uhat_f%v`. When a derived type of `IR` is defined for a given `beta`, SVD of $\{\hat{U}_l(\mathrm{i}\nu_n)\}$ for a given `beta` is performed to define `IR%uhat_f%inv_s`, `IR%uhat_f%ut`, and `IR%uhat_f%v`, which are used in subroutines `fit_matsubara_f` and `evaluate_matsubara_f`. The basis functions on fermionic sampling Matsubara frequencies is SVDecomposed in advance  as follows:

$$
\begin{align*}
\hat{U}_l(\mathrm{i}\nu_n) &= A_{nl} \\
&=\sum_{r,r'}U_{nr} \Sigma_{rr'} (V^\mathrm{T})_{r'l} \\
(A &= U \Sigma V^\mathrm{T}),
\end{align*}
$$

with

$$
\Sigma_{rr'} = \sigma_r\delta_{r,r'}~  (\sigma_r \neq 0).
$$

If the following fitting problem is given for a certain fermionic function $G(\mathrm{i}\nu_n)$ defined on Matsubara frequencies,

$$
\sum_{l}\hat{U}_l(\mathrm{i}\nu_n)G_l\approx G(\mathrm{i}\nu_n),
$$

you can solve the problem using `fit_matsubara_f` as follows:

$$
G_l\approx \sum_{r, n}V_{lr}\Sigma^+_{rr}(U^\mathrm{T})_{rn}G(\mathrm{i}\nu_n).
$$

`IR%uhat_f%a`, `IR%uhat_f%ut`, and `IR%uhat_f%v` are 2-dimensional `COMPLEX(KIND(0D0))` arrays corresponding to $A$, $U^\mathrm{T}$, and $V$, respectively. `IR%uhat_f%inv_s` is the 1-dimensional `DOUBLE PRECISION` array storing the components of the diagonal matrix $\Sigma^+$, namely $\{1/\sigma_r\}$. `IR%uhat_f%ns` is the size of  `IR%uhat_f%inv_s`. `IR%uhat_f%m` and `IR%uhat_f%n` equal to `IR%nfreq_f` and `IR%size`, respectively.

### `TYPE(DecomposedMatrix):: IR%uhat_b`
It refers to the derived type of `DecomposedMatrix` which contains `IR%uhat_b%a`, `IR%uhat_b%inv_s`, `IR%uhat_b%ut`, and `IR%uhat_b%v`, which are used in subroutines `fit_matsubara_b` and `evaluate_matsubara_b`. `IR%u%a` corresponds to $\{\hat{U}_l(\mathrm{i}\nu_n)\}$ defined on bosonic Matsubara frequencies. This derived type is similar to `IR%uhat_f`, so please see the description of `IR%uhat_f` as well.

### `TYPE(DecomposedMatrix):: IR%u`
It refers to the derived type of `DecomposedMatrix` which contains `IR%u%a`, `IR%u%inv_s`, `IR%u%ut`, and `IR%u%v`, which are used in subroutines `fit_matsubara_b` and `evaluate_matsubara_b`. `IR%u%a` corresponds to $\{U_l(\tau_m)\}$ on sampling points of imaginary time. This derived type is similar to `IR%uhat_f`, so please see the description of `IR%uhat_f` as well.


## Subroutines

### `SUBROUTINE set_beta`
The subroutine (re)sets the value of `beta`.
`IR%u_data`, `IR%uhat_f_data`, and `IR%uhat_b_data` store the dimensionless forms of the IR-basis functions. This subroutine replaces the dimensionless variable `IR%x` with the `IR%beta` of them. This subroutine also do SVD of the basis functions.

The usage is

```fortran
call set_beta(obj, beta)
```

where `beta` is a `DOUBLE PRECISION` variable and `obj` is the derived type of "`IR`". In this subroutine, the objects in `obj` depending on `beta` are updated.


### `SUBROUTINE fit_matsubara_f`
The subroutine fits a set of expansion coefficients $G_l$ to a given fermionic function $G(\mathrm{i}\nu_n)$ on sampling Matsubara frequencies by using SVD.

$$
\begin{align*}
G_l = {\mathop{\rm argmin}\limits}_{G_l}\left|G(\mathrm{i}\nu_n) - \sum_{l}\hat{U}_l(\mathrm{i}\nu_n)G_l \right|^2
\end{align*}
$$

The usage is

```fortran
call fit_matsubara_f(obj, g_in, g_out)
```

where `g_in` and `g_out` correspond to $G(\mathrm{i}\nu_n)$ and $G_l$, respectively. `obj` is the derived type of "`IR`", and `g_in` and `g_out` are 2-dimensional `COMPLEX(KIND(0D0))` arrays. The inputs are `obj` and `g_in` and the output is `g_out`. Before calling this subroutine, you should reshape the array of $G(\mathrm{i}\nu_n)$ to a 2-dimensional array whose last axis corresponds to $l$ and allocate `g_out` with appropriate shape.
That is, `g_in` and `g_out` should be allocated so as to have shapes of `(**, obj%nfreq_f)` and `(**, obj%size)`, respectively.

### `SUBROUTINE fit_matsubara_b`
The subroutine fits a set of expansion coefficients $G_l$ to a given bosonic function $G(\mathrm{i}\nu_n)$ on sampling Matsubara frequencies by using SVD.

$$
\begin{align*}
G_l = {\mathop{\rm argmin}\limits}_{G_l}\left|G(\mathrm{i}\nu_n) - \sum_{l}\hat{U}_l(\mathrm{i}\nu_n)G_l \right|^2
\end{align*}
$$

The usage is

```fortran
call fit_matsubara_b(obj, g_in, g_out)
```

where `g_in` and `g_out` correspond to $G(\mathrm{i}\nu_n)$ and $G_l$, respectively. `obj` is the derived type of "`IR`", and `g_in` and `g_out` are 2-dimensional `COMPLEX(KIND(0D0))` arrays. The inputs are `obj` and `g_in` and the output is `g_out`. Before calling this subroutine, you should reshape the array of $G(\mathrm{i}\nu_n)$ to a 2-dimensional array whose last axis corresponds to $l$ and allocate `g_out` with appropriate shape.
That is, `g_in` and `g_out` should be allocated so as to have shapes of `(**, obj%nfreq_b)` and `(**, obj%size)`, respectively.

### `SUBROUTINE fit_tau`
The subroutine fits a set of expansion coefficients $G_l$ to a given imaginary-time function $G(\tau_m)$ on sampling points by using SVD.

$$
\begin{align*}
G_l = {\mathop{\rm argmin}\limits}_{G_l}\left|G(\tau_m) - \sum_{l}U_l(\tau_m)G_l \right|^2
\end{align*}
$$

The usage is

```fortran
call fit_tau(obj, g_in, g_out)
```

where `g_in` and `g_out` correspond to $G(\tau_m)$ and $G_l$, respectively. `obj` is the derived type of "`IR`", and `g_in` and `g_out` are 2-dimensional `COMPLEX(KIND(0D0))` arrays. The inputs are `obj` and `g_in` and the output is `g_out`. Before calling this subroutine, you should reshape the array of $G(\tau_m)$ to a 2-dimensional array whose last axis corresponds to $l$ and allocate `g_out` with appropriate shape.
That is, `g_in` and `g_out` should be allocated so as to have shapes of `(**, obj%ntau)` and `(**, obj%size)`, respectively.

### `SUBROUTINE evaluate_matsubara_f`
This subroutine reconstructs a fermionic function $G(\mathrm{i}\nu_n)$ on sampling Matsubara frequencies from a given set of expansion coefficients $G_l$ as follows:

$$
\begin{align*}
G(\mathrm{i}\nu_n) = \sum_{l}\hat{U}_l(\mathrm{i}\nu_n)G_l
\end{align*}
$$

The usage is

```fortran
call evaluate_matsubara_f(obj, g_in, g_out)
```

where `g_in` and `g_out` correspond to $G_l$ and $G(\mathrm{i}\nu_n)$, respectively. `obj` is the derived type of "`IR`", and `g_in` and `g_out` are 2-dimensional `COMPLEX(KIND(0D0))` arrays. The inputs are `obj` and `g_in` and the output is `g_out`. Before calling this subroutine, you should reshape the array of $G_l$ to a 2-dimensional array whose last axis corresponds to $\mathrm{i}\nu_n$ and allocate `g_out` with appropriate shape.
That is, `g_in` and `g_out` should be allocated so as to have shapes of `(**, obj%size)` and  `(**, obj%nfreq_f)`, respectively.

### `SUBROUTINE evaluate_matsubara_b`
This subroutine reconstructs a bosonic function $G(\mathrm{i}\nu_n)$ on sampling Matsubara frequencies from a given set of expansion coefficients $G_l$ as follows:

$$
\begin{align*}
G(\mathrm{i}\nu_n) = \sum_{l}\hat{U}_l(\mathrm{i}\nu_n)G_l
\end{align*}
$$

The usage is

```fortran
call evaluate_matsubara_b(obj, g_in, g_out)
```

where `g_in` and `g_out` correspond to $G_l$ and $G(\mathrm{i}\nu_n)$, respectively. `obj` is the derived type of "`IR`", and `g_in` and `g_out` are 2-dimensional `COMPLEX(KIND(0D0))` arrays. The inputs are `obj` and `g_in` and the output is `g_out`. Before calling this subroutine, you should reshape the array of $G_l$ to a 2-dimensional array whose last axis corresponds to $\mathrm{i}\nu_n$ and allocate `g_out` with appropriate shape.
That is, `g_in` and `g_out` should be allocated so as to have shapes of `(**, obj%size)` and  `(**, obj%nfreq_b)`, respectively.

### `SUBROUTINE evaluate_tau`
This subroutine reconstructs a function $G(\tau_m)$ on sampling Matsubara frequencies from a given set of expansion coefficients $G_l$ as follows:

$$
\begin{align*}
G(\tau_m) = \sum_{l}U_l(\tau_m)G_l
\end{align*}
$$

The usage is

```fortran
call evaluate_tau(obj, g_in, g_out)
```

where `g_in` and `g_out` correspond to $G_l$ and $G(\tau_m)$, respectively. `obj` is the derived type of "`IR`", and `g_in` and `g_out` are 2-dimensional `COMPLEX(KIND(0D0))` arrays. The inputs are `obj` and `g_in` and the output is `g_out`. Before calling this subroutine, you should reshape the array of $G_l$ to a 2-dimensional array whose last axis corresponds to $\mathrm{i}\nu_n$ and allocate `g_out` with appropriate shape.
That is, `g_in` and `g_out` should be allocated so as to have shapes of `(**, obj%size)` and  `(**, obj%ntau)`, respectively.
