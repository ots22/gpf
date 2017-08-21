# GPF

GPF is a small Fortran library for Gaussian process regression.  It currently
implements value predictions with dense Gaussian processes, and projected-process
approximate Gaussian processes (described by [Rasmussen and Williams, 2006, chapter 8](http://www.gaussianprocess.org/gpml/chapters/RW8.pdf)).

## Installation (Linux)

GPF has the following external dependencies

* LAPACK
* NLopt (http://ab-initio.mit.edu/wiki/index.php/NLopt)
* makedepf90

GPF is written in standard Fortran 2008.  A makefile is provided for gfortran on linux. 
It has been tested with gfortran 6.2, although earlier versions may also work.

Issuing
```sh
$ make all
```
will build the library (lib/libgpf.a) and tests.  The tests can be run with
```sh
$ make test
```

## Usage

### Examples

Some example programs using the library can be found in the directory `programs`.  See these 
directories for further details.
* gp_train: train (and optionally optimize) a sparse GP from data files
* gp_predict: read in a GP file and produce outputs for points read from standard input

### Linking with the library
Ensure that the module files `gp.mod`, `gp_dense.mod` and `gp_sparse.mod` are in the
include path of your compiler.  Link with the static library `libgpf.a` by providing 
the `-lgpf` flag (gfortran).

### Constructing a Gaussian process object from data
Include the module with either `use m_gp_dense` or `use m_gp_sparse`. These modules provide the types `gp_dense` (full Gaussian process) and `gp_sparse` (approximate Gaussian process using  A `gp` object can be constructed
as follows:
```f90
type(gp_dense) :: my_gp
type(cov_sqexp) :: cf
type(noise_value_only) :: nm
my_gp = DenseGP(nu=[1.d-9], theta=[1.4_dp], x=x(1:N,:), obs_type=obs_type(1:N), 
                t=t(1:N), CovFunction=cf, NoiseModel=nm)
```
where
...

### Reading a `gp` object from a file

```f90
DenseGP('in.gp')
```
or
```f90
SparseGP('in.gp')
```

### Writing a `gp` object to a file

```f90
my_gp%write_out('out.gp')
```

### Maximum liklihood optimization of the hyperparameters

```f90
use gp_optim
...
log_lik_optim(my_gp, lbounds, ubounds, optimize_max_iter, optimize_ftol)
```
#### Parameters:
* `my_gp`: an object of class BaseGP (either SparseGP or DenseGP)
* `lbounds` and `ubounds`: real arrays of length equal to the number of noise hyperparemeters
plus the number of covariance hyperparameters, representing lower and upper bounds of
the hyperparameters in the optimization.
* `optimize_max_iter`: maximum number of iterations before stopping.
* `optimize_ftol`: stop the optimization when an optimization step changes the log likelihood by less than this amount.

### Making predictions

A prediction from the process at a coordinate `x` can be obtained as
```f90
y = my_gp%predict(x, obs_type)
```
* `x`: a real array of the dimension of the process
* `obs_type`: an integer representing the type of the observation. A value of `0` gives a prediction of the value.
A value _i_ with _i > 1_, gives the predicted partial derivative of _y_ with respect to the _i_ th component of _x_.

### Covariance functions

### Noise models

### C bindings





