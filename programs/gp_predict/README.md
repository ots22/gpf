# gp_predict

Read a Gaussian process from a file `in.gp` in the current directory,
and compute predicted function values from coordinates read from
standard input.

#### Input (STDIN)

```
x_1 x_2 ... x_n obs_type
```

The _n_ floating point numbers `x_1` to `x_n` are an input coordinate,
where `n` must agree with the dimension of the Gaussian process read
from the file.

The integer `obs_type` is 0 to obtain a prediction of the underlying
function value _y_.  A value _i_ with 1 <= _i_ <= _n_ gives the predicted
partial derivative of _y_ with respect to `x_i`.

#### Output (STDOUT)

```
x_1 x_2 ... x_n obs_type y
```

`x_1` to `x_n` and obs_type are echoed from the input, and `y` is the
predicted function value or the partial derivative (depending on
obs_type as above)
