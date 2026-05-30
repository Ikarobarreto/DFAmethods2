# rho Detrended Multiple-Correlation Coefficient

Calculates DMC²

## Usage

``` r
dmc2(data, dpo = 1, int = T, np = 91, overlap = T)
```

## Arguments

- data:

  is a matrix of time series

- dpo:

  detrending polynomial order

- int:

  logical. if TRUE integration process will be applied.

- np:

  number of point scales.

- overlap:

  logical. if TRUE overlapping windows will be applied.

## Value

Scale s, Multiple Detrended Correlation

## References

Zebende, G. F. and da Silva Filho, A. M. (2018). Detrended multiple
cross-correlation coefficient. *Physica A*, 510, 91-97.

Wang, F., Xu, J. and Fan, Q. (2021). Statistical test for detrended
multiple cross-correlation coefficient. *Communications in Nonlinear
Science and Numerical Simulation*, 99, 105781.

## See also

[`vignette("DFAmethods2")`](https://ikarobarreto.github.io/DFAmethods2/articles/DFAmethods2.md)
