# rho Detrended Cross-Correlation Coefficient

Calculates rho DCCA

## Usage

``` r
rhodcca(data, dpo = 1, int = T, np = 91, overlap = T)
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

Scale s, Rho-DCCA

## References

Podobnik, B. and Stanley, H. E. (2008). Detrended cross-correlation
analysis: a new method for analyzing two nonstationary time series.
*Physical Review Letters*, 100(8), 084102.

Zebende, G. F. (2011). DCCA cross-correlation coefficient: quantifying
level of cross-correlation. *Physica A*, 390(4), 614-618.

## See also

[`vignette("DFAmethods2")`](https://ikarobarreto.github.io/DFAmethods2/articles/DFAmethods2.md)
