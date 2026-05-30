# rho Detrended Partial Cross-Correlation Coefficient

Calculates DPCCA

## Usage

``` r
rhodpcca(data, dpo = 1, int = T, np = 91, overlap = T)
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

Scale s, Rho DPCCA

## References

Yuan, N. et al. (2015). Detrended partial-cross-correlation analysis: a
new method for analyzing correlations in complex system. *Scientific
Reports*, 5, 8143.

## See also

[`rhodcca`](https://ikarobarreto.github.io/DFAmethods2/reference/rhodcca.md),
[`vignette("DFAmethods2")`](https://ikarobarreto.github.io/DFAmethods2/articles/DFAmethods2.md)
