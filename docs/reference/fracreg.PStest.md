# Podobnik-Shen Test

Calculates Podobnik-Shen Test for Beta-DFA = 0

## Usage

``` r
fracreg.PStest(data, B = 100, dpo = 1, int = T, np = 91, overlap = T)
```

## Arguments

- data:

  is a matrix of time series

- B:

  number of surrogate series

- dpo:

  detrending polynomial order

- int:

  logical. if TRUE integration process will be applied.

- np:

  number of point scales.

- overlap:

  logical. if TRUE overlapping windows will be applied.

## Value

A matrix with scale-wise Beta-DFA and critic region of Podobnik-Shen
Test.

## References

Podobnik, B., Jiang, Z.-Q., Zhou, W.-X. and Stanley, H. E. (2011).
Statistical tests for power-law cross-correlated processes. *Physical
Review E*, 84(6), 066118.

Shen, C. (2015). A new detrended semipartial cross-correlation analysis.
*Physics Letters A*, 379(44), 2962-2969.

## See also

[`fracreg.Ktest`](https://ikarobarreto.github.io/DFAmethods2/reference/fracreg.Ktest.md),
[`fracreg.IUTest`](https://ikarobarreto.github.io/DFAmethods2/reference/fracreg.IUTest.md),
[`vignette("DFAmethods2")`](https://ikarobarreto.github.io/DFAmethods2/articles/DFAmethods2.md)
