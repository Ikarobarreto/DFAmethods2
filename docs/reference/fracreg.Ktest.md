# Kristoufek Test

Calculates Kristoufek Test for Beta-DFA = Beta-OLS

## Usage

``` r
fracreg.Ktest(data, B = 100, dpo = 1, int = T, np = 91, overlap = T)
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

A matrix with scale-wise Beta-DFA and critic region of Kristoufek Test.

## References

Kristoufek, L. (2015). Detrended fluctuation analysis as a regression
framework: estimating dependence at different scales. *Physical Review
E*, 91(2), 022802.

## See also

[`fracreg.PStest`](https://ikarobarreto.github.io/DFAmethods2/reference/fracreg.PStest.md),
[`fracreg.IUTest`](https://ikarobarreto.github.io/DFAmethods2/reference/fracreg.IUTest.md),
[`vignette("DFAmethods2")`](https://ikarobarreto.github.io/DFAmethods2/articles/DFAmethods2.md)
