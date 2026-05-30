# Standardized Beta DFA

Calculates Standardized Beta DFA

## Usage

``` r
sbdfa(data, dpo = 1, int = T, np = 91, overlap = T)
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

Scale s, Standardized Beta DFA estimates

## References

Barreto, I. D. C., Dore, L. H., Stosic, T. and Stosic, B. D. (2021).
Extending DFA-based multiple linear regression inference: application to
acoustic impedance models. *Physica A*, 582, 126259.

Kristoufek, L. (2015). Detrended fluctuation analysis as a regression
framework: estimating dependence at different scales. *Physical Review
E*, 91(2), 022802.

## See also

[`betadfa`](https://ikarobarreto.github.io/DFAmethods2/reference/betadfa.md),
[`vignette("DFAmethods2")`](https://ikarobarreto.github.io/DFAmethods2/articles/DFAmethods2.md)
