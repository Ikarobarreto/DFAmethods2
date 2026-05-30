# Beta DFA

Calculates Beta DFA

## Usage

``` r
betadfa(data, dpo = 1, int = T, np = 91, overlap = T)
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

Scale s, Beta DFA estimates

## References

Kristoufek, L. (2015). Detrended fluctuation analysis as a regression
framework: estimating dependence at different scales. *Physical Review
E*, 91(2), 022802.

## See also

[`sbdfa`](https://ikarobarreto.github.io/DFAmethods2/reference/sbdfa.md),
[`fracreg`](https://ikarobarreto.github.io/DFAmethods2/reference/fracreg.md),
[`vignette("DFAmethods2")`](https://ikarobarreto.github.io/DFAmethods2/articles/DFAmethods2.md)
