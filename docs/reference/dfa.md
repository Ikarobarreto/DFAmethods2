# Detrended Fluctuation Analysis

Calculates DFA

## Usage

``` r
dfa(data, dpo = 1, int = T, np = 91, overlap = T)
```

## Arguments

- data:

  is a vector of time series

- dpo:

  detrending polynomial order

- int:

  logical. if TRUE integration process will be applied.

- np:

  number of point scales.

- overlap:

  logical. if TRUE overlapping windows will be applied.

## Value

Scale s, Detrended Fluctuation Function F

## References

Peng, C.-K. et al. (1994). Mosaic organization of DNA nucleotides.
*Physical Review E*, 49(2), 1685-1689.

## See also

[`vignette("DFAmethods2")`](https://ikarobarreto.github.io/DFAmethods2/articles/DFAmethods2.md)
