# Detrended Fluctuation Analysis

Calculates DFA

## Usage

``` r
dfa(data, dpo = 1, int = TRUE, np = 91, overlap = TRUE)
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

[`vignette("DFATools")`](https://ikarobarreto.github.io/DFATools/articles/DFATools.md)

## Examples

``` r
set.seed(1)
x <- cumsum(rnorm(300))
dfa(x, np = 20)
#> # A tibble: 20 × 2
#>        s      F
#>    <int>  <dbl>
#>  1    10   1.94
#>  2    12   3.29
#>  3    13   4.17
#>  4    14   5.19
#>  5    15   6.37
#>  6    16   7.74
#>  7    17   9.30
#>  8    18  11.1 
#>  9    19  13.1 
#> 10    20  15.3 
#> 11    21  17.8 
#> 12    23  23.7 
#> 13    25  30.8 
#> 14    27  39.3 
#> 15    29  49.2 
#> 16    31  60.8 
#> 17    34  81.5 
#> 18    37 106.  
#> 19    40 136.  
#> 20    43 170.  
```
