# rho Detrended Cross-Correlation Coefficient

Calculates rho DCCA

## Usage

``` r
rhodcca(data, dpo = 1, int = TRUE, np = 91, overlap = TRUE)
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

[`vignette("DFATools")`](https://ikarobarreto.github.io/DFATools/articles/DFATools.md)

## Examples

``` r
set.seed(1)
d <- data.frame(x = cumsum(rnorm(300)), y = cumsum(rnorm(300)))
rhodcca(d, np = 20)
#> # A tibble: 20 × 2
#>        s  DCCA12
#>    <int>   <dbl>
#>  1    10  0.0430
#>  2    12  0.0614
#>  3    13  0.0700
#>  4    14  0.0790
#>  5    15  0.0891
#>  6    16  0.100 
#>  7    17  0.113 
#>  8    18  0.126 
#>  9    19  0.140 
#> 10    20  0.154 
#> 11    21  0.168 
#> 12    23  0.196 
#> 13    25  0.222 
#> 14    27  0.247 
#> 15    29  0.271 
#> 16    31  0.294 
#> 17    34  0.326 
#> 18    37  0.356 
#> 19    40  0.384 
#> 20    43 NA     
```
