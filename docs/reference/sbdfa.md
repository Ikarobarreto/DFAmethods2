# Standardized Beta DFA

Calculates Standardized Beta DFA

## Usage

``` r
sbdfa(data, dpo = 1, int = TRUE, np = 91, overlap = TRUE)
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

[`betadfa`](https://ikarobarreto.github.io/DFATools/reference/betadfa.md),
[`vignette("DFATools")`](https://ikarobarreto.github.io/DFATools/articles/DFATools.md)

## Examples

``` r
set.seed(1)
d <- data.frame(y = cumsum(rnorm(300)), x1 = cumsum(rnorm(300)),
                x2 = cumsum(rnorm(300)))
sbdfa(d, np = 20)
#> # A tibble: 20 × 3
#>        s      x1      x2
#>    <int>   <dbl>   <dbl>
#>  1    10  0.0536 -0.0703
#>  2    12  0.0761 -0.0766
#>  3    13  0.0874 -0.0829
#>  4    14  0.0997 -0.0912
#>  5    15  0.114  -0.101 
#>  6    16  0.129  -0.112 
#>  7    17  0.147  -0.125 
#>  8    18  0.166  -0.139 
#>  9    19  0.186  -0.155 
#> 10    20  0.207  -0.171 
#> 11    21  0.228  -0.187 
#> 12    23  0.270  -0.220 
#> 13    25  0.311  -0.252 
#> 14    27  0.350  -0.283 
#> 15    29  0.386  -0.313 
#> 16    31  0.421  -0.341 
#> 17    34  0.467  -0.379 
#> 18    37  0.508  -0.412 
#> 19    40  0.543  -0.440 
#> 20    43 NA      NA     
```
