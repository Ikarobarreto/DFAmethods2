# rho Detrended Partial Cross-Correlation Coefficient

Calculates DPCCA

## Usage

``` r
rhodpcca(data, dpo = 1, int = TRUE, np = 91, overlap = TRUE)
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

[`rhodcca`](https://ikarobarreto.github.io/DFATools/reference/rhodcca.md),
[`vignette("DFATools")`](https://ikarobarreto.github.io/DFATools/articles/DFATools.md)

## Examples

``` r
set.seed(1)
d <- data.frame(x = cumsum(rnorm(300)), y = cumsum(rnorm(300)),
                z = cumsum(rnorm(300)))
rhodpcca(d, np = 20)
#> # A tibble: 20 × 4
#>        s DPCCA12 DPCCA13 DPCCA23
#>    <int>   <dbl>   <dbl>   <dbl>
#>  1    10  0.0531 -0.0696   0.154
#>  2    12  0.0748 -0.0753   0.195
#>  3    13  0.0856 -0.0813   0.215
#>  4    14  0.0973 -0.0891   0.233
#>  5    15  0.111  -0.0982   0.251
#>  6    16  0.125  -0.109    0.269
#>  7    17  0.142  -0.121    0.285
#>  8    18  0.160  -0.135    0.302
#>  9    19  0.178  -0.149    0.318
#> 10    20  0.198  -0.164    0.333
#> 11    21  0.217  -0.180    0.347
#> 12    23  0.256  -0.211    0.373
#> 13    25  0.294  -0.242    0.396
#> 14    27  0.330  -0.272    0.418
#> 15    29  0.365  -0.302    0.437
#> 16    31  0.398  -0.331    0.454
#> 17    34  0.443  -0.372    0.475
#> 18    37  0.484  -0.410    0.493
#> 19    40  0.522  -0.444    0.508
#> 20    43 NA      NA       NA    
```
