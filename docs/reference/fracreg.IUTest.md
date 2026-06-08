# Intersection-Union Test

Calculates Intersection-Union Test for Beta-DFA = 0 or Beta-DFA =
Beta-OLS

## Usage

``` r
fracreg.IUTest(data, B = 100, dpo = 1, int = TRUE, np = 91, overlap = TRUE)
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

A matrix with scale-wise Beta-DFA and critic region of Kristoufek Test
and Podobnik-Shen Test.

## References

Barreto, I. D. C., Dore, L. H., Stosic, T. and Stosic, B. D. (2021).
Extending DFA-based multiple linear regression inference: application to
acoustic impedance models. *Physica A*, 582, 126259.

## See also

[`fracreg.PStest`](https://ikarobarreto.github.io/DFATools/reference/fracreg.PStest.md),
[`fracreg.Ktest`](https://ikarobarreto.github.io/DFATools/reference/fracreg.Ktest.md),
[`vignette("DFATools")`](https://ikarobarreto.github.io/DFATools/articles/DFATools.md)

## Examples

``` r
set.seed(1)
d <- data.frame(y = cumsum(rnorm(250)), x = cumsum(rnorm(250)))
# \donttest{
fracreg.IUTest(d, B = 20, np = 15)
#> $PStest
#> # A tibble: 15 × 4
#>      bet1  slci1 suci1     s
#>     <dbl>  <dbl> <dbl> <int>
#>  1 -0.224 -0.259 0.215    10
#>  2 -0.244 -0.272 0.225    11
#>  3 -0.259 -0.282 0.232    12
#>  4 -0.271 -0.289 0.236    13
#>  5 -0.279 -0.294 0.239    14
#>  6 -0.285 -0.298 0.240    15
#>  7 -0.289 -0.303 0.240    16
#>  8 -0.292 -0.307 0.239    17
#>  9 -0.296 -0.317 0.234    19
#> 10 -0.300 -0.332 0.226    21
#> 11 -0.302 -0.362 0.241    23
#> 12 -0.302 -0.400 0.259    25
#> 13 -0.298 -0.457 0.296    28
#> 14 -0.292 -0.510 0.347    31
#> 15 -0.284 -0.557 0.393    34
#> 
#> $Ktest
#> # A tibble: 15 × 4
#>      bet1   klci1 kuci1     s
#>     <dbl>   <dbl> <dbl> <int>
#>  1 -0.224 -0.0252 0.439    10
#>  2 -0.244 -0.0310 0.452    11
#>  3 -0.259 -0.0332 0.462    12
#>  4 -0.271 -0.0319 0.470    13
#>  5 -0.279 -0.0280 0.481    14
#>  6 -0.285 -0.0251 0.494    15
#>  7 -0.289 -0.0276 0.506    16
#>  8 -0.292 -0.0309 0.517    17
#>  9 -0.296 -0.0405 0.539    19
#> 10 -0.300 -0.0685 0.559    21
#> 11 -0.302 -0.101  0.576    23
#> 12 -0.302 -0.131  0.591    25
#> 13 -0.298 -0.169  0.608    28
#> 14 -0.292 -0.201  0.620    31
#> 15 -0.284 -0.226  0.628    34
#> 
# }
```
