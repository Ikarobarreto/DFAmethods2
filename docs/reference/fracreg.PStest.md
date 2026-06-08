# Podobnik-Shen Test

Calculates Podobnik-Shen Test for Beta-DFA = 0

## Usage

``` r
fracreg.PStest(data, B = 100, dpo = 1, int = TRUE, np = 91, overlap = TRUE)
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

[`fracreg.Ktest`](https://ikarobarreto.github.io/DFATools/reference/fracreg.Ktest.md),
[`fracreg.IUTest`](https://ikarobarreto.github.io/DFATools/reference/fracreg.IUTest.md),
[`vignette("DFATools")`](https://ikarobarreto.github.io/DFATools/articles/DFATools.md)

## Examples

``` r
set.seed(1)
d <- data.frame(y = cumsum(rnorm(250)), x = cumsum(rnorm(250)))
# \donttest{
fracreg.PStest(d, B = 20, np = 15)
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
# }
```
