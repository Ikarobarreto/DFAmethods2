# Kristoufek Test

Calculates Kristoufek Test for Beta-DFA = Beta-OLS

## Usage

``` r
fracreg.Ktest(data, B = 100, dpo = 1, int = TRUE, np = 91, overlap = TRUE)
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

[`fracreg.PStest`](https://ikarobarreto.github.io/DFATools/reference/fracreg.PStest.md),
[`fracreg.IUTest`](https://ikarobarreto.github.io/DFATools/reference/fracreg.IUTest.md),
[`vignette("DFATools")`](https://ikarobarreto.github.io/DFATools/articles/DFATools.md)

## Examples

``` r
set.seed(1)
d <- data.frame(y = cumsum(rnorm(250)), x = cumsum(rnorm(250)))
# \donttest{
fracreg.Ktest(d, B = 20, np = 15)
#> # A tibble: 15 × 4
#>      bet1   klci1 kuci1     s
#>     <dbl>   <dbl> <dbl> <int>
#>  1 -0.224 -0.0536 0.571    10
#>  2 -0.244 -0.0570 0.582    11
#>  3 -0.259 -0.0568 0.588    12
#>  4 -0.271 -0.0626 0.590    13
#>  5 -0.279 -0.0672 0.590    14
#>  6 -0.285 -0.0715 0.587    15
#>  7 -0.289 -0.0760 0.582    16
#>  8 -0.292 -0.0807 0.576    17
#>  9 -0.296 -0.108  0.566    19
#> 10 -0.300 -0.147  0.560    21
#> 11 -0.302 -0.186  0.563    23
#> 12 -0.302 -0.226  0.577    25
#> 13 -0.298 -0.282  0.600    28
#> 14 -0.292 -0.332  0.625    31
#> 15 -0.284 -0.375  0.660    34
# }
```
