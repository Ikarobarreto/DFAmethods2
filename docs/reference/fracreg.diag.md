# Scale-wise collinearity diagnostics for the fractal regression

Returns scale-dependent multicollinearity diagnostics for the DFA-based
multiple regression: the variance inflation factors (VIF), the condition
number (`kappa`) of the predictors' scale-wise correlation matrix, and
the scale-wise adjusted coefficient of determination. These reveal
multicollinearity that depends on the time scale.

## Usage

``` r
fracreg.diag(data, dpo = 1, int = TRUE, np = 91, overlap = TRUE, abs = FALSE)
```

## Arguments

- data:

  a matrix or data frame of time series; the first column is the
  response and the remaining columns are the predictors.

- dpo:

  detrending polynomial order.

- int:

  logical. If TRUE the integration process is applied.

- np:

  number of point scales.

- overlap:

  logical. If TRUE overlapping windows are applied.

- abs:

  logical. If TRUE the absolute detrended covariance is used.

## Value

A tibble with the scale `s`, one `VIF_*` column per predictor, the
condition number `kappa` and the adjusted `R2adj`, one row per scale.

## References

Tilfani, O., Kristoufek, L., Ferreira, P. and El Boukfaoui, M. Y.
(2022). Heterogeneity in economic relationships: scale dependence
through the multivariate fractal regression. *Physica A*, 588, 126530.

## See also

[`fracreg`](https://ikarobarreto.github.io/DFATools/reference/fracreg.md),
[`vignette("DFATools")`](https://ikarobarreto.github.io/DFATools/articles/DFATools.md)

## Examples

``` r
set.seed(1)
d <- data.frame(y = cumsum(rnorm(300)), x1 = cumsum(rnorm(300)),
                x2 = cumsum(rnorm(300)))
fracreg.diag(d, np = 20)
#> # A tibble: 20 × 5
#>        s VIF_x1 VIF_x2 kappa   R2adj
#>    <int>  <dbl>  <dbl> <dbl>   <dbl>
#>  1    10   1.02   1.02  1.36 -0.0669
#>  2    12   1.04   1.04  1.47 -0.0806
#>  3    13   1.05   1.05  1.53 -0.0874
#>  4    14   1.05   1.05  1.59 -0.0954
#>  5    15   1.06   1.06  1.64 -0.0981
#>  6    16   1.07   1.07  1.70 -0.109 
#>  7    17   1.08   1.08  1.75 -0.112 
#>  8    18   1.09   1.09  1.80 -0.115 
#>  9    19   1.10   1.10  1.85 -0.118 
#> 10    20   1.11   1.11  1.90 -0.108 
#> 11    21   1.11   1.11  1.94 -0.111 
#> 12    23   1.13   1.13  2.02 -0.103 
#> 13    25   1.14   1.14  2.08 -0.0940
#> 14    27   1.15   1.15  2.13 -0.0865
#> 15    29   1.16   1.16  2.17 -0.0824
#> 16    31   1.16   1.16  2.19 -0.0846
#> 17    34   1.16   1.16  2.19 -0.0779
#> 18    37   1.16   1.16  2.17 -0.0172
#> 19    40   1.15   1.15  2.13 -0.0264
#> 20    43   1.14   1.14  2.08 -0.0726
```
