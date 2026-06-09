# Wild bootstrap inference for the fractal regression

Wild-bootstrap confidence intervals and significance for the
scale-dependent regression coefficients. It resamples the per-box
detrended-moment scores (the same scores used by
`fracreg(variance = "hc")`), so it is robust to heteroscedasticity and –
with dependent weights – to dependence between boxes, and it does not
recompute the DFA for each replicate.

## Usage

``` r
fracreg.WB(
  data,
  B = 999,
  weights = c("dependent", "rademacher", "mammen"),
  bandwidth = NULL,
  dpo = 1,
  int = TRUE,
  np = 91,
  min_boxes = 15,
  abs = FALSE
)
```

## Arguments

- data:

  a matrix or data frame; first column the response, the rest the
  predictors.

- B:

  number of bootstrap replicates.

- weights:

  multiplier-weight scheme: `"dependent"` (Shao 2010), `"rademacher"` or
  `"mammen"`.

- bandwidth:

  kernel bandwidth for `"dependent"` weights; defaults to \\T_s^{1/3}\\
  per scale.

- dpo:

  detrending polynomial order.

- int:

  logical. If TRUE the integration process is applied.

- np:

  number of point scales.

- min_boxes:

  minimum number of non-overlapping boxes \\T_s = \lfloor N/s\rfloor\\
  required for inference; at scales below this floor the bootstrap is
  skipped and the interval is returned as `NA` (default 15). Score-based
  inference treats the boxes as disjoint sampling units; `fracreg.WB()`
  therefore takes no `overlap` argument by design.

- abs:

  logical. If TRUE the absolute detrended covariance is used.

## Value

A tibble with the scale `s` and, per predictor, the estimate (`beta_*`),
the lower/upper interval bounds (`lower_*`, `upper_*`) and the p-value
(`p_*`), one row per scale.

## Details

For each scale \\s\\ with \\T_s\\ boxes, scores \\r_v\\ and inverse
\\F\_{XX}(s)^{-1}\\, each replicate is \\\hat\beta^\*\_b(s) =
\hat\beta(s) + F\_{XX}(s)^{-1}\\(1/T_s)\sum_v r_v w\_{v,b}\\. The
interval is the 2.5% / 97.5% quantiles of \\\hat\beta^\*\_b(s)\\; the
two-sided p-value tests \\H_0\\: \beta_j(s) = 0\\. The `"dependent"`
weights (Shao 2010) keep coverage when the boxes are dependent (e.g.
strong long memory), where the analytic t and the i.i.d. weights can
fail.

## References

Shao, X. (2010). The dependent wild bootstrap. *Journal of the American
Statistical Association*, 105(489), 218-235.

Mammen, E. (1993). Bootstrap and wild bootstrap for high dimensional
linear models. *The Annals of Statistics*, 21(1), 255-285.

Tilfani, O., Kristoufek, L., Ferreira, P. and El Boukfaoui, M. Y.
(2022). Heterogeneity in economic relationships: scale dependence
through the multivariate fractal regression. *Physica A*, 588, 126530.

## See also

[`fracreg`](https://ikarobarreto.github.io/DFATools/reference/fracreg.md),
[`fracreg.IUTest`](https://ikarobarreto.github.io/DFATools/reference/fracreg.IUTest.md),
[`vignette("DFATools")`](https://ikarobarreto.github.io/DFATools/articles/DFATools.md)

## Examples

``` r
set.seed(1)
d <- data.frame(y = cumsum(rnorm(300)), x1 = cumsum(rnorm(300)),
                x2 = cumsum(rnorm(300)))
# \donttest{
fracreg.WB(d, B = 199, np = 15)
#> Warning: fracreg.WB(): N = 300 < 500; DFA-based inference may be unreliable (Likens et al., 2019).
#> Warning: fracreg.WB(): scales {21, 23, 26, 29, 32, ...} have fewer than min_boxes = 15 non-overlapping boxes; their standard errors and intervals were set to NA. Reduce the maximum scale or increase N.
#> # A tibble: 15 × 9
#>        s beta_x1 lower_x1 upper_x1    p_x1 beta_x2 lower_x2  upper_x2    p_x2
#>    <int>   <dbl>    <dbl>    <dbl>   <dbl>   <dbl>    <dbl>     <dbl>   <dbl>
#>  1    10  0.111   -0.130     0.344  0.362  -0.0824   -0.272  0.101     0.412 
#>  2    11  0.0775  -0.0613    0.236  0.271  -0.431    -0.644 -0.209     0     
#>  3    12  0.178   -0.0670    0.468  0.201  -0.247    -0.456 -0.0171    0.0302
#>  4    13  0.0562  -0.349     0.501  0.714  -0.0205   -0.261  0.287     0.945 
#>  5    14  0.285   -0.0259    0.542  0.0603 -0.205    -0.600  0.111     0.241 
#>  6    15  0.0189  -0.366     0.363  0.894  -0.240    -0.605  0.126     0.332 
#>  7    17  0.116   -0.201     0.515  0.503  -0.0982   -0.422  0.205     0.663 
#>  8    19  0.375    0.0779    0.632  0.0201 -0.174    -0.391 -0.000576  0.0503
#>  9    21  0.266   NA        NA     NA      -0.0509   NA     NA        NA     
#> 10    23  0.380   NA        NA     NA      -0.0292   NA     NA        NA     
#> 11    26  0.617   NA        NA     NA      -0.342    NA     NA        NA     
#> 12    29  0.419   NA        NA     NA      -0.0880   NA     NA        NA     
#> 13    32  0.580   NA        NA     NA      -0.287    NA     NA        NA     
#> 14    36  0.455   NA        NA     NA      -0.597    NA     NA        NA     
#> 15    40  0.986   NA        NA     NA      -0.763    NA     NA        NA     
# }
```
