# Multiple Fractal Regression

Calculates the scale-dependent (DFA-based) multiple linear regression:
the coefficients, their variance and confidence intervals at each scale,
together with scale-wise collinearity diagnostics.

## Usage

``` r
fracreg(
  data,
  dpo = 1,
  int = TRUE,
  np = 91,
  overlap = FALSE,
  variance = c("inv_corrected", "inv", "marginal", "hc", "none"),
  H_eps = NULL,
  min_boxes = 15,
  abs = FALSE
)
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

  logical. If TRUE overlapping windows are used. Defaults to FALSE:
  score-based inference treats the boxes as sampling units and requires
  them disjoint. With `overlap = TRUE` only point estimates are returned
  (`variance = "none"`).

- variance:

  coefficient-variance estimator: `"inv_corrected"` (default) the
  memory-corrected inverse form
  \\F^2\_\varepsilon(s)\[F\_{XX}(s)^{-1}\]\_{jj}/(T_s-k)\cdot(2\widehat
  H+1)^2/21\\ (Barreto et al. 2026, Eq. M4.2'); `"inv"` the uncorrected
  inverse (Tilfani et al. 2022); `"marginal"` the legacy
  \\1/F^2\_{X_j}(s)\\ form (Shen 2015), which under-covers under
  collinearity; `"hc"` the heteroscedasticity-consistent sandwich; or
  `"none"` for point estimates only.

- H_eps:

  optional numeric. A pre-computed DFA exponent of the regression error
  to use in the memory correction; if `NULL` (default) it is estimated
  from the OLS residual.

- min_boxes:

  minimum number of non-overlapping boxes \\T_s = \lfloor N/s\rfloor\\
  required for inference at a scale (default 15). Scales with \\T_s \<
  \\ `min_boxes` return `NA` standard errors, confidence limits and
  p-values, with a single warning listing them. A warning is also issued
  when \\N \< 500\\ (Likens et al. 2019).

- abs:

  logical. If TRUE the absolute detrended covariance is used in the
  cross-correlation step (more robust to outliers).

## Value

A list with, among others: scale `s`, detrended fluctuation `F`, `DCCA`,
`DPCCA`, beta estimates `BDFA`, standardized betas `BSDFA`, residual
variance `UDFA`, coefficient variance `VDFA`, multiple correlation
`DMC2`, `R2DFA`, confidence limits `UCIB`/`LCIB`, `p.value`, critical
value `TC`, the scale-wise diagnostics `VIF`, condition number `kappa`
and adjusted `R2adj`, the per-series DFA exponent `alpha`, the residual
DFA exponent `H_resid` and memory factor `kappa_factor` used in the
correction, the chosen `variance_method`, the effective degrees of
freedom `df_eff` and the box counts `Ts`.

## Details

The variance of the scale-dependent coefficients follows Tilfani et al.
(2022, Eq. 25), \\\mathrm{var}(\hat\beta_j(s)) = F^2\_\varepsilon(s)
\[F\_{XX}(s)^{-1}\]\_{jj}\\, i.e. it uses the full inverse of the
detrended covariance matrix of the predictors, consistent with how the
coefficients themselves are estimated. Under collinearity this differs
from the legacy "marginal" form \\F^2\_\varepsilon(s) / F^2\_{X_j}(s)\\
(which under-covers); the two coincide for orthogonal predictors. The
default `variance = "inv_corrected"` multiplies the inverse form by the
memory-correction factor \\(2\widehat H + 1)^2 / 21\\ (Barreto et al.
2026, Eq. M4.2'), with \\\widehat H\\ the DFA exponent of the OLS
residual, restoring nominal coverage under polynomial detrending;
`variance = "inv"` omits the factor and `variance = "marginal"`
reproduces the legacy form.

The variance is normalised by the residual degrees of freedom \\T_s -
k\\, where \\T_s = \lfloor N/s \rfloor\\ is the number of
non-overlapping boxes at scale \\s\\ and \\k\\ the number of predictors;
the same \\T_s - k\\ is the degrees of freedom of the `t` quantile.
\\T_s\\ counts disjoint boxes regardless of `overlap`. Scales with \\T_s
\le k\\ return `NA` limits with a warning. The analytic interval can
under-cover under strong long-range dependence (estimated DFA exponent
above \\3/4\\, the Hermite-Rosenblatt threshold); a warning is then
issued and
[`fracreg.WB`](https://ikarobarreto.github.io/DFATools/reference/fracreg.WB.md)
(dependent wild bootstrap) should be preferred.

## References

Barreto, I. D. C., Dore, L. H., Stosic, T. and Stosic, B. D. (2021).
Extending DFA-based multiple linear regression inference: application to
acoustic impedance models. *Physica A*, 582, 126259.

Tilfani, O., Kristoufek, L., Ferreira, P. and El Boukfaoui, M. Y.
(2022). Heterogeneity in economic relationships: scale dependence
through the multivariate fractal regression. *Physica A*, 588, 126530.

Shen, C. (2015). A new detrended semipartial cross-correlation analysis.
*Physics Letters A*, 379(44), 2962-2969.

Kristoufek, L. (2015). Detrended fluctuation analysis as a regression
framework: estimating dependence at different scales. *Physical Review
E*, 91(2), 022802.

## See also

[`fracreg.diag`](https://ikarobarreto.github.io/DFATools/reference/fracreg.diag.md),
[`betadfa`](https://ikarobarreto.github.io/DFATools/reference/betadfa.md),
[`effsizeDFA`](https://ikarobarreto.github.io/DFATools/reference/effsizeDFA.md),
[`vignette("DFATools")`](https://ikarobarreto.github.io/DFATools/articles/DFATools.md)

## Examples

``` r
set.seed(1)
x1 <- rnorm(400); x2 <- rnorm(400)        # stationary predictors
d <- data.frame(y = 0.7 * x1 - 0.5 * x2 + rnorm(400), x1 = x1, x2 = x2)
fit <- fracreg(d, dpo = 1, int = TRUE, np = 20, overlap = FALSE)
#> Warning: fracreg(): N = 400 < 500; DFA-based inference may be unreliable (Likens et al., 2019).
#> Warning: fracreg(): scales {29, 32, 35, 39, 43, ...} have fewer than min_boxes = 15 non-overlapping boxes (T_s = floor(N/s)); their standard errors and intervals were set to NA. Reduce the maximum scale or increase N.
round(fit$BDFA[, 1, 10], 2)               # coefficients at the 10th scale
#> [1]  0.66 -0.34
```
