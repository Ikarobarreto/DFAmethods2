# Multiple Fractal Regression

Calculates the scale-dependent (DFA-based) multiple linear regression:
the coefficients, their variance and confidence intervals at each scale,
together with scale-wise collinearity diagnostics.

## Usage

``` r
fracreg(
  data,
  dpo,
  int,
  np = 91,
  overlap = TRUE,
  vcov = c("inverse", "marginal", "HC"),
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

  logical. If TRUE overlapping windows are applied.

- vcov:

  coefficient-variance estimator: `"inverse"` (default; full
  \\F\_{XX}(s)^{-1}\\, Tilfani et al. 2022), `"marginal"` (legacy
  \\1/F^2\_{X_j}(s)\\, Shen 2015), or `"HC"` (experimental
  heteroscedasticity-consistent sandwich built from the per-box
  detrended moment scores; `overlap = FALSE` is recommended and the
  finite-sample normalisation is still being calibrated).

- abs:

  logical. If TRUE the absolute detrended covariance is used in the
  cross-correlation step (more robust to outliers).

## Value

A list with, among others: scale `s`, detrended fluctuation `F`, `DCCA`,
`DPCCA`, beta estimates `BDFA`, standardized betas `BSDFA`, residual
variance `UDFA`, coefficient variance `VDFA`, multiple correlation
`DMC2`, `R2DFA`, confidence limits `UCIB`/`LCIB`, `p.value`, critical
value `TC`, the scale-wise diagnostics `VIF`, condition number `kappa`
and adjusted `R2adj`, and the per-series DFA exponent `alpha` (a
long-memory proxy).

## Details

The variance of the scale-dependent coefficients follows Tilfani et al.
(2022, Eq. 25), \\\mathrm{var}(\hat\beta_j(s)) = F^2\_\varepsilon(s)
\[F\_{XX}(s)^{-1}\]\_{jj}\\, i.e. it uses the full inverse of the
detrended covariance matrix of the predictors, consistent with how the
coefficients themselves are estimated. Under collinearity this differs
from the legacy "marginal" form \\F^2\_\varepsilon(s) / F^2\_{X_j}(s)\\
(which under-covers); the two coincide for orthogonal predictors. Set
`vcov = "marginal"` to reproduce the legacy behaviour.

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
round(fit$BDFA[, 1, 10], 2)               # coefficients at the 10th scale
#> [1]  0.66 -0.34
```
