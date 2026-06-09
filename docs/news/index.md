# Changelog

## DFATools 1.0.0

First public release. `DFATools` is the continuation, under a new name,
of the development package previously called `DFAmethods2`; it collects
Detrended Fluctuation Analysis and related scale-dependent methods in
one place.

### Methods

- Detrended fluctuation analysis (`dfa`) and detrended cross-correlation
  analysis, with the detrended cross-correlation and partial
  cross-correlation coefficients (`rhodcca`, `rhodpcca`) and the
  detrended multiple correlation (`dmc2`).
- Scale-dependent (fractal) multiple regression
  [`fracreg()`](https://ikarobarreto.github.io/DFATools/reference/fracreg.md),
  with scale-dependent standardized coefficients, the scale-wise effect
  size (`effsizeDFA`) and the smooth coefficient profile (`betadfa`,
  `sbdfa`).
- Significance tests for the scale-dependent coefficients: Shen-Podobnik
  (`fracreg.PStest`), Kristoufek (`fracreg.Ktest`) and the
  intersection-union test (`fracreg.IUTest`). The scale-dependent
  standardized coefficients, the scale-wise effect size and the
  intersection-union test follow Barreto et al.
  2021. <doi:10.1016/j.physa.2021.126259>.

### Inference for `fracreg()`

- **Default coefficient variance is now the memory-corrected inverse
  form** (Barreto et al. 2026, Eq. M4.2’):
  `var(beta_j(s)) = F2_eps(s) * [F_XX(s)^-1]_jj / (T_s - k) * (2 H_hat + 1)^2 / 21`,
  where `H_hat` is the DFA exponent of the OLS residual (estimated
  internally, or supplied via `H_eps`). The factor restores nominal CI
  coverage under polynomial detrending; the uncorrected inverse
  over-covers and the legacy marginal form under-covers under
  collinearity. The estimator is selected with
  `variance = c("inv_corrected", "inv", "marginal", "hc", "none")` (this
  replaces the former `vcov` argument).
- **`overlap` now defaults to `FALSE`.** Score-based inference treats
  the boxes as disjoint sampling units; `overlap = TRUE` returns point
  estimates only (`variance = "none"`) with a message.
- [`fracreg()`](https://ikarobarreto.github.io/DFATools/reference/fracreg.md)
  gains output fields `H_resid` (residual DFA exponent), `kappa_factor`
  (the memory factor), `variance_method`, `df_eff` (= `T_s - k`) and
  `Ts` (= `floor(N / s)`). The adjusted R-squared now uses `T_s` (Eq.
  M3.1).
- The variance and the matching `t` quantile are normalised by the
  residual degrees of freedom `T_s - k`, where `T_s = floor(N / s)` is
  the number of non-overlapping boxes and `k` the number of predictors
  (the count of disjoint boxes, regardless of `overlap`). Scales with
  `T_s <= k` return `NA` limits with a warning.
- [`fracreg()`](https://ikarobarreto.github.io/DFATools/reference/fracreg.md)
  now returns the per-series DFA exponent `$alpha` and warns when it
  exceeds 3/4 in any series, where the analytic interval can under-cover
  under strong long-range dependence (Hermite-Rosenblatt threshold) and
  the dependent wild bootstrap
  ([`fracreg.WB()`](https://ikarobarreto.github.io/DFATools/reference/fracreg.WB.md))
  is preferred.
- [`fracreg()`](https://ikarobarreto.github.io/DFATools/reference/fracreg.md)
  now reports **two-sided** p-values (`$p.value`) for the
  scale-dependent coefficients, `2 * (1 - pt(|t|, T_s - k))`.
- `vcov = "HC"` is an experimental heteroscedasticity-consistent
  (sandwich) estimator built from per-box detrended-moment scores;
  `overlap = FALSE` is recommended.
- Scale-wise collinearity diagnostics are returned: the variance
  inflation factors `$VIF`, the condition number `$kappa` of the
  predictors’ scale-wise correlation matrix, and the adjusted
  coefficient of determination `$R2adj` (Tilfani et al. 2022, Eq. 26).
  [`fracreg.diag()`](https://ikarobarreto.github.io/DFATools/reference/fracreg.diag.md)
  returns them as a tidy tibble.
- [`fracreg.WB()`](https://ikarobarreto.github.io/DFATools/reference/fracreg.WB.md)
  gives wild-bootstrap confidence intervals and p-values for the
  scale-dependent coefficients, robust to heteroscedasticity and (with
  the default `"dependent"` weights, Shao 2010) to dependence between
  boxes. It resamples the per-box moment scores, so it does not
  recompute the DFA per replicate; the i.i.d. weights reproduce the
  `vcov = "HC"` sandwich variance.
- `abs` uses the absolute detrended covariance in the cross-correlation
  step (more robust to outliers).
