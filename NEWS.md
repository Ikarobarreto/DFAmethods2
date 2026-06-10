# DFATools 1.0.0

First public release. `DFATools` is the continuation, under a new name, of the
development package previously called `DFAmethods2`; it collects Detrended
Fluctuation Analysis and related scale-dependent methods in one place.

## Methods

* Detrended fluctuation analysis (`dfa`) and detrended cross-correlation
  analysis, with the detrended cross-correlation and partial cross-correlation
  coefficients (`rhodcca`, `rhodpcca`) and the detrended multiple correlation
  (`dmc2`).
* Scale-dependent (fractal) multiple regression `fracreg()`, with
  scale-dependent standardized coefficients, the scale-wise effect size
  (`effsizeDFA`) and the smooth coefficient profile (`betadfa`, `sbdfa`).
* Significance tests for the scale-dependent coefficients: Shen-Podobnik
  (`fracreg.PStest`), Kristoufek (`fracreg.Ktest`) and the intersection-union
  test (`fracreg.IUTest`). The scale-dependent standardized coefficients, the
  scale-wise effect size and the intersection-union test follow Barreto et al.
  (2021) <doi:10.1016/j.physa.2021.126259>.

## Inference for `fracreg()`

* **Default coefficient variance is now the memory-corrected inverse form**
  (Barreto et al. 2026, Eq. M4.2'):
  `var(beta_j(s)) = F2_eps(s) * [F_XX(s)^-1]_jj / (T_s - k) * (2 H_hat + 1)^2 / 21`,
  where `H_hat` is the DFA exponent of the OLS residual (estimated internally,
  or supplied via `H_eps`). The factor restores nominal CI coverage under
  polynomial detrending; the uncorrected inverse over-covers and the legacy
  marginal form under-covers under collinearity. The estimator is selected with
  `variance = c("inv_corrected", "inv", "marginal", "hc", "none")` (this replaces
  the former `vcov` argument).
* **`overlap` now defaults to `FALSE`.** Score-based inference treats the boxes
  as disjoint sampling units; `overlap = TRUE` returns point estimates only
  (`variance = "none"`) with a message.
* `fracreg()` gains output fields `H_resid` (residual DFA exponent),
  `c_factor` (the memory multiplier `c(H) = (2H+1)^2 / 21`, the inverse of the
  bilinear ratio `kappa(s, H)` defined in Barreto et al. 2026; renamed from
  the earlier `kappa_factor`, which named the same quantity by its inverse),
  `variance_method`, `df_eff` (= `T_s - k`; the package uses `T_s - k` rather
  than the `T_s - k - 1` of Tilfani et al. 2022 because the scale-wise
  intercept is derived from the slopes and does not consume a degree of
  freedom)
  and `Ts` (= `floor(N / s)`). The adjusted R-squared now uses `T_s` (Eq. M3.1).
* The variance and the matching `t` quantile are normalised by the residual
  degrees of freedom `T_s - k`, where `T_s = floor(N / s)` is the number of
  non-overlapping boxes and `k` the number of predictors (the count of disjoint
  boxes, regardless of `overlap`). Scales with `T_s <= k` return `NA` limits
  with a warning.
* `fracreg()` now returns the per-series DFA exponent `$alpha` and warns when it
  exceeds 3/4 in any series, where the analytic interval can under-cover under
  strong long-range dependence (Hermite-Rosenblatt threshold) and the dependent
  wild bootstrap (`fracreg.WB()`) is preferred.
* `fracreg()` now reports **two-sided** p-values (`$p.value`) for the
  scale-dependent coefficients, `2 * (1 - pt(|t|, T_s - k))`.
* New `min_boxes` argument (default 15) in `fracreg()` and `fracreg.WB()`:
  scales whose non-overlapping box count `T_s = floor(N/s)` falls below this
  floor return `NA` standard errors, confidence limits and p-values, listed in
  a single warning. A separate warning is issued when `N < 500`
  (Likens et al. 2019).
* `fracreg.WB()` no longer accepts an `overlap` argument; score-based inference
  requires disjoint boxes, so the function is non-overlapping by construction
  (paper M8). Calls that pass `overlap = ...` will fail with the standard
  "unused argument" error.
* New `variance = "inv_theoretical"` and `auto_select_kappa` argument.
  The package now ships an internal lookup table `kappa_th_table`
  (Barreto et al. 2026 Table A.1) of the theoretical bilinear ratio
  `kappa_th(s, H)`. The new method `variance = "inv_theoretical"` applies
  the per-scale factor `1 / kappa_th(s, H_resid)` bilinearly interpolated
  from this table -- accurate over the full calibrated H regime, not just
  the closed-form's `[0.5, 0.95]` range. `auto_select_kappa = TRUE` (the
  default) makes `variance = "inv_corrected"` use the table only at scales
  with `T_s > 500`, where the closed form's absorbed `1/T_s` term has
  become negligible and the constant `21` would over-shrink; set FALSE to
  force the closed form everywhere. The returned `$c_factor` field is now
  a per-scale vector reflecting the factor actually applied.
* `fracreg()` reserves a new `variance = "wildboot"` placeholder for an
  analytical dependent-wild-bootstrap variant planned for DFATools >= 1.1
  (cf. Barreto et al. 2026 P2 programme); calling it today errors with a
  pointer to `fracreg.WB(weights = 'dependent')`, the implemented wild
  bootstrap built on the same per-box moment scores.
* `fracreg()` now issues a graded advisory on the residual exponent: a
  message when `H_resid < 0.5` (below the calibration range of the closed-
  form factor); the existing Hermite-Rosenblatt warning at `H_resid > 0.75`;
  a message when `H_resid >= 0.95` (close to the non-stationary regime);
  and a stronger saturation warning at `H_resid >= 0.99` flagging possible
  non-stationary regime or omitted long-memory variable.
* Reference list extended with Hu et al. 2001 (PRE 64:011114), Kantelhardt et
  al. 2002 (Physica A 316:87), Kwapien et al. 2015 (PRE 92:052815), Sikora et
  al. 2020 (PRE 101:032114) and Cavalcanti 2019 (PhD thesis, UFRPE).
* **`dfa()` now follows the Peng convention** explicitly. `$F` is the
  root-mean-squared fluctuation \eqn{F(s) = \sqrt{F^2(s)}}, so the DFA exponent
  is the slope of `log(F)` against `log(s)` directly. The previous output
  (the mean *squared* fluctuation, exposed as `$F`) is now in `$F2`
  (`sqrt($F2) == $F` by construction). The estimated exponent itself is
  returned in the new `$alpha` field. The return shape changed from a tibble
  to a list. Code that took `slope(log(dfa()$F))` and labelled it
  \eqn{\alpha} was reporting `2 * alpha`; that code now gives the correct
  Peng exponent. Plot helpers (`plotdfa()`) and the internal residual-
  exponent estimator in `fracreg()` have been updated accordingly.
* `vcov = "HC"` is an experimental heteroscedasticity-consistent (sandwich)
  estimator built from per-box detrended-moment scores; `overlap = FALSE` is
  recommended.
* Scale-wise collinearity diagnostics are returned: the variance inflation
  factors `$VIF`, the condition number `$kappa` of the predictors' scale-wise
  correlation matrix, and the adjusted coefficient of determination `$R2adj`
  (Tilfani et al. 2022, Eq. 26). `fracreg.diag()` returns them as a tidy tibble.
* `fracreg.WB()` gives wild-bootstrap confidence intervals and p-values for the
  scale-dependent coefficients, robust to heteroscedasticity and (with the
  default `"dependent"` weights, Shao 2010) to dependence between boxes. It
  resamples the per-box moment scores, so it does not recompute the DFA per
  replicate; the i.i.d. weights reproduce the `vcov = "HC"` sandwich variance.
* `abs` uses the absolute detrended covariance in the cross-correlation step
  (more robust to outliers).
