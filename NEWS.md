# DFAmethods2 0.1.3

## Behaviour change

* `fracreg()` now estimates the coefficient variance from the **full inverse**
  of the detrended covariance matrix of the predictors,
  `var(beta_j(s)) = F2_eps(s) * [F_XX(s)^-1]_jj` (Tilfani et al. 2022, Eq. 25),
  consistent with how the coefficients themselves are estimated. The previous
  "marginal" form, `F2_eps(s) / F2_Xj(s)`, under-covers under multicollinearity;
  the two coincide for orthogonal predictors. Set `vcov = "marginal"` to restore
  the previous behaviour. This changes the confidence limits (`UCIB`, `LCIB`),
  `p.value` and `TC` when the predictors are correlated.

## Bug fixes

* Fixed a numerical bug in the residual detrended fluctuation (`UDFA`) and hence
  in `R2DFA` and `VDFA2`: `polyfit()`'s incremental `covar0` cache key
  (`pboxsize`) was `static` and persisted across the per-scale `setup()` calls,
  which reallocate `covar0` without re-zeroing it, so the residual DFA in
  `fracreg()` accumulated onto uninitialised memory and returned nonsensical
  (even negative) values. `pboxsize` is now reset by `setup()`/`setup_dcca()`.
  The main confidence intervals were not affected (they use F2_eps = F2_Y(1-DMC)).

## New features

* `fracreg()` gains a `vcov = c("inverse", "marginal")` argument (default
  `"inverse"`) and an `abs` argument to use the absolute detrended covariance in
  the cross-correlation step (more robust to outliers).
* `fracreg()` now also returns scale-wise collinearity diagnostics: the variance
  inflation factors `$VIF`, the condition number `$kappa` of the predictors'
  scale-wise correlation matrix, and the adjusted coefficient of determination
  `$R2adj` (Tilfani et al. 2022, Eq. 26).
* New `fracreg.diag()` returns those scale-wise diagnostics as a tidy tibble.

# DFAmethods2 0.1.2

* Addressed CRAN review: use `TRUE`/`FALSE` instead of `T`/`F`; added executable
  examples to all Rd files; replaced `print()`/`cat()` console output with
  `message()`/`warning()`/`stop()`.

# DFAmethods2 0.1.1

* Fixed a segmentation fault on 64-bit Linux: the box-size array was passed from
  R as an integer vector but declared `long *` in C.

# DFAmethods2 0.1.0

* First release.
