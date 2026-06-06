## Update

This is a minor update (0.1.3).

* The coefficient-variance estimator in `fracreg()` now uses the full inverse of
  the predictors' detrended covariance matrix (Tilfani et al. 2022, Eq. 25),
  consistent with the coefficient estimates. This changes the confidence limits,
  p-values and critical values for correlated predictors; the previous behaviour
  is available via `vcov = "marginal"`. Documented in NEWS.md.
* Adds scale-wise collinearity diagnostics (`$VIF`, `$kappa`, `$R2adj`) to
  `fracreg()` and a new `fracreg.diag()` function.

## Test environments
* local Windows 10 install, R 4.5.1
* win-builder (R-devel)
* GitHub Actions: Ubuntu (R-devel/release/oldrel), macOS, Windows

## R CMD check results

0 errors | 0 warnings | 1 note

The note flags possibly misspelled words in DESCRIPTION (author/researcher
surnames and the standard "et al." citation abbreviation), which are spelled
correctly.
