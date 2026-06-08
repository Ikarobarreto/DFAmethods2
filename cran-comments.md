## Submission

This is a new submission.

`DFATools` collects Detrended Fluctuation Analysis and related scale-dependent
methods (detrended cross-correlation and partial cross-correlation coefficients,
detrended multiple regression with scale-wise inference, and the associated
significance tests). It is the continuation, under a new name, of a package
previously developed under the name `DFAmethods2` (never released on CRAN).

## Test environments
* local Windows 10 install, R 4.5.1
* win-builder (R-devel and R-release)
* GitHub Actions: Ubuntu (R-devel/release/oldrel), macOS, Windows

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.
* The note flags possibly misspelled words in DESCRIPTION (author/researcher
  surnames and the standard "et al." citation abbreviation), which are spelled
  correctly.

The package contains compiled C code; the box-size arrays are passed between R
and C as `int` and the package registration uses `R_init_DFATools`, both checked
on 64-bit Linux via GitHub Actions.
