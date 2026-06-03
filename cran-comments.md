## Resubmission

This is a resubmission addressing the reviewer comments (Benjamin Altmann) on
version 0.1.1:

* Replaced all uses of `T`/`F` with `TRUE`/`FALSE` in function arguments and
  bodies, and renamed the `F` variable used inside `plotdfa()`/`plotdcca()`.
* Added small executable examples to the Rd files of all exported functions, to
  illustrate their use and enable automatic testing.
* Replaced `print()`/`cat()` console messages with `warning()`/`stop()` in
  `plotdcca()` and `plotrdcca()`.

## Test environments
* local Windows 10 install, R 4.5.1
* win-builder (R-devel)
* GitHub Actions: Ubuntu (R-devel/release/oldrel), macOS, Windows

## R CMD check results

0 errors | 0 warnings | 1 note

* New submission.

The note also flags possibly misspelled words in DESCRIPTION (Barreto,
Kristoufek, Podobnik, Shen, "et", "al"); these are author/researcher surnames
and the standard "et al." citation abbreviation, and are spelled correctly.
