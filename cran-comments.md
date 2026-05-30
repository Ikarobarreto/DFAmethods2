## Resubmission

This is a resubmission. The previous submission (0.1.0) passed the incoming
checks on Windows but failed on 64-bit Linux. In this version I have:

* Fixed a memory error (segmentation fault on 64-bit Linux). The box-size array
  was passed from R via `.C()` as an integer vector but declared as `long *` in
  C; on LP64 platforms `long` is 8 bytes while an R integer is 4 bytes, which
  corrupted memory. The C code now uses `int *`, matching `as.integer()`.
* Removed dead C code that wrote to `stderr` (the "compiled code ... 'stderr'"
  note from the Debian check).

## Test environments
* local Windows 10 install, R 4.5.1
* win-builder (R-devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* New submission.

The note also flags possibly misspelled words in DESCRIPTION (Barreto,
Kristoufek, Podobnik, Shen, "et", "al"); these are author/researcher surnames
and the standard "et al." citation abbreviation, and are spelled correctly.
