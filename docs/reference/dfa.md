# Detrended Fluctuation Analysis

Calculates the detrended fluctuation function of a single series.

## Usage

``` r
dfa(data, dpo = 1, int = TRUE, np = 91, overlap = TRUE)
```

## Arguments

- data:

  a numeric vector or single-column matrix.

- dpo:

  detrending polynomial order (default 1).

- int:

  logical; if TRUE the input is integrated into the profile (the
  standard use for stationary inputs).

- np:

  number of scales (box sizes).

- overlap:

  logical; if TRUE overlapping windows are used.

## Value

A list with the scale vector `s`, the fluctuation function `F`, the
squared fluctuation `F2` and the estimated DFA exponent `alpha`.

## Details

The C primitive computes the mean squared fluctuation \\F^2(s)\\ (the
average of the within-box detrended residual variance over the boxes of
size \\s\\). The return uses the conventional Peng et al. (1994) form:

- `$F = sqrt(F^2)` – the root mean-squared fluctuation, so that
  \\\alpha\\ is the slope of \\\log F(s)\\ vs \\\log s\\;

- `$F2 = F^2` – the squared fluctuation, the legacy quantity consumed
  internally by the package (`rhodcca`, `fracreg`, ...) and useful for
  combining DFA values across series;

- `$alpha` – the estimated DFA exponent (= Hurst exponent for
  self-similar processes), the slope of \\\log F^2(s)/2\\ against \\\log
  s\\ over all positive scales.

## References

Peng, C.-K., Buldyrev, S. V., Havlin, S., Simons, M., Stanley, H. E. and
Goldberger, A. L. (1994). Mosaic organization of DNA nucleotides.
*Physical Review E*, 49(2), 1685-1689.

## See also

[`plotdfa`](https://ikarobarreto.github.io/DFATools/reference/plotdfa.md),
[`vignette("DFATools")`](https://ikarobarreto.github.io/DFATools/articles/DFATools.md)

## Examples

``` r
set.seed(1)
x <- cumsum(rnorm(300))           # random walk: alpha ~ 1.5
fy <- dfa(x, np = 20)
round(fy$alpha, 3)
#> [1] 1.546
```
