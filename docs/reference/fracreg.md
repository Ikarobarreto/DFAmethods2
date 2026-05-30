# Multiple Fractal Regression

Calculates Fractal Regression.

## Usage

``` r
fracreg(data, dpo, int, np = 91, overlap = T)
```

## Arguments

- data:

  is a matrix of time series.

- dpo:

  Detrending polynomial order.

- int:

  logical. if TRUE integration process will be applied.

- np:

  number of point scales.

- overlap:

  logical. if TRUE overlapping windows will be applied.

## Value

Scale s, Detrended Fluctuation Function F, Rho-DCCA, Rho DPCCA, Beta DFA
estimates, Standardized Beta DFA estimates, DFA Residuals, DFA Variance,
DFA Upper and Lower confidence interval, Multiple Detrended Correlation,
DFA R² , DFA p-value and DFA Calculated T statistics.

## References

Barreto, I. D. C., Dore, L. H., Stosic, T. and Stosic, B. D. (2021).
Extending DFA-based multiple linear regression inference: application to
acoustic impedance models. *Physica A*, 582, 126259.

Shen, C. (2015). A new detrended semipartial cross-correlation analysis.
*Physics Letters A*, 379(44), 2962-2969.

Kristoufek, L. (2015). Detrended fluctuation analysis as a regression
framework: estimating dependence at different scales. *Physical Review
E*, 91(2), 022802.

## See also

[`betadfa`](https://ikarobarreto.github.io/DFAmethods2/reference/betadfa.md),
[`effsizeDFA`](https://ikarobarreto.github.io/DFAmethods2/reference/effsizeDFA.md),
[`vignette("DFAmethods2")`](https://ikarobarreto.github.io/DFAmethods2/articles/DFAmethods2.md)
