# Intersection-Union Test

Calculates Intersection-Union Test for Beta-DFA = 0 or Beta-DFA =
Beta-OLS

## Usage

``` r
fracreg.IUTest(data, B = 100, dpo = 1, int = T, np = 91, overlap = T)
```

## Arguments

- data:

  is a matrix of time series

- B:

  number of surrogate series

- dpo:

  detrending polynomial order

- int:

  logical. if TRUE integration process will be applied.

- np:

  number of point scales.

- overlap:

  logical. if TRUE overlapping windows will be applied.

## Value

A matrix with scale-wise Beta-DFA and critic region of Kristoufek Test
and Podobnik-Shen Test.

## References

Barreto, I. D. C., Dore, L. H., Stosic, T. and Stosic, B. D. (2021).
Extending DFA-based multiple linear regression inference: application to
acoustic impedance models. *Physica A*, 582, 126259.

## See also

[`fracreg.PStest`](https://ikarobarreto.github.io/DFAmethods2/reference/fracreg.PStest.md),
[`fracreg.Ktest`](https://ikarobarreto.github.io/DFAmethods2/reference/fracreg.Ktest.md),
[`vignette("DFAmethods2")`](https://ikarobarreto.github.io/DFAmethods2/articles/DFAmethods2.md)
