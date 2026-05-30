# f² scale-wise effect sizes

Calculates f² scale-wise effect sizes

## Usage

``` r
effsizeDFA(data, dpo = 1, int = T, np = 91, overlap = T)
```

## Arguments

- data:

  is a matrix of time series

- dpo:

  detrending polynomial order

- int:

  logical. if TRUE integration process will be applied.

- np:

  number of point scales.

- overlap:

  logical. if TRUE overlapping windows will be applied.

## Value

Scale s, f² scale-wise effect sizes

## References

Barreto, I. D. C., Dore, L. H., Stosic, T. and Stosic, B. D. (2021).
Extending DFA-based multiple linear regression inference: application to
acoustic impedance models. *Physica A*, 582, 126259.

Cohen, J. (1988). *Statistical Power Analysis for the Behavioral
Sciences*, 2nd ed. Lawrence Erlbaum Associates.

## See also

[`fracreg`](https://ikarobarreto.github.io/DFAmethods2/reference/fracreg.md),
[`vignette("DFAmethods2")`](https://ikarobarreto.github.io/DFAmethods2/articles/DFAmethods2.md)
