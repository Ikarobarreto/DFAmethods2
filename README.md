
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DFAmethods2

<!-- badges: start -->

<!-- badges: end -->

**DFAmethods2** is a toolbox of Detrended Fluctuation Analysis (DFA)
methods for measuring long-range correlations, cross-correlations and
regression in **nonstationary** time series, where classical
autocorrelation and ordinary regression are biased by trends. It
provides:

- **Detrended Fluctuation Analysis** — `dfa()`
- **Detrended Cross-Correlation Analysis** and the **rho-DCCA**
  coefficient — `rhodcca()`
- **Detrended Partial Cross-Correlation** coefficient — `rhodpcca()`
- **Detrended Multiple Cross-Correlation** coefficient (DMC) — `dmc2()`
- **Detrended fractal regression**: scale-wise and standardized
  coefficients — `betadfa()`, `sbdfa()`, `fracreg()`
- **Scale-wise** $f^2$ **effect sizes** — `effsizeDFA()`
- **Significance tests**: Podobnik-Shen, Kristoufek and
  intersection-union — `fracreg.PStest()`, `fracreg.Ktest()`,
  `fracreg.IUTest()`

The scale-dependent standardized coefficients, the scale-wise effect
size and the intersection-union test follow Barreto et al. (2021)
[doi:10.1016/j.physa.2021.126259](https://doi.org/10.1016/j.physa.2021.126259).

## Installation

``` r
# install.packages("devtools")
devtools::install_github("Ikarobarreto/DFAmethods2")
```

## Example

A response `y` driven by a covariate `x2`, both random-walk-like:

``` r
library(DFAmethods2)
#> Registered S3 method overwritten by 'quantmod':
#>   method            from
#>   as.zoo.data.frame zoo

set.seed(1)
x2  <- cumsum(rnorm(800))
y   <- 0.7 * x2 + cumsum(rnorm(800))
dat <- data.frame(y = y, x2 = x2)

# DFA fluctuation function of y
head(dfa(y, np = 30))
#> # A tibble: 6 × 2
#>       s     F
#>   <int> <dbl>
#> 1    10  3.50
#> 2    11  4.65
#> 3    12  6.03
#> 4    13  7.68
#> 5    14  9.62
#> 6    15 11.9

# scale-wise rho-DCCA cross-correlation coefficient
head(rhodcca(dat, np = 30))
#> # A tibble: 6 × 2
#>       s DCCA12
#>   <int>  <dbl>
#> 1    10  0.487
#> 2    11  0.471
#> 3    12  0.457
#> 4    13  0.445
#> 5    14  0.434
#> 6    15  0.426
```

See `vignette("DFAmethods2")` for the full theory and worked examples.
