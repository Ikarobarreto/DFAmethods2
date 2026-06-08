# Plot DFA

Plot of Detrended Fluctuation Analysis

## Usage

``` r
plotdfa(dfa, seg = FALSE, point = NULL, main = NULL)
```

## Arguments

- dfa:

  is a dfa object

- seg:

  logical. If TRUE, alpha DFA will be calculated in 2 segments

- point:

  indicate in which point segmented alpha DFA should be calculated

- main:

  plot title

## Value

a plot of Detrended Fluctuation Analysis.

## Examples

``` r
set.seed(1)
x <- cumsum(rnorm(300))
plotdfa(dfa(x, np = 20))
#> `geom_smooth()` using formula = 'y ~ x'
```
