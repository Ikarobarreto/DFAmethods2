# Plot rho-DCCA

Plot of Detrended Cross Correlation Coefficient

## Usage

``` r
plotrdcca(rdcca, var)
```

## Arguments

- rdcca:

  is a rhodcca object

- var:

  character. Indicate which pair in rho dcca object you want to plot.

## Value

a plot of Detrended Cross Correlation Analysis.

## Examples

``` r
set.seed(1)
d <- data.frame(x = cumsum(rnorm(300)), y = cumsum(rnorm(300)))
plotrdcca(rhodcca(d, np = 20), var = "12")
```
