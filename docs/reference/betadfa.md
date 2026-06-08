# Beta DFA

Calculates Beta DFA

## Usage

``` r
betadfa(data, dpo = 1, int = TRUE, np = 91, overlap = TRUE)
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

Scale s, Beta DFA estimates

## References

Kristoufek, L. (2015). Detrended fluctuation analysis as a regression
framework: estimating dependence at different scales. *Physical Review
E*, 91(2), 022802.

## See also

[`sbdfa`](https://ikarobarreto.github.io/DFATools/reference/sbdfa.md),
[`fracreg`](https://ikarobarreto.github.io/DFATools/reference/fracreg.md),
[`vignette("DFATools")`](https://ikarobarreto.github.io/DFATools/articles/DFATools.md)

## Examples

``` r
set.seed(1)
d <- data.frame(y = cumsum(rnorm(300)), x1 = cumsum(rnorm(300)),
                x2 = cumsum(rnorm(300)))
betadfa(d, np = 20)
#> # A tibble: 20 × 3
#>        s     x1      x2
#>    <int>  <dbl>   <dbl>
#>  1    10 0.0479 -0.0690
#>  2    12 0.0687 -0.0774
#>  3    13 0.0795 -0.0849
#>  4    14 0.0914 -0.0943
#>  5    15 0.105  -0.105 
#>  6    16 0.121  -0.118 
#>  7    17 0.139  -0.132 
#>  8    18 0.158  -0.147 
#>  9    19 0.179  -0.163 
#> 10    20 0.201  -0.180 
#> 11    21 0.224  -0.197 
#> 12    23 0.271  -0.231 
#> 13    25 0.317  -0.264 
#> 14    27 0.362  -0.296 
#> 15    29 0.407  -0.325 
#> 16    31 0.450  -0.352 
#> 17    34 0.510  -0.389 
#> 18    37 0.566  -0.419 
#> 19    40 0.617  -0.444 
#> 20    43 0.663  -0.463 
```
