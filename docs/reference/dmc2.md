# rho Detrended Multiple-Correlation Coefficient

Calculates DMC²

## Usage

``` r
dmc2(data, dpo = 1, int = TRUE, np = 91, overlap = TRUE)
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

Scale s, Multiple Detrended Correlation

## References

Zebende, G. F. and da Silva Filho, A. M. (2018). Detrended multiple
cross-correlation coefficient. *Physica A*, 510, 91-97.

Wang, F., Xu, J. and Fan, Q. (2021). Statistical test for detrended
multiple cross-correlation coefficient. *Communications in Nonlinear
Science and Numerical Simulation*, 99, 105781.

## See also

[`vignette("DFATools")`](https://ikarobarreto.github.io/DFATools/articles/DFATools.md)

## Examples

``` r
set.seed(1)
d <- data.frame(y = cumsum(rnorm(300)), x1 = cumsum(rnorm(300)),
                x2 = cumsum(rnorm(300)))
dmc2(d, np = 20)
#> # A tibble: 20 × 2
#>       sn     dmc2
#>    <int>    <dbl>
#>  1    10  0.00668
#>  2    12  0.00943
#>  3    13  0.0115 
#>  4    14  0.0141 
#>  5    15  0.0175 
#>  6    16  0.0219 
#>  7    17  0.0273 
#>  8    18  0.0337 
#>  9    19  0.0413 
#> 10    20  0.0499 
#> 11    21  0.0595 
#> 12    23  0.0811 
#> 13    25  0.105  
#> 14    27  0.131  
#> 15    29  0.158  
#> 16    31  0.187  
#> 17    34  0.230  
#> 18    37  0.273  
#> 19    40  0.316  
#> 20    43 NA      
```
