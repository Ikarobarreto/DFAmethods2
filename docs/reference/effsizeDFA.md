# f² scale-wise effect sizes

Calculates f² scale-wise effect sizes

## Usage

``` r
effsizeDFA(data, dpo = 1, int = TRUE, np = 91, overlap = TRUE)
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

[`fracreg`](https://ikarobarreto.github.io/DFATools/reference/fracreg.md),
[`vignette("DFATools")`](https://ikarobarreto.github.io/DFATools/articles/DFATools.md)

## Examples

``` r
set.seed(1)
d <- data.frame(y = cumsum(rnorm(300)), x1 = cumsum(rnorm(300)),
                x2 = cumsum(rnorm(300)))
effsizeDFA(d, np = 20)
#> # A tibble: 20 × 6
#>        s         R2      R2_x1       R2_x2      f2_x1      f2_x2
#>    <int>      <dbl>      <dbl>       <dbl>      <dbl>      <dbl>
#>  1    10  0.0000446  0.0000150  0.00000341  0.0000296  0.0000412
#>  2    12  0.0000889  0.0000149  0.0000143   0.0000740  0.0000746
#>  3    13  0.000132   0.0000175  0.0000241   0.000114   0.000108 
#>  4    14  0.000200   0.0000221  0.0000389   0.000178   0.000161 
#>  5    15  0.000307   0.0000288  0.0000630   0.000278   0.000244 
#>  6    16  0.000478   0.0000387  0.000102    0.000440   0.000376 
#>  7    17  0.000743   0.0000525  0.000162    0.000691   0.000581 
#>  8    18  0.00114    0.0000713  0.000253    0.00107    0.000887 
#>  9    19  0.00171    0.0000962  0.000381    0.00161    0.00133  
#> 10    20  0.00249    0.000128   0.000558    0.00237    0.00194  
#> 11    21  0.00354    0.000168   0.000794    0.00339    0.00276  
#> 12    23  0.00658    0.000274   0.00147     0.00635    0.00514  
#> 13    25  0.0110     0.000415   0.00244     0.0107     0.00867  
#> 14    27  0.0171     0.000604   0.00374     0.0168     0.0136   
#> 15    29  0.0250     0.000845   0.00542     0.0248     0.0201   
#> 16    31  0.0348     0.00115    0.00749     0.0349     0.0283   
#> 17    34  0.0529     0.00176    0.0113      0.0540     0.0440   
#> 18    37  0.0747     0.00256    0.0160      0.0780     0.0635   
#> 19    40  0.0997     0.00354    0.0217      0.107      0.0866   
#> 20    43 NA         NA         NA          NA         NA        
```
