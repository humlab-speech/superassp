# Per-GCI LPC analysis and inverse filtering R wrapper around the Rcpp implementation of `calc_residual.m`. If GCI vector is too short, falls back to pure R.

Per-GCI LPC analysis and inverse filtering R wrapper around the Rcpp
implementation of `calc_residual.m`. If GCI vector is too short, falls
back to pure R.

## Usage

``` r
calc_residual_r(x, x_lpc, ord_lpc, GCI)
```
