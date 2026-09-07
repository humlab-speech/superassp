# Autocorrelation LPC analysis Computes LPC coefficients and residual energy for a windowed frame using the autocorrelation method (Levinson-Durbin).

Autocorrelation LPC analysis Computes LPC coefficients and residual
energy for a windowed frame using the autocorrelation method
(Levinson-Durbin).

## Usage

``` r
lpcauto_frame(s, p)
```

## Arguments

- s:

  Signal frame (numeric vector), should already be windowed

- p:

  LPC order

## Value

List with `ar` (LPC coefficients, first = 1) and `e` (residual energy)
