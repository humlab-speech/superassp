# Column name behind an aesthetic mapping

Column name behind an aesthetic mapping

## Usage

``` r
.assp_aes_name(aes)
```

## Arguments

- aes:

  Aesthetic from a
  [`ggplot2::aes()`](https://ggplot2.tidyverse.org/reference/aes.html)
  mapping (a quosure), or `NULL`.

## Value

Character. The column name when the aesthetic is a single column, `NULL`
for anything else (constants, calls like `log(F1_Hz)`), where ggplot2's
own label is the right one.
