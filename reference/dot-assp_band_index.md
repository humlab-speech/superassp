# Split expanded column names back into coefficient bands

[`as.data.frame.AsspDataObj()`](https://humlab-speech.github.io/superassp/reference/AsspDataObj.md)
expands a matrix track into one column per coefficient (`DFT_dB_1`,
`DFT_dB_2`, ...). Grouping those names back lets a spectral track keep a
single band identity with a coefficient index.

## Usage

``` r
.assp_band_index(cols)
```

## Arguments

- cols:

  Character vector. Column names of the track table.

## Value

list with `band` (character) and `bin` (integer, `NA` outside a band).
