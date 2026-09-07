# Convert ERB-rate Scale to Frequency

Converts values from the ERB-rate scale to frequency in Hz. This is the
inverse of
[`ucnv_hz_to_erb`](https://humlab-speech.github.io/superassp/reference/ucnv_hz_to_erb.md).

## Usage

``` r
ucnv_erb_to_hz(erb, method = c("glasberg1990", "moore1983"), as_units = NULL)
```

## Arguments

- erb:

  Numeric vector or units object with ERB-rate values.

- method:

  Character string specifying which inverse formula to use. Must match
  the forward conversion method.

- as_units:

  Logical. If TRUE, returns a units object with "Hz" units.

## Value

If `as_units = TRUE`: a units object with "Hz" units. Otherwise, a
numeric vector of frequencies in Hz.

## See also

[`ucnv_hz_to_erb`](https://humlab-speech.github.io/superassp/reference/ucnv_hz_to_erb.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Round trip
erb <- ucnv_hz_to_erb(1000)
ucnv_erb_to_hz(erb)  # Should return ~1000
} # }
```
