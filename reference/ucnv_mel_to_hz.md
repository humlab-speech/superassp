# Convert Mel Scale to Frequency

Converts values from the mel scale to frequency in Hz. This is the
inverse of
[`ucnv_hz_to_mel`](https://humlab-speech.github.io/superassp/reference/ucnv_hz_to_mel.md).

## Usage

``` r
ucnv_mel_to_hz(mel, method = c("htk", "slaney"), as_units = NULL)
```

## Arguments

- mel:

  Numeric vector or units object with mel values.

- method:

  Character string specifying which inverse formula to use.

- as_units:

  Logical. If TRUE, returns a units object with "Hz" units.

## Value

If `as_units = TRUE`: a units object with "Hz" units. Otherwise, a
numeric vector of frequencies in Hz.

## See also

[`ucnv_hz_to_mel`](https://humlab-speech.github.io/superassp/reference/ucnv_hz_to_mel.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Round trip
mel <- ucnv_hz_to_mel(1000)
ucnv_mel_to_hz(mel)  # Should return ~1000
} # }
```
