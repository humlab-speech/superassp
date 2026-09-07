# Convert Frequency to Mel Scale

Converts frequency values from Hz to the mel scale, a perceptual pitch
scale where equal distances sound equally different to human listeners.

## Usage

``` r
ucnv_hz_to_mel(freq, method = c("htk", "slaney"), as_units = NULL)
```

## Arguments

- freq:

  Numeric vector or units object with frequency values in Hz.

- method:

  Character string specifying the conversion formula to use. Options
  are:

  - `"htk"` (default): HTK formula using base-10 logarithm

  - `"slaney"`: Slaney's formula (linear below 1000 Hz)

- as_units:

  Logical. If TRUE, returns a units object with "mel" units.

## Value

If `as_units = TRUE`: a units object with "mel" units. Otherwise, a
numeric vector of mel values.

## Details

Two formulas are available:

**HTK formula** (default): \$\$mel = 2595 \cdot \log\_{10}(1 +
f/700)\$\$

This gives exactly 1000 mels at 1000 Hz.

**Slaney formula**: Linear below 1000 Hz, logarithmic above. Used in
some audio processing libraries.

## References

(O'Shaughnessy 1987)

## Examples

``` r
if (FALSE) { # \dontrun{
# 1000 Hz = 1000 mels (approximately)
ucnv_hz_to_mel(1000)

# Vector conversion
ucnv_hz_to_mel(c(100, 500, 1000, 2000, 4000))

# With units
library(units)
freq <- set_units(1000, Hz)
ucnv_hz_to_mel(freq)
} # }
```
