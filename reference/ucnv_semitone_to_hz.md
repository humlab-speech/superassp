# Convert Semitones to Frequency

Converts semitone values to frequency in Hz relative to a reference
frequency. This is the inverse of
[`ucnv_hz_to_semitone`](https://humlab-speech.github.io/superassp/reference/ucnv_hz_to_semitone.md).

## Usage

``` r
ucnv_semitone_to_hz(
  semitone,
  ref_freq = NULL,
  ref_source = c("A4", "UEP83", "Praat"),
  as_units = NULL
)
```

## Arguments

- semitone:

  Numeric vector or units object with semitone values.

- ref_freq:

  Reference frequency in Hz. Default is NULL, which uses the reference
  determined by `ref_source`. Can also be a units object. If specified,
  overrides `ref_source`.

- ref_source:

  Character string specifying the reference standard. Options are:

  - `"UEP83"`: UEP 1983 standard (110 Hz = A2)

  - `"Praat"`: Praat convention (100 Hz)

  - `"A4"` (default): Concert pitch A4 = 440 Hz

  Ignored if `ref_freq` is explicitly specified.

- as_units:

  Logical. If TRUE, returns a units object with "Hz" units.

## Value

If `as_units = TRUE`: a units object with "Hz" units. Otherwise, a
numeric vector of frequencies in Hz.

## Details

The conversion formula is: \$\$f = f\_{ref} \cdot 2^{ST/12}\$\$

See
[`ucnv_hz_to_semitone`](https://humlab-speech.github.io/superassp/reference/ucnv_hz_to_semitone.md)
for details on reference standards.

## See also

[`ucnv_hz_to_semitone`](https://humlab-speech.github.io/superassp/reference/ucnv_hz_to_semitone.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Using different reference standards
ucnv_semitone_to_hz(12, ref_source = "UEP83")   # 220 Hz (octave above 110)
ucnv_semitone_to_hz(12, ref_source = "Praat")   # 200 Hz (octave above 100)
ucnv_semitone_to_hz(12, ref_source = "A4")      # 880 Hz (octave above 440)

# Musical scale from A4
ucnv_semitone_to_hz(0:12, ref_freq = 440)
} # }
```
