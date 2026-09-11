# Convert Frequency to Semitones

Converts frequency values to semitones relative to a reference
frequency. In equal temperament, each semitone represents a frequency
ratio of 2^(1/12).

## Usage

``` r
ucnv_hz_to_semitone(
  freq,
  ref_freq = NULL,
  ref_source = c("A4", "UEP83", "Praat"),
  as_units = NULL
)
```

## Arguments

- freq:

  Numeric vector or units object with frequency values in Hz.

- ref_freq:

  Reference frequency in Hz. Default is NULL, which uses the reference
  determined by `ref_source`. Can also be a units object. If specified,
  overrides `ref_source`.

- ref_source:

  Character string specifying the reference standard. Options are:

  - `"UEP83"`: UEP 1983 standard (110 Hz = A2/A1 in Helmholtz notation)

  - `"Praat"`: Praat convention (100 Hz arbitrary reference)

  - `"A4"` (default): Concert pitch A4 = 440 Hz

  Ignored if `ref_freq` is explicitly specified.

- as_units:

  Logical. If TRUE, returns a units object with "semitone" units.

## Value

If `as_units = TRUE`: a units object with "semitone" units. Otherwise, a
numeric vector of semitone values.

## Details

The conversion formula is: \$\$ST = 12 \cdot \log_2(f / f\_{ref})\$\$

where f is the frequency and f_ref is the reference frequency.

**Reference Standards:**

- **UEP 1983** (Schutte & Seidner): Uses 110 Hz (A2, or A1 in Helmholtz
  notation) as the reference for voice range profiles (phonetograms).
  This standard is commonly used in clinical phoniatrics and voice
  assessment.

- **Praat**: Uses 100 Hz as an arbitrary reference frequency for
  semitone-based pitch analysis. This provides a convenient round number
  for relative pitch measurements in speech analysis.

- **A4** (default): Uses 440 Hz (A4, concert pitch) as reference,
  following standard musical tuning conventions.

This gives the number of semitones above (positive) or below (negative)
the reference frequency.

## References

(Schutte and Seidner 1983)

## Examples

``` r
if (FALSE) { # \dontrun{
# Using different reference standards
ucnv_hz_to_semitone(220, ref_source = "UEP83")   # 12 ST above 110 Hz
ucnv_hz_to_semitone(200, ref_source = "Praat")   # 12 ST above 100 Hz
ucnv_hz_to_semitone(880, ref_source = "A4")      # 12 ST above 440 Hz

# Explicit reference frequency (overrides ref_source)
ucnv_hz_to_semitone(880, ref_freq = 440)

# Musical notes relative to A4
ucnv_hz_to_semitone(c(440, 494, 523.25, 587.33, 659.25, 698.46, 783.99, 880))

# With units
library(units)
freq <- set_units(880, Hz)
ucnv_hz_to_semitone(freq, ref_source = "A4")
} # }
```
