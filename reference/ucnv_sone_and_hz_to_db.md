# Convert Sone and Hz to dB

Convenience function to convert loudness (sone) and frequency (Hz) to
sound pressure level (dB), combining ISO 532 and ISO 226:2023
conversions.

## Usage

``` r
ucnv_sone_and_hz_to_db(sone, freq_hz, method = c("zwicker", "moore-glasberg"))
```

## Arguments

- sone:

  Numeric; loudness in sones. Must be positive.

- freq_hz:

  Numeric; frequency in Hz (20-12500 Hz)

- method:

  Character; ISO 532 method for sone-to-phon conversion. One of
  `"zwicker"` (default) or `"moore-glasberg"`.

## Value

Numeric vector of sound pressure level in dB SPL (re 20 μPa)

## Details

This function combines two conversions:

1.  sone → phon using ISO 532 Zwicker or Moore-Glasberg method

2.  (phon, frequency) → SPL using ISO 226:2023 equal-loudness contours

Equivalent to:
`ucnv_phon_and_hz_to_db(ucnv_sone_to_phon(sone, method), freq_hz)`

## See also

[`ucnv_sone_to_phon()`](https://humlab-speech.github.io/superassp/reference/ucnv_sone_to_phon.md),
[`ucnv_phon_and_hz_to_db()`](https://humlab-speech.github.io/superassp/reference/ucnv_phon_and_hz_to_db.md)

## Examples

``` r
# At 1 kHz reference
ucnv_sone_and_hz_to_db(1, 1000)  # Returns 40 dB (by definition)
#> [1] 40
ucnv_sone_and_hz_to_db(2, 1000)  # Returns ~50 dB (2× loudness ≈ +10 dB)
#> [1] 50

# Same loudness at different frequencies requires different SPL
ucnv_sone_and_hz_to_db(2, 100)   # Higher dB (low freq needs more SPL)
#> [1] 71.36646
ucnv_sone_and_hz_to_db(2, 1000)  # Reference
#> [1] 50
ucnv_sone_and_hz_to_db(2, 4000)  # Lower dB (high freq more sensitive)
#> [1] 47.29502
```
