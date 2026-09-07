# Convert dB and Hz Directly to Sone

Convenience function to convert sound pressure level (dB) and frequency
(Hz) directly to loudness (sone), combining ISO 226:2023 and ISO 532
conversions.

## Usage

``` r
ucnv_db_and_hz_to_sone(
  spl_db,
  freq_hz,
  method = c("zwicker", "moore-glasberg")
)
```

## Arguments

- spl_db:

  Numeric; sound pressure level in dB SPL (re 20 μPa)

- freq_hz:

  Numeric; frequency in Hz (20-12500 Hz)

- method:

  Character; ISO 532 method for phon-to-sone conversion. One of
  `"zwicker"` (default) or `"moore-glasberg"`.

## Value

Numeric vector of loudness values in sones

## Details

This function combines two conversions:

1.  (SPL, frequency) → phon using ISO 226:2023 equal-loudness contours

2.  phon → sone using ISO 532 Zwicker or Moore-Glasberg method

Equivalent to:
`ucnv_phon_to_sone(ucnv_db_and_hz_to_phon(spl_db, freq_hz), method)`

## See also

[`ucnv_db_and_hz_to_phon()`](https://humlab-speech.github.io/superassp/reference/ucnv_db_and_hz_to_phon.md),
[`ucnv_phon_to_sone()`](https://humlab-speech.github.io/superassp/reference/ucnv_phon_to_sone.md)

## Examples

``` r
# At 1 kHz reference
ucnv_db_and_hz_to_sone(40, 1000)  # Returns 1.0 (by definition)
#> [1] 1
ucnv_db_and_hz_to_sone(50, 1000)  # Returns ~2.0 (+10 dB ≈ 2× loudness)
#> [1] 2

# Low frequency requires more SPL for same loudness
ucnv_db_and_hz_to_sone(60, 100)   # Lower sone value
#> [1] 0.610957
ucnv_db_and_hz_to_sone(60, 1000)  # Higher sone value
#> [1] 4
```
