# Convert Phon to Sone

Convert loudness level (phon) to loudness (sone) using ISO 532 methods.

## Usage

``` r
ucnv_phon_to_sone(phon, method = c("zwicker", "moore-glasberg"))
```

## Arguments

- phon:

  Numeric; loudness level in phons. Must be non-negative.

- method:

  Character; conversion method. One of:

  - `"zwicker"` (default): ISO 532-1 Zwicker method with piecewise
    formula

  - `"moore-glasberg"`: ISO 532-2 Moore-Glasberg method with lookup
    table

## Value

Numeric vector of loudness values in sones

## Details

The Zwicker method (ISO 532-1) uses:

- For phon \< 40: `sone = (phon / 40)^(1/0.35)`

- For phon ≥ 40: `sone = 2^((phon - 40) / 10)`

The Moore-Glasberg method (ISO 532-2) uses a lookup table with 23
reference points and log-linear interpolation.

By definition, 1 sone = 40 phons (a 1 kHz tone at 40 dB SPL).

### Vectorization

This function is fully vectorized and can process vectors of phon values
efficiently.

## See also

[`ucnv_sone_to_phon()`](https://humlab-speech.github.io/superassp/reference/ucnv_sone_to_phon.md),
[`ucnv_db_and_hz_to_phon()`](https://humlab-speech.github.io/superassp/reference/ucnv_db_and_hz_to_phon.md),
[`ucnv_phon_and_hz_to_db()`](https://humlab-speech.github.io/superassp/reference/ucnv_phon_and_hz_to_db.md)

## Examples

``` r
# Reference value
ucnv_phon_to_sone(40)  # Returns 1.0 (by definition)
#> [1] 1

# Doubling property: +10 phon ≈ 2× loudness
ucnv_phon_to_sone(50)  # Returns ~2.0
#> [1] 2
ucnv_phon_to_sone(60)  # Returns ~4.0
#> [1] 4

# Low levels (below 40 phon)
ucnv_phon_to_sone(20)  # Returns ~0.15
#> [1] 0.1380112

# Vectorized
ucnv_phon_to_sone(c(20, 40, 60, 80))
#> [1]  0.1380112  1.0000000  4.0000000 16.0000000

# Moore-Glasberg method
ucnv_phon_to_sone(40, method = "moore-glasberg")
#> [1] 1
```
