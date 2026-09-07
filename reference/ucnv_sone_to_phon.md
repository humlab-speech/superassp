# Convert Sone to Phon

Convert loudness (sone) to loudness level (phon) using ISO 532 methods.

## Usage

``` r
ucnv_sone_to_phon(sone, method = c("zwicker", "moore-glasberg"))
```

## Arguments

- sone:

  Numeric; loudness in sones. Must be positive.

- method:

  Character; conversion method. One of:

  - `"zwicker"` (default): ISO 532-1 Zwicker method with piecewise
    formula

  - `"moore-glasberg"`: ISO 532-2 Moore-Glasberg method with lookup
    table

## Value

Numeric vector of loudness level values in phons

## Details

The Zwicker method (ISO 532-1) uses:

- For sone \< 1: `phon = 40 x sone^0.35`

- For sone \>= 1: `phon = 40 + 10 x log2(sone)`

The Moore-Glasberg method (ISO 532-2) uses a lookup table with 23
reference points and log-linear interpolation.

By definition, 1 sone = 40 phons (a 1 kHz tone at 40 dB SPL).

### Vectorization

This function is fully vectorized and can process vectors of sone values
efficiently.

## See also

[`ucnv_phon_to_sone()`](https://humlab-speech.github.io/superassp/reference/ucnv_phon_to_sone.md),
[`ucnv_db_and_hz_to_phon()`](https://humlab-speech.github.io/superassp/reference/ucnv_db_and_hz_to_phon.md),
[`ucnv_phon_and_hz_to_db()`](https://humlab-speech.github.io/superassp/reference/ucnv_phon_and_hz_to_db.md)

## Examples

``` r
# Reference value
ucnv_sone_to_phon(1)  # Returns 40 (by definition)
#> [1] 40

# Doubling property: 2x loudness ~= +10 phon
ucnv_sone_to_phon(2)  # Returns ~50
#> [1] 50
ucnv_sone_to_phon(4)  # Returns ~60
#> [1] 60

# Low levels (below 1 sone)
ucnv_sone_to_phon(0.15)  # Returns ~20
#> [1] 20.59169

# Vectorized
ucnv_sone_to_phon(c(0.15, 1, 4, 16))
#> [1] 20.59169 40.00000 60.00000 80.00000

# Moore-Glasberg method
ucnv_sone_to_phon(1, method = "moore-glasberg")
#> [1] 40
```
