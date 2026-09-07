# Convert Loudness Level (Phon) and Frequency to Sound Pressure Level

Converts loudness level (phon) at a given frequency to sound pressure
level (dB) according to ISO 226:2023 equal-loudness-level contours.

## Usage

``` r
ucnv_phon_and_hz_to_db(phon, freq_hz)
```

## Arguments

- phon:

  Numeric vector; loudness level in phon (20 to 90/80 phon)

- freq_hz:

  Numeric vector; frequency in Hz (20 to 12500 Hz)

## Value

Numeric vector; sound pressure level in dB (re 20 μPa)

## Details

This function implements Formula (1) from ISO 226:2023 Section 4.1.

**Loudness level (phon)** represents the perceived loudness of a sound.
This function calculates the sound pressure level required at a given
frequency to achieve a specified loudness level.

**Valid ranges**:

- Frequency: 20-12500 Hz (1/3-octave standard frequencies)

- Phon: 20-90 phon (20-4000 Hz), 20-80 phon (5000-12500 Hz)

- Below 20 phon: Informative only (near hearing threshold)

- Above limits: Limited experimental data

**Interpolation**: Frequencies between standard 1/3-octave values are
interpolated on a log-frequency scale.

**Vectorization**: Both phon and freq_hz can be vectors. If both are
vectors, they must be the same length, or one must be length 1
(recycled).

## References

( )

## See also

[`ucnv_db_and_hz_to_phon`](https://humlab-speech.github.io/superassp/reference/ucnv_db_and_hz_to_phon.md)
for the inverse conversion

## Examples

``` r
if (FALSE) { # \dontrun{
# 40 phon at 1000 Hz = 40 dB (by definition)
ucnv_phon_and_hz_to_db(40, 1000)

# Same loudness (40 phon) at different frequencies
ucnv_phon_and_hz_to_db(40, 100)   # Low frequency needs higher SPL
ucnv_phon_and_hz_to_db(40, 1000)  # Reference: 40 dB
ucnv_phon_and_hz_to_db(40, 4000)  # High frequency needs lower SPL

# Round-trip conversion check
spl_original <- 60
freq <- 1000
phon_calculated <- ucnv_db_and_hz_to_phon(spl_original, freq)
spl_recovered <- ucnv_phon_and_hz_to_db(phon_calculated, freq)
# spl_recovered ≈ spl_original

# Calculate equal-loudness contour (40 phon)
frequencies <- c(100, 200, 500, 1000, 2000, 4000, 8000)
spls <- ucnv_phon_and_hz_to_db(40, frequencies)
plot(frequencies, spls, log = "x", type = "b",
     xlab = "Frequency (Hz)", ylab = "SPL (dB)",
     main = "40 Phon Equal-Loudness Contour")
} # }
```
