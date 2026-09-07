# Spectral moments (CoG, SD, skewness, kurtosis)

Extracts the four spectral moments — center of gravity, standard
deviation, skewness, and kurtosis — as time-series tracks via Praat's
spectrogram in pladdrr. Spectral moments characterize spectral shape and
are widely used for consonant place-of-articulation and voice quality
analysis.

## Usage

``` r
trk_spectral_moments(listOfFiles, ...)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- windowLength:

  Numeric. Spectrogram analysis window length in seconds. Default 0.005
  s.

- maximum_frequency:

  Numeric. Highest frequency included in moment calculations in Hz. Set
  to 0 to use the Nyquist frequency. Default 0.

- time_step:

  Numeric. Frame shift in seconds; sets output frame rate (1 / time_step
  Hz). Default 0.005 s (200 Hz).

- frequency_step:

  Numeric. Frequency resolution of the spectrogram in Hz. Default 20 Hz.

- power:

  Numeric. Exponent applied to the amplitude spectrum before computing
  moments. Default 2 (power spectrum).

- windowShape:

  Character. Window shape for audio extraction. Default `"Gaussian1"`.

- relativeWidth:

  Numeric. Relative width of the extraction window. Default 1.0.

- toFile:

  Logical. If `TRUE`, write SSFF output files and return the paths
  written (invisibly). If `FALSE`, return an `AsspDataObj`. Default
  `TRUE`.

- explicitExt:

  Character. Output file extension. Default `"spm"`.

## Value

If `toFile = FALSE`: an `AsspDataObj` with tracks:

- `cog`:

  REAL32, Hz, n_frames x 1. Spectral center of gravity (spectral mean).

- `sd`:

  REAL32, Hz, n_frames x 1. Spectral standard deviation (spread).

- `skewness`:

  REAL32, dimensionless, n_frames x 1. Spectral skewness (asymmetry).

- `kurtosis`:

  REAL32, dimensionless, n_frames x 1. Spectral kurtosis (peakedness).

Frame rate: `1 / time_step` Hz (default 200 Hz). If `toFile = TRUE`:
character vector of output file paths, returned invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
# Extract spectral moments
result <- trk_spectral_moments("speech.wav", toFile = FALSE)

# Access tracks
plot(result$cog, type = "l", main = "Center of Gravity")

# Write to file
trk_spectral_moments("speech.wav", toFile = TRUE)
} # }
```
