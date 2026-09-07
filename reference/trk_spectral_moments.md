# Spectral moments (CoG, SD, skewness, kurtosis)

Extracts the four spectral moments — center of gravity, standard
deviation, skewness, and kurtosis — as time-series tracks via Praat's
spectrogram in pladdrr. Spectral moments characterize spectral shape and
are widely used for consonant place-of-articulation and voice quality
analysis.

## Usage

``` r
trk_spectral_moments(listOfFiles, beginTime = 0, endTime = 0, windowLength = 0.005, maximum_frequency = 0, time_step = 0.005, frequency_step = 20, power = 2, windowShape = "Gaussian1", relativeWidth = 1, toFile = TRUE, explicitExt = "spm", outputDirectory = NULL, verbose = TRUE)
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

- beginTime:

  Start time for the extracted portion in seconds. Default: NULL
  (beginning of signal). Note: uses `beginTime`/`endTime` (seconds)
  matching DSP function conventions, unlike
  [`read_audio()`](https://humlab-speech.github.io/superassp/reference/read_audio.md)
  which uses `begin`/`end`.

- endTime:

  The end time of the section of the sound files that should be analysed
  (in seconds). Use 0 for end of file.

- outputDirectory:

  The directory where the slice file should be stored. If not defiled
  (NULL), the sparse slice file will placed in the same folder as the
  media file.

- verbose:

  Logical. Show a progress bar (sequential path) or a progress-aware
  parallel apply (`pbapply`/`pbmcapply`, if installed).

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
