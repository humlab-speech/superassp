# Pitch tracking via Praat's subharmonic summation (SHS) method

Tracks fundamental frequency (F0) using Praat's subharmonic summation
(SHS) algorithm via pladdrr. SHS works in the frequency domain by
summing subharmonically compressed spectral representations and is
robust to spectral-domain interference; prefer it over AC/CC for speech
with strong vocal fry or very low F0.

## Usage

``` r
trk_pitch_shs(listOfFiles, beginTime = 0, endTime = 0, time_step = 0.01, minimum_f0 = 50, maximum_f0 = 500, maximum_frequency_components = 1250, maximum_number_of_subharmonics = 15, number_of_candidates = 15, compression_factor = 0.84, number_of_points_per_octave = 48, windowShape = "Gaussian1", relativeWidth = 1, toFile = TRUE, explicitExt = "psh", outputDirectory = NULL, verbose = TRUE)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- beginTime:

  Numeric. Start of analysis window in seconds. Default 0 (file start).

- endTime:

  Numeric. End of analysis window in seconds. Default 0 (file end).

- time_step:

  Numeric. Frame shift in seconds; sets output frame rate (1 / time_step
  Hz). Default 0.01 s (100 Hz).

- minimum_f0:

  Numeric. Lower F0 bound in Hz. Default 50 Hz.

- maximum_f0:

  Numeric. Upper F0 bound (ceiling) in Hz. Default 500 Hz.

- maximum_frequency_components:

  Numeric. Highest frequency included in the spectral analysis in Hz.
  Default 1250 Hz.

- maximum_number_of_subharmonics:

  Integer. Maximum number of subharmonic levels to sum. Default 15.

- number_of_candidates:

  Integer. Maximum pitch candidates per frame. Default 15.

- compression_factor:

  Numeric. Spectral compression factor applied during subharmonic
  summation. Default 0.84.

- number_of_points_per_octave:

  Integer. Spectral resolution in the log-frequency representation
  (points per octave). Default 48.

- windowShape:

  Character. Window shape applied to audio before loading. Default
  `"Gaussian1"`.

- relativeWidth:

  Numeric. Relative width of the window. Default 1.0.

- toFile:

  Logical. If `TRUE`, write SSFF output files and return the count
  written (invisibly). If `FALSE`, return an `AsspDataObj`. Default
  `TRUE`.

- explicitExt:

  Character. Output file extension. Default `"psh"`.

- outputDirectory:

  Character. Directory for output files. `NULL` (default) writes
  alongside the input file.

- verbose:

  Logical. Print per-file progress. Default `TRUE`.

## Value

If `toFile = FALSE`: an `AsspDataObj` with track:

- `F0`:

  REAL32, Hz, n_frames x 1. Fundamental frequency. 0 encodes unvoiced
  frames.

Frame rate: `1 / time_step` Hz (default 100 Hz). If `toFile = TRUE`:
integer count of files written, returned invisibly.

## See also

[`trk_pitch_cc`](https://humlab-speech.github.io/superassp/reference/trk_pitch_cc.md),
[`trk_pitch_ac`](https://humlab-speech.github.io/superassp/reference/trk_pitch_ac.md),
[`trk_pitch_spinet`](https://humlab-speech.github.io/superassp/reference/trk_pitch_spinet.md)

## Examples

``` r
if (FALSE) { # \dontrun{
trk_pitch_shs(
  system.file("samples", "sustained", "a1.wav", package = "superassp"),
  toFile = FALSE
)
} # }
```
