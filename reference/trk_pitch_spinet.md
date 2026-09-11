# Pitch tracking via Praat's SPINET method

Tracks fundamental frequency (F0) using Praat's SPINET (SPectral
INtegration and Evaluation of Temporal patterns) algorithm via pladdrr.
SPINET uses a bank of gammatone filters followed by temporal
integration, mimicking auditory processing; it can detect F0 in
conditions where time-domain methods fail.

## Usage

``` r
trk_pitch_spinet(
  listOfFiles,
  beginTime = 0,
  endTime = 0,
  time_step = 0.005,
  window_length = 0.04,
  minimum_filter_frequency = 70,
  maximum_filter_frequency = 5000,
  number_of_filters = 250,
  maximum_f0 = 500,
  number_of_candidates = 15,
  windowShape = "Gaussian1",
  relativeWidth = 1,
  toFile = TRUE,
  explicitExt = "psp",
  outputDirectory = NULL,
  verbose = TRUE
)
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
  Hz). Default 0.005 s (200 Hz).

- window_length:

  Numeric. Duration of the integration window in seconds. Default 0.04
  s.

- minimum_filter_frequency:

  Numeric. Center frequency of the lowest gammatone filter in Hz.
  Default 70 Hz.

- maximum_filter_frequency:

  Numeric. Center frequency of the highest gammatone filter in Hz.
  Default 5000 Hz.

- number_of_filters:

  Integer. Number of gammatone filters in the bank. Default 250.

- maximum_f0:

  Numeric. Upper F0 bound (ceiling) in Hz. Default 500 Hz.

- number_of_candidates:

  Integer. Maximum pitch candidates per frame. Default 15.

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

  Character. Output file extension. Default `"psp"`.

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

Frame rate: `1 / time_step` Hz (default 200 Hz). If `toFile = TRUE`:
integer count of files written, returned invisibly.

## Details

SPINET is slower than AC/CC/SHS due to the gammatone filterbank but can
outperform them on atypical voices (very low F0, strong noise). No lower
F0 bound is exposed — the minimum detectable F0 is determined by
`minimum_filter_frequency`.

## See also

[`trk_pitch_cc`](https://humlab-speech.github.io/superassp/reference/trk_pitch_cc.md),
[`trk_pitch_ac`](https://humlab-speech.github.io/superassp/reference/trk_pitch_ac.md),
[`trk_pitch_shs`](https://humlab-speech.github.io/superassp/reference/trk_pitch_shs.md)

## Examples

``` r
if (FALSE) { # \dontrun{
trk_pitch_spinet(
  system.file("samples", "sustained", "a1.wav", package = "superassp"),
  toFile = FALSE
)
} # }
```
