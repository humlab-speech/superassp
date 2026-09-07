# Pitch tracking via Praat's cross-correlation method

Tracks fundamental frequency (F0) using Praat's cross-correlation (CC)
pitch algorithm via pladdrr. The CC method correlates short-term
waveforms across time and is more robust to noise than autocorrelation
at the cost of slightly higher compute time.

## Usage

``` r
trk_pitch_cc(listOfFiles, ...)
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

- minimum_f0:

  Numeric. Lower F0 bound in Hz. Default 75 Hz.

- maximum_f0:

  Numeric. Upper F0 bound (ceiling) in Hz. Default 600 Hz.

- very_accurate:

  Logical. Use slower, higher-accuracy candidate search. Default `TRUE`.

- number_of_candidates:

  Integer. Maximum pitch candidates per frame. Default 15.

- silence_threshold:

  Numeric. Frames with amplitude below this fraction of the global
  maximum are treated as silent (0–1). Default 0.03.

- voicing_threshold:

  Numeric. Minimum strength for a frame to be voiced (0–1). Default
  0.45.

- octave_cost:

  Numeric. Penalty per octave above `minimum_f0` to discourage
  high-frequency candidates. Default 0.01.

- octave_jump_cost:

  Numeric. Penalty for octave jumps between adjacent frames. Default
  0.35.

- voiced_voiceless_cost:

  Numeric. Penalty for voiced/unvoiced transitions. Default 0.14.

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

  Character. Output file extension. Default `"pcc"`.

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

## See also

[`trk_pitch_ac`](https://humlab-speech.github.io/superassp/reference/trk_pitch_ac.md),
[`trk_pitch_shs`](https://humlab-speech.github.io/superassp/reference/trk_pitch_shs.md),
[`trk_pitch_spinet`](https://humlab-speech.github.io/superassp/reference/trk_pitch_spinet.md)

## Examples

``` r
if (FALSE) { # \dontrun{
trk_pitch_cc(
  system.file("samples", "sustained", "a1.wav", package = "superassp"),
  toFile = FALSE
)
} # }
```
