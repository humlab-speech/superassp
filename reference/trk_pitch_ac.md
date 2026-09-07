# Pitch tracking via Praat's autocorrelation method

Tracks fundamental frequency (F0) using Praat's autocorrelation (AC)
pitch algorithm via pladdrr. The AC method is Praat's standard pitch
tracker and is well-suited to clean speech; prefer `trk_pitch_cc` for
noisier signals.

## Usage

``` r
trk_pitch_ac(listOfFiles, beginTime = 0, endTime = 0, time_step = 0.005, minimum_f0 = 75, maximum_f0 = 600, very_accurate = TRUE, number_of_candidates = 15, silence_threshold = 0.03, voicing_threshold = 0.45, octave_cost = 0.01, octave_jump_cost = 0.35, voiced_voiceless_cost = 0.14, windowShape = "Gaussian1", relativeWidth = 1, toFile = TRUE, explicitExt = "pac", outputDirectory = NULL, verbose = TRUE)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- explicitExt:

  Character. Output file extension. Default `"pac"`.

- beginTime:

  Start time for the extracted portion in seconds. Default: NULL
  (beginning of signal). Note: uses `beginTime`/`endTime` (seconds)
  matching DSP function conventions, unlike
  [`read_audio()`](https://humlab-speech.github.io/superassp/reference/read_audio.md)
  which uses `begin`/`end`.

- endTime:

  The end time of the section of the sound files that should be analysed
  (in seconds). Use 0 for end of file.

- time_step:

  Numeric. Frame shift in seconds; sets output frame rate (1 / time_step
  Hz). Set to 0 for Praat's automatic choice. Default 0.

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

  Numeric. Voicing decision threshold (0–1). Default 0.3 (more
  permissive than RAPT). Increase toward 0.5 to reduce false voiced
  frames.

- octave_cost:

  Numeric. Penalty per octave above `minimum_f0` to discourage
  high-frequency candidates. Default 0.01.

- octave_jump_cost:

  Numeric. Penalty for octave jumps between adjacent frames. Default
  0.35.

- voiced_voiceless_cost:

  Numeric. Penalty for voiced/unvoiced transitions. Default 0.14.

- windowShape:

  Character. Window shape applied to the extracted audio segment.
  Default `"Gaussian1"`.

- relativeWidth:

  Numeric. Relative width of the extraction window. Default 1.0.

- toFile:

  Logical. If `TRUE`, write SSFF output files and return the count
  written. If `FALSE`, return an `AsspDataObj` (single file only).
  Default `TRUE`.

- outputDirectory:

  The directory where the slice file should be stored. If not defiled
  (NULL), the sparse slice file will placed in the same folder as the
  media file.

- verbose:

  Logical. Show a progress bar (sequential path) or a progress-aware
  parallel apply (`pbapply`/`pbmcapply`, if installed).

## Value

If `toFile = FALSE`: an `AsspDataObj` with track:

- `F0`:

  REAL32, Hz, n_frames x 1. Fundamental frequency. 0 encodes unvoiced
  frames.

Frame rate: `1 / time_step` Hz (default 200 Hz). If `toFile = TRUE`:
integer count of files written, returned invisibly.

## See also

[`trk_pitch_cc`](https://humlab-speech.github.io/superassp/reference/trk_pitch_cc.md),
[`trk_pitch_shs`](https://humlab-speech.github.io/superassp/reference/trk_pitch_shs.md),
[`trk_pitch_spinet`](https://humlab-speech.github.io/superassp/reference/trk_pitch_spinet.md)

## Examples

``` r
if (FALSE) { # \dontrun{
trk_pitch_ac(
  system.file("samples", "sustained", "a1.wav", package = "superassp"),
  toFile = FALSE
)
} # }
```
