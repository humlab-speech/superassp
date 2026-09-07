# Pitch tracking via Praat's autocorrelation method

Tracks fundamental frequency (F0) using Praat's autocorrelation (AC)
pitch algorithm via pladdrr. The AC method is Praat's standard pitch
tracker and is well-suited to clean speech; prefer `trk_pitch_cc` for
noisier signals.

## Usage

``` r
trk_pitch_ac(listOfFiles, ...)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- explicitExt:

  Character. Output file extension. Default `"pac"`.

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
