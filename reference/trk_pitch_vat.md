# Track fundamental frequency using SRH via the voiceanalysis package

Returns F0, voiced/unvoiced decisions, and SRH amplitude per frame using
the SRH algorithm (Drugman and Alwan 2011) . Prefer this over
[`trk_pitch_srh`](https://humlab-speech.github.io/superassp/reference/trk_pitch_srh.md)
when MATLAB-VAT parity matters.

## Usage

``` r
trk_pitch_vat(listOfFiles, ...)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- minF:

  Minimum F0 in Hz (default 50).

- maxF:

  Maximum F0 in Hz (default 500).

## Value

If `toFile = FALSE`: an `AsspDataObj` with tracks:

- `f0`:

  REAL32, fundamental frequency in Hz, n_frames × 1. Zero indicates
  unvoiced frames.

- `vad`:

  REAL32, voiced/unvoiced decision (0 = unvoiced, 1 = voiced), n_frames
  × 1.

- `srh_val`:

  REAL32, SRH amplitude per frame, n_frames × 1.

Frame rate: 100 Hz (fixed 10 ms hop). Audio is resampled to 16 kHz
internally to match the MATLAB pipeline. If `toFile = TRUE`: invisibly
returns the count of files written.

## Details

Bit-faithful Rcpp port of the Summation of Residual Harmonics pitch
tracker from the original Kane MATLAB Voice Analysis Toolkit. The two
implementations differ in framing and smoothing details; use
`trk_pitch_srh` for the native superassp version.

## References

(Drugman and Alwan 2011) (Kane and Gobl 2013)

## See also

[`trk_pitch_srh`](https://humlab-speech.github.io/superassp/reference/trk_pitch_srh.md)

## Examples

``` r
if (FALSE) { # \dontrun{
trk_pitch_vat(
  system.file("samples", "sustained", "a1.wav", package = "superassp"),
  toFile = FALSE
)
} # }
```
