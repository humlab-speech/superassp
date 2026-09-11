# Track fundamental frequency using SRH via the voiceanalysis package

Returns F0, voiced/unvoiced decisions, and SRH amplitude per frame using
the SRH algorithm (Drugman and Alwan 2011) . Prefer this over
[`trk_pitch_srh`](https://humlab-speech.github.io/superassp/reference/trk_pitch_srh.md)
when MATLAB-VAT parity matters.

## Usage

``` r
trk_pitch_vat(
  listOfFiles,
  beginTime = 0,
  endTime = 0,
  minF = 50,
  maxF = 500,
  toFile = TRUE,
  explicitExt = "f0v",
  outputDirectory = NULL,
  verbose = TRUE
)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- minF:

  Minimum F0 in Hz (default 50).

- maxF:

  Maximum F0 in Hz (default 500).

- beginTime:

  Start time for the extracted portion in seconds. Default: NULL
  (beginning of signal). Note: uses `beginTime`/`endTime` (seconds)
  matching DSP function conventions, unlike
  [`read_audio()`](https://humlab-speech.github.io/superassp/reference/read_audio.md)
  which uses `begin`/`end`.

- endTime:

  The end time of the section of the sound files that should be analysed
  (in seconds). Use 0 for end of file.

- toFile:

  Logical. If `TRUE`, write SSFF output files and return the count
  written. If `FALSE`, return an `AsspDataObj` (single file only).
  Default `TRUE`.

- explicitExt:

  By default, a character "d" will be prepended to the file name suffix
  when writing the output to file. The user can also specify an explicit
  extension which will be used instead.

- outputDirectory:

  The directory where the slice file should be stored. If not defiled
  (NULL), the sparse slice file will placed in the same folder as the
  media file.

- verbose:

  Logical. Show a progress bar (sequential path) or a progress-aware
  parallel apply (`pbapply`/`pbmcapply`, if installed).

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
