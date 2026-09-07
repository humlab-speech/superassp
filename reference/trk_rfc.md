# Track LP reflection coefficients

Linear Prediction analysis of audio signals using the autocorrelation
method and Durbin recursion, implemented in the *libassp* C library
(Scheffers 2012) . Returns per-frame RMS amplitudes and reflection
coefficients. Use `trk_rfc` when a lattice-filter or PARCOR
representation of the vocal tract is needed.

## Usage

``` r
trk_rfc(listOfFiles = NULL,
  beginTime = 0.0,
  centerTime = FALSE,
  endTime = 0.0,
  windowShift = 5.0,
  windowSize = 20.0,
  effectiveLength = TRUE,
  window = 'BLACKMAN',
  analysisOrder = NULL,
  preemphasis = -0.95,
  toFile = FALSE,
  explicitExt = NULL,
  outputDirectory = NULL,
  assertLossless = NULL,
  logToFile = FALSE,
  keepConverted = FALSE,
  convertOverwrites = FALSE,
  verbose = TRUE)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- beginTime:

  Numeric. Start of analysis window in seconds. Default 0 (file start).

- centerTime:

  Numeric or logical. Single-frame analysis time point in seconds;
  overrides `beginTime`, `endTime`, and `windowShift`. Default `FALSE`.

- endTime:

  Numeric. End of analysis window in seconds. Default 0 (file end).

- windowShift:

  Numeric. Frame shift in milliseconds; sets output frame rate (1000 /
  windowShift Hz). Default 5 ms.

- windowSize:

  Numeric. Analysis window size in milliseconds. Default 20 ms.

- effectiveLength:

  Logical. Make window size effective rather than exact. Default
  `FALSE`.

- window:

  Character. Analysis window function type. Default `"BLACKMAN"`. See
  [AsspWindowTypes](https://humlab-speech.github.io/superassp/reference/AsspWindowTypes.md).

- analysisOrder:

  Integer. LP order; `NULL` or 0 defaults to sample rate in kHz + 3.
  Default `NULL`.

- preemphasis:

  Numeric. Pre-emphasis factor. Default -0.95.

- toFile:

  Logical. If `TRUE`, write SSFF output files and return the count
  written (invisibly). If `FALSE`, return an `AsspDataObj`. Default
  `FALSE`.

- explicitExt:

  Character. Output file extension override. Default `NULL` (uses
  `"rfc"`).

- outputDirectory:

  Character. Directory for output files. `NULL` (default) writes
  alongside the input file.

- assertLossless:

  Character vector of additional file extensions to treat as losslessly
  encoded.

- logToFile:

  Logical. Write processing log to a file in `outputDirectory` rather
  than the console. Default `FALSE`.

- keepConverted:

  Logical. Retain intermediate transcoded files. Default `FALSE`.

- convertOverwrites:

  Logical. Allow transcoding to overwrite existing files. Default
  `FALSE`.

- verbose:

  Logical. Print per-file progress. Default `FALSE`.

## Value

If `toFile = FALSE`: an `AsspDataObj` with tracks:

- `RMS[dB]`:

  REAL32, dB, n_frames x 1. RMS amplitude of the input frame.

- `gain[dB]`:

  REAL32, dB, n_frames x 1. RMS amplitude of the LP residual.

- `RFC`:

  REAL32, dimensionless, n_frames x `analysisOrder` columns. Reflection
  (PARCOR) coefficients, one per LP order.

Frame rate: `1000 / windowShift` Hz (default 200 Hz). If
`toFile = TRUE`: integer count of files written, returned invisibly.

## Details

`analysisOrder = NULL` selects an order equal to sample rate in kHz + 3.
Pre-emphasis is applied before LP analysis and affects the residual
gain.

## References

Scheffers M (2012). “Advanced Speech Signal Processor.”
<https://sourceforge.net/projects/libassp/files/libassp/>.

## See also

wrassp::rfcana
[AsspWindowTypes](https://humlab-speech.github.io/superassp/reference/AsspWindowTypes.md)
[av::av_audio_convert](https://docs.ropensci.org/av//reference/encoding.html)

## Author

Raphael Winkelmann

Lasse Bombien

Fredrik Nylén

## Examples

``` r
path2wav <- list.files(system.file("samples", "sustained", package = "superassp"),
                       pattern = glob2rx("a1.wav"), full.names = TRUE)
res <- trk_rfc(path2wav, toFile = FALSE)
#> Applying `rfcana()` to 1 recording
matplot(seq(0, n_records(res) - 1) / sample_rate(res) +
          attr(res, "startTime"),
        res$RFC, type = "l",
        xlab = "time (s)", ylab = "Reflection coefficient values")
```
