# Track LP-derived vocal tract area function coefficients

Linear Prediction analysis of audio signals using the autocorrelation
method and Durbin recursion, implemented in the *libassp* C library
(Scheffers 2012) . Returns per-frame RMS amplitudes and vocal tract area
function (ARF) coefficients derived from the reflection coefficients.
Useful for vocal tract modelling applications.

## Usage

``` r
trk_arf(listOfFiles = NULL,
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

  Start time for the extracted portion in seconds. Default: NULL
  (beginning of signal). Note: uses `beginTime`/`endTime` (seconds)
  matching DSP function conventions, unlike
  [`read_audio()`](https://humlab-speech.github.io/superassp/reference/read_audio.md)
  which uses `begin`/`end`.

- centerTime:

  Numeric or logical. Single-frame analysis time point in seconds;
  overrides `beginTime`, `endTime`, and `windowShift`. Default `FALSE`.

- endTime:

  The end time of the section of the sound files that should be analysed
  (in seconds). Use 0 for end of file.

- windowShift:

  Numeric. Frame shift in milliseconds; sets output frame rate
  (`1000 / windowShift` Hz). Default 5.0 ms (200 Hz). Must be strictly
  less than 32 ms (the 512-sample analysis window at 16 kHz). Values
  other than the training default (5 ms) may slightly reduce accuracy.

- windowSize:

  Numeric. Smoothing filter window size in milliseconds, applied to both
  median (periodicity) and mean (F0) post-processing filters. Default 15
  ms.

- effectiveLength:

  Logical. Make window size effective rather than exact. Default
  `FALSE`.

- window:

  Character. Analysis window function type. Default `"BLACKMAN"`. See
  [AsspWindowTypes](https://humlab-speech.github.io/superassp/reference/AsspWindowTypes.md)
  for supported types.

- analysisOrder:

  Integer. Number of lag coefficients per frame. `0` sets order to
  sample rate in kHz + 3 (e.g. 19 for 16 kHz audio). Default 0.

- preemphasis:

  Numeric. Pre-emphasis factor (-1 \<= val \<= 0); default is
  sample-rate- and nominalF1-dependent.

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

  Logical. Show a progress bar (sequential path) or a progress-aware
  parallel apply (`pbapply`/`pbmcapply`, if installed).

## Value

If `toFile = FALSE`: an `AsspDataObj` with tracks:

- `RMS[dB]`:

  REAL32, dB, n_frames x 1. RMS amplitude of the input frame.

- `gain[dB]`:

  REAL32, dB, n_frames x 1. RMS amplitude of the LP residual.

- `ARF`:

  REAL32, dimensionless area ratios, n_frames x `analysisOrder` columns.
  Area function coefficients derived from reflection coefficients.

Frame rate: `1000 / windowShift` Hz (default 200 Hz). If
`toFile = TRUE`: integer count of files written, returned invisibly.

## Details

ARF coefficients are derived from the LP reflection coefficients via the
standard area ratio transformation. See `trk_rfc` for parameter details.

## References

Scheffers M (2012). “Advanced Speech Signal Processor.”
<https://sourceforge.net/projects/libassp/files/libassp/>.

## See also

[wrassp::rfcana](https://rdrr.io/pkg/wrassp/man/rfcana.html)
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
res <- trk_arf(path2wav, toFile = FALSE)
#> Applying `arfana()` to 1 recording
matplot(seq(0, n_records(res) - 1) / sample_rate(res) +
          attr(res, "startTime"),
        res$ARF, type = "l",
        xlab = "time (s)", ylab = "Area function")
```
