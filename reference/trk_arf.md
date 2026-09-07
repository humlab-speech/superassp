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
res <- trk_arf(path2wav, toFile = FALSE)
#> Applying `arfana()` to 1 recording
matplot(seq(0, n_records(res) - 1) / sample_rate(res) +
          attr(res, "startTime"),
        res$ARF, type = "l",
        xlab = "time (s)", ylab = "Area function")
```
