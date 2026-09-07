# Track LP-derived log area ratios

Linear Prediction analysis of audio signals using the autocorrelation
method and Durbin recursion, implemented in the *libassp* C library
(Scheffers 2012) . Returns per-frame RMS amplitudes and log area ratio
(LAR) coefficients. LARs are a log-domain reparameterisation of
reflection coefficients that can be more numerically stable near the
unit circle.

## Usage

``` r
trk_lar(listOfFiles = NULL,
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

- `LAR`:

  REAL32, dimensionless, n_frames x `analysisOrder` columns. Log area
  ratios derived from reflection coefficients.

Frame rate: `1000 / windowShift` Hz (default 200 Hz). If
`toFile = TRUE`: integer count of files written, returned invisibly.

## Details

LAR coefficients are defined as log((1 + k) / (1 - k)) where k is the
reflection coefficient. See `trk_rfc` for parameter details.

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
res <- trk_lar(path2wav, toFile = FALSE)
#> Applying `larana()` to 1 recording
matplot(seq(0, n_records(res) - 1) / sample_rate(res) +
          attr(res, "startTime"),
        res$LAR, type = "l",
        xlab = "time (s)", ylab = "Log area ratios")
```
