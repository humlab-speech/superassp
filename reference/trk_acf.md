# Track short-term autocorrelation function

Computes the short-term autocorrelation function (ACF) of audio signals
using the *libassp* C library (Scheffers 2012) . Useful as a front-end
feature for voicing detection and LP-based analysis. Prefer `trk_acf`
over manual lag computation when frame-synchronous output in SSFF format
is needed.

## Usage

``` r
trk_acf(listOfFiles, beginTime = 0, centerTime = FALSE, endTime = 0, windowShift = 5, windowSize = 20, effectiveLength = TRUE, window = "BLACKMAN", analysisOrder = 0, energyNormalization = FALSE, lengthNormalization = FALSE, toFile = TRUE, explicitExt = "acf", outputDirectory = NULL, assertLossless = NULL, logToFile = FALSE, keepConverted = FALSE, convertOverwrites = FALSE, verbose = TRUE)
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

  Logical. Make window size effective rather than exact. Default `TRUE`.

- window:

  Character. Analysis window function type. Default `"BLACKMAN"`. See
  [AsspWindowTypes](https://humlab-speech.github.io/superassp/reference/AsspWindowTypes.md)
  for supported types.

- analysisOrder:

  Integer. Number of lag coefficients per frame. `0` sets order to
  sample rate in kHz + 3 (e.g. 19 for 16 kHz audio). Default 0.

- energyNormalization:

  Logical. Compute energy-normalised ACF. Default `FALSE`.

- lengthNormalization:

  Logical. Compute length-normalised ACF. Default `FALSE`.

- toFile:

  Logical. If `TRUE`, write SSFF output files and return the count
  written (invisibly). If `FALSE`, return an `AsspDataObj`. Default
  `TRUE`.

- explicitExt:

  Character. Output file extension. Default `"acf"`.

- outputDirectory:

  Character. Directory for output files. `NULL` (default) writes
  alongside the input file.

- verbose:

  Logical. Print per-file progress. Default `TRUE`.

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

## Value

If `toFile = FALSE`: an `AsspDataObj` with track:

- `ACF`:

  REAL32, dimensionless, n_frames x `analysisOrder` columns.
  Autocorrelation coefficients at lags 0 … analysisOrder-1.

Frame rate: `1000 / windowShift` Hz (default 200 Hz). If
`toFile = TRUE`: integer count of files written, returned invisibly.

## Details

`analysisOrder = 0` selects an order equal to the sample rate in kHz +
3. Energy normalisation divides each frame's ACF by its lag-0 value.
Length normalisation divides by frame length. Both can be combined.

## References

Scheffers M (2012). “Advanced Speech Signal Processor.”
<https://sourceforge.net/projects/libassp/files/libassp/>.

## See also

[wrassp::acfana](https://rdrr.io/pkg/wrassp/man/acfana.html)

[AsspWindowTypes](https://humlab-speech.github.io/superassp/reference/AsspWindowTypes.md)

[av::av_audio_convert](https://docs.ropensci.org/av//reference/encoding.html)

## Examples

``` r
# get path to audio file
path2wav <- list.files(system.file("samples","sustained", package = "superassp"), pattern = glob2rx("a1.wav"), full.names = TRUE)

# calculate short-term autocorrelation
res <- trk_acf(path2wav, toFile=FALSE)
#> Applying `method(trk_acf, class_character)()` to 1 recording

# plot short-term autocorrelation values
matplot(seq(0, n_records(res) - 1) / sample_rate(res) +
        attr(res, 'startTime'),
        res$ACF,
        type='l',
        xlab='time (s)',
        ylab='Short-term autocorrelation values')

```
