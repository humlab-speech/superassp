# Differentiate an audio waveform

Applies a first-order finite-difference filter to audio signals using
the *libassp* C library (Scheffers 2012) . Forward, backward, and
central difference modes are supported. Useful for pre-emphasising high
frequencies or computing the derivative of a waveform before further
analysis.

## Usage

``` r
trk_afdiff(listOfFiles, ...)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- computeBackwardDifference:

  Logical. Use backward difference instead of forward. Default `FALSE`.

- computeCentralDifference:

  Logical. Use central difference instead of forward. Default `FALSE`.

- channel:

  Integer. Audio channel to process (1-based). Default `1L`.

- beginTime:

  Numeric. Start of analysis window in seconds. Default 0 (file start).

- endTime:

  Numeric. End of analysis window in seconds. Default 0 (file end).

- toFile:

  Logical. If `TRUE`, write SSFF output files and return the count
  written (invisibly). If `FALSE`, return an `AsspDataObj`. Default
  `TRUE`.

- explicitExt:

  Character. Output file extension. Default `"dif"`.

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

  Logical. Print per-file progress. Default `TRUE`.

## Value

If `toFile = FALSE`: an `AsspDataObj` with track name preserved from
libassp output (typically the same label as the input audio channel),
containing INT16 or REAL32 differentiated sample values. If
`toFile = TRUE`: integer count of files written, returned invisibly.

## Details

At most one of `computeBackwardDifference` or `computeCentralDifference`
should be `TRUE`. If both are `FALSE` (default), forward difference is
used.

## References

Scheffers M (2012). “Advanced Speech Signal Processor.”
<https://sourceforge.net/projects/libassp/files/libassp/>.

## Author

Fredrik Nylén

## Examples

``` r
path2wav <- list.files(system.file("samples", "sustained", package = "superassp"),
                       pattern = glob2rx("a1.wav"), full.names = TRUE)
res <- trk_afdiff(path2wav, toFile = FALSE)
#> Applying `method(trk_afdiff, class_character)()` to 1 recording
```
