# Track short-term zero-crossing rate

Computes the average of the short-term positive and negative
zero-crossing rates of audio signals using the *libassp* C library
(Scheffers 2012) . ZCR is a simple, fast measure correlated with
spectral centroid and useful for voicing detection and fricative
classification.

## Usage

``` r
trk_zcr(listOfFiles, beginTime = 0, centerTime = FALSE, endTime = 0, windowShift = 5, windowSize = 25, toFile = TRUE, explicitExt = "zcr", outputDirectory = NULL, assertLossless = NULL, logToFile = FALSE, convertOverwrites = FALSE, keepConverted = FALSE, verbose = TRUE)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- windowSize:

  Numeric. Analysis window size in milliseconds. Default 25 ms.

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

- convertOverwrites:

  Logical. Allow transcoding to overwrite existing files. Default
  `FALSE`.

- keepConverted:

  Logical. Retain intermediate transcoded files. Default `FALSE`.

- verbose:

  Logical. Show a progress bar (sequential path) or a progress-aware
  parallel apply (`pbapply`/`pbmcapply`, if installed).

## Value

If `toFile = FALSE`: an `AsspDataObj` with track:

- `ZCR[Hz]`:

  REAL32, Hz, n_frames x 1 column. Average zero-crossing rate (positive
  and negative crossings combined) per frame, expressed as a rate in Hz.

Frame rate: `1000 / windowShift` Hz (default 200 Hz). If
`toFile = TRUE`: integer count of files written, returned invisibly.

## Details

The ZCR is reported in Hz (crossings per second), averaged over the
positive and negative zero-crossing rates within each analysis window.

## See also

[wrassp::zcrana](https://rdrr.io/pkg/wrassp/man/zcrana.html)

## Author

Raphael Winkelmann

Lasse Bombien

Fredrik Nylén

## Examples

``` r
# get path to audio file
path2wav <- list.files(system.file("samples","sustained", package = "superassp"), pattern = glob2rx("a1.wav"), full.names = TRUE)

# calculate zcr values
res <- trk_zcr(path2wav, toFile=FALSE)
#> Applying `method(trk_zcr, class_character)()` to 1 recording

# plot zcr values
plot(seq(0, n_records(res) - 1) / sample_rate(res) +
      attr(res, 'startTime'),
    res[["ZCR[Hz]"]],
    type='l',
    xlab='time (s)',
    ylab='Zero Crossing Rates (Hz)')

```
