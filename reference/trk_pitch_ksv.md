# Track fundamental frequency using the KSV periodicity detector

Estimates the fundamental frequency f_(o) using the Schäfer-Vincent
periodicity detection algorithm (Schäfer-Vincent 1983) implemented in
the *libassp* C library (Scheffers 2012) . This extremum-based method is
fast and works directly on the waveform without spectral analysis.

## Usage

``` r
trk_pitch_ksv(
  listOfFiles,
  beginTime = 0,
  endTime = 0,
  windowShift = 5,
  gender = "u",
  maxF = 600,
  minF = 50,
  minAmp = 50,
  maxZCR = 3000,
  toFile = FALSE,
  explicitExt = "fo",
  outputDirectory = NULL,
  assertLossless = NULL,
  logToFile = FALSE,
  convertOverwrites = FALSE,
  keepConverted = FALSE,
  verbose = TRUE
)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- beginTime:

  Numeric. Start of analysis window in seconds. Default 0 (file start).

- endTime:

  Numeric. End of analysis window in seconds. Default 0 (file end).

- windowShift:

  Numeric. Frame shift in milliseconds; sets output frame rate (1000 /
  windowShift Hz). Default 5 ms.

- gender:

  Character. Gender-specific f_(o) search range: `"f"` (female, 80–640
  Hz), `"m"` (male, 50–400 Hz), `"u"` (unknown, 50–600 Hz, default).

- maxF:

  Numeric. Maximum f_(o) in Hz. Default 600.

- minF:

  Numeric. Minimum f_(o) in Hz. Default 50.

- minAmp:

  Numeric. Minimum waveform amplitude threshold for voiced frames.
  Default 50.

- maxZCR:

  Numeric. Maximum zero-crossing rate in Hz for voicing detection.
  Default 3000.0.

- toFile:

  Logical. If `TRUE`, write SSFF output files and return the count
  written (invisibly). If `FALSE`, return an `AsspDataObj`. Default
  `FALSE`.

- explicitExt:

  Character. Output file extension. Default `"fo"`.

- outputDirectory:

  Character. Directory for output files. `NULL` (default) writes
  alongside the input file.

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

  Logical. Print per-file progress. Default `FALSE`.

## Value

If `toFile = FALSE`: an `AsspDataObj` with track:

- `fo[Hz]`:

  REAL32, Hz, n_frames x 1 column. Estimated fundamental frequency; 0
  indicates unvoiced frames.

Frame rate: `1000 / windowShift` Hz (default 200 Hz). If
`toFile = TRUE`: integer count of files written, returned invisibly.

## Details

`gender` sets the default `minF`/`maxF` search range but is overridden
by explicit `minF`/`maxF` values. `minAmp` and `maxZCR` control voicing
detection independently of the pitch range.

## References

Schäfer-Vincent K (1983). “Pitch Period Detection and Chaining: Method
and Evaluation.” *Phonetica*, **40**(3), 177–202. ISSN 0031-8388.
[doi:10.1159/000261691](https://doi.org/10.1159/000261691) .  
  
Scheffers M (2012). “Advanced Speech Signal Processor.”
<https://sourceforge.net/projects/libassp/files/libassp/>.

## See also

[wrassp::ksvF0](https://rdrr.io/pkg/wrassp/man/ksvF0.html)

## Author

Raphael Winkelmann

Lasse Bombien

Fredrik Nylén

## Examples

``` r
# get path to audio file
path2wav <- list.files(
   system.file("samples", "sustained", package = "superassp"),
   pattern = glob2rx("a1.wav"), full.names = TRUE)

# calculate fundamental frequency contour
res <- trk_pitch_ksv(path2wav, toFile=FALSE)
#> Applying `method(trk_pitch_ksv, class_character)()` to 1 recording

# plot the fundamental frequency contour
plot(seq(0, n_records(res) - 1) / sample_rate(res) +
      attr(res, 'startTime'),
    res[["fo[Hz]"]],
    type='l',
    xlab='time (s)',
    ylab=expression(paste(f[o]," frequency (Hz)")))

```
