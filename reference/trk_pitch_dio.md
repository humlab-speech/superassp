# DIO Pitch Tracking (C++ implementation)

Extract F0 using the DIO algorithm from the WORLD vocoder (via SPTK).
DIO is designed for high-quality pitch extraction for speech synthesis
applications.

## Usage

``` r
trk_pitch_dio(listOfFiles, beginTime = 0, endTime = 0, windowShift = 10, minF = 60, maxF = 400, voicing_threshold = 0.1, toFile = TRUE, explicitExt = "f0", outputDirectory = NULL, verbose = TRUE, parallel = NULL, n_cores = NULL)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- voicing_threshold:

  Voicing threshold (default: 0.1, valid range: 0.02-0.2 for WORLD/DIO)

- parallel:

  Logical. Use parallel processing for multiple files. `NULL` (default)
  enables automatically for 2+ files.

- n_cores:

  Integer. Number of cores for parallel processing. `NULL` (default)
  uses `detectCores() - 1`.

- beginTime:

  Start time for the extracted portion in seconds. Default: NULL
  (beginning of signal). Note: uses `beginTime`/`endTime` (seconds)
  matching DSP function conventions, unlike
  [`read_audio()`](https://humlab-speech.github.io/superassp/reference/read_audio.md)
  which uses `begin`/`end`.

- endTime:

  The end time of the section of the sound files that should be analysed
  (in seconds). Use 0 for end of file.

- windowShift:

  Numeric. Frame shift in milliseconds; sets output frame rate
  (`1000 / windowShift` Hz). Default 5.0 ms (200 Hz). Must be strictly
  less than 32 ms (the 512-sample analysis window at 16 kHz). Values
  other than the training default (5 ms) may slightly reduce accuracy.

- minF:

  Numeric. Minimum F0 in Hz for the internal pitch estimator. Lower
  values allow lower-pitched voices but may increase false positives.
  Default 40.0 Hz.

- maxF:

  Numeric. Maximum F0 in Hz to treat as voiced. Default 400 Hz (speech).
  Must be \<= 2093.75 Hz (model maximum; C7). For music, use 2093.75.

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

If toFile=TRUE, returns the number of successfully processed files. If
toFile=FALSE, returns AsspDataObj or list of AsspDataObj objects.

## Examples

``` r
wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")

# \donttest{
f0 <- trk_pitch_dio(wav, toFile = FALSE, verbose = FALSE)
head(as.data.frame(f0))
#>   frame_time f0
#> 1       0.00  0
#> 2       0.01  0
#> 3       0.02  0
#> 4       0.03  0
#> 5       0.04  0
#> 6       0.05  0

# Custom F0 range
trk_pitch_dio(wav, minF = 80, maxF = 350, toFile = FALSE, verbose = FALSE)
#> In-memory Assp Data Object
#> Format: SSFF (binary)
#> 404 records at 100 Hz
#> Duration: 4.040000 s
#> Number of tracks: 1 
#>   f0 (1 fields)
# }
```
