# Track fundamental frequency using SWIPE (Sawtooth Waveform Inspired Pitch Estimator)

Extracts F0 by matching the speech spectrum against sawtooth waveform
templates across candidate F0 values via SPTK. SWIPE is particularly
effective for noisy or challenging recording conditions where
cross-correlation trackers (RAPT, Snack) are less reliable. Its default
voicing threshold (0.3) is lower than RAPT's (0.6), yielding more voiced
frames.

## Usage

``` r
trk_pitch_swipe(
  listOfFiles,
  beginTime = 0,
  endTime = 0,
  windowShift = 10,
  minF = 60,
  maxF = 400,
  voicing_threshold = 0.3,
  toFile = TRUE,
  explicitExt = "f0",
  outputDirectory = NULL,
  verbose = TRUE,
  parallel = NULL,
  n_cores = NULL
)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- voicing_threshold:

  Numeric. Voicing decision threshold (0–1). Default 0.3 (more
  permissive than RAPT). Increase toward 0.5 to reduce false voiced
  frames.

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

If `toFile = FALSE`: an `AsspDataObj` with track:

- `f0`:

  REAL32, fundamental frequency in Hz, n_frames × 1. Zero indicates
  unvoiced frames.

Frame rate: `1000 / windowShift` Hz (default 100 Hz). If
`toFile = TRUE`: integer count of files written, returned invisibly.

## Examples

``` r
wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")

# \donttest{
f0 <- trk_pitch_swipe(wav, toFile = FALSE, verbose = FALSE)
head(as.data.frame(f0))
#>   frame_time f0
#> 1       0.00  0
#> 2       0.01  0
#> 3       0.02  0
#> 4       0.03  0
#> 5       0.04  0
#> 6       0.05  0

# Custom range and voicing threshold
trk_pitch_swipe(wav, minF = 100, maxF = 500, voicing_threshold = 0.4,
                toFile = FALSE, verbose = FALSE)
#> In-memory Assp Data Object
#> Format: SSFF (binary)
#> 404 records at 100 Hz
#> Duration: 4.040000 s
#> Number of tracks: 1 
#>   f0 (1 fields)
# }
```
