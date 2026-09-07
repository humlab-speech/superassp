# Track fundamental frequency using REAPER (Robust Epoch And Pitch EstimatoR)

Extracts F0 and glottal closure instant (epoch) times simultaneously
using REAPER from SPTK. REAPER's EpochTracker uses two-pass
correlation + dynamic programming for joint epoch and pitch estimation.
When epoch times are also needed, this is more efficient than running a
separate pitchmark detector. For epoch-only output, see
[`trk_pitchmark_reaper`](https://humlab-speech.github.io/superassp/reference/trk_pitchmark_reaper.md).

## Usage

``` r
trk_pitch_reaper(listOfFiles, beginTime = 0, endTime = 0, windowShift = 10, minF = 60, maxF = 400, voicing_threshold = 0.9, toFile = TRUE, explicitExt = "f0", outputDirectory = NULL, verbose = TRUE, parallel = NULL, n_cores = NULL)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- voicing_threshold:

  Numeric. Voicing decision threshold (0–1; higher = more conservative).
  Default 0.9.

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

Additionally, the following attributes are set on the returned object:
`epochs` (numeric vector of GCI times in seconds), `n_epochs` (integer
count), `polarity` (signal polarity estimate). Frame rate:
`1000 / windowShift` Hz (default 100 Hz). If `toFile = TRUE`: integer
count of files written, returned invisibly.

## Examples

``` r
wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")

# \donttest{
# F0 plus glottal-epoch information
result <- trk_pitch_reaper(wav, toFile = FALSE, verbose = FALSE)
head(as.data.frame(result))
#>   frame_time f0
#> 1       0.00  0
#> 2       0.01  0
#> 3       0.02  0
#> 4       0.03  0
#> 5       0.04  0
#> 6       0.05  0
epochs <- attr(result, "epochs")
# }
```
