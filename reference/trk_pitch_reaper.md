# Track fundamental frequency using REAPER (Robust Epoch And Pitch EstimatoR)

Extracts F0 and glottal closure instant (epoch) times simultaneously
using REAPER from SPTK. REAPER's EpochTracker uses two-pass
correlation + dynamic programming for joint epoch and pitch estimation.
When epoch times are also needed, this is more efficient than running a
separate pitchmark detector. For epoch-only output, see
[`trk_pitchmark_reaper`](https://humlab-speech.github.io/superassp/reference/trk_pitchmark_reaper.md).

## Usage

``` r
trk_pitch_reaper(listOfFiles, ...)
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
#> Warning: Package 'units' not available. Skipping unit assignment.
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
