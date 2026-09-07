# Track fundamental frequency using SWIPE (Sawtooth Waveform Inspired Pitch Estimator)

Extracts F0 by matching the speech spectrum against sawtooth waveform
templates across candidate F0 values via SPTK. SWIPE is particularly
effective for noisy or challenging recording conditions where
cross-correlation trackers (RAPT, Snack) are less reliable. Its default
voicing threshold (0.3) is lower than RAPT's (0.6), yielding more voiced
frames.

## Usage

``` r
trk_pitch_swipe(listOfFiles, ...)
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
#> Warning: Package 'units' not available. Skipping unit assignment.
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
