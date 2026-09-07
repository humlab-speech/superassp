# Harvest Pitch Tracking (C++ implementation)

Extract F0 (fundamental frequency) using the Harvest algorithm from
WORLD vocoder (via SPTK). This is a high-performance C++ implementation
that is 2-3x faster than the Python version and requires no Python
dependencies.

Harvest is designed to be robust and accurate for speech analysis, with
good performance even on noisy signals.

All input media formats are supported via the av package, including
video files from which audio will be automatically extracted.

## Usage

``` r
trk_pitch_harvest(listOfFiles, ...)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- windowShift:

  Frame shift in milliseconds (default: 10.0)

- minF:

  Minimum F0 in Hz (default: 60.0)

- maxF:

  Maximum F0 in Hz (default: 400.0)

- voicing_threshold:

  Voicing threshold (default: 0.1, valid range: 0.02-0.2 for
  WORLD/Harvest)

- toFile:

  Write results to file (default: TRUE)

- explicitExt:

  Output file extension (default: "f0")

- parallel:

  Logical. Use parallel processing for multiple files. `NULL` (default)
  enables automatically for 2+ files.

- n_cores:

  Integer. Number of cores for parallel processing. `NULL` (default)
  uses `detectCores() - 1`.

## Value

If toFile=TRUE, returns the number of successfully processed files. If
toFile=FALSE, returns AsspDataObj or list of AsspDataObj objects.

## Examples

``` r
wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")

# \donttest{
f0 <- trk_pitch_harvest(wav, toFile = FALSE, verbose = FALSE)
head(as.data.frame(f0))
#> Warning: Package 'units' not available. Skipping unit assignment.
#>   frame_time f0
#> 1       0.00  0
#> 2       0.01  0
#> 3       0.02  0
#> 4       0.03  0
#> 5       0.04  0
#> 6       0.05  0

# Custom F0 range
trk_pitch_harvest(wav, minF = 75, maxF = 300, toFile = FALSE, verbose = FALSE)
#> In-memory Assp Data Object
#> Format: SSFF (binary)
#> 404 records at 100 Hz
#> Duration: 4.040000 s
#> Number of tracks: 1 
#>   f0 (1 fields)
# }
```
