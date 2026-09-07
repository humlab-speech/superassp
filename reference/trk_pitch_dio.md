# DIO Pitch Tracking (C++ implementation)

Extract F0 using the DIO algorithm from the WORLD vocoder (via SPTK).
DIO is designed for high-quality pitch extraction for speech synthesis
applications.

## Usage

``` r
trk_pitch_dio(listOfFiles, ...)
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

## Value

If toFile=TRUE, returns the number of successfully processed files. If
toFile=FALSE, returns AsspDataObj or list of AsspDataObj objects.

## Examples

``` r
wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")

# \donttest{
f0 <- trk_pitch_dio(wav, toFile = FALSE, verbose = FALSE)
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
trk_pitch_dio(wav, minF = 80, maxF = 350, toFile = FALSE, verbose = FALSE)
#> In-memory Assp Data Object
#> Format: SSFF (binary)
#> 404 records at 100 Hz
#> Duration: 4.040000 s
#> Number of tracks: 1 
#>   f0 (1 fields)
# }
```
