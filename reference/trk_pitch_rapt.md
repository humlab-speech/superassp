# Track fundamental frequency using RAPT (Robust Algorithm for Pitch Tracking)

Extracts F0 using the RAPT dynamic-programming pitch tracker from SPTK.
RAPT is a normalized cross-correlation tracker with robust
voiced/unvoiced decisions and is a reliable general-purpose choice for
speech with moderate noise. Prefer SWIPE for cleaner but noisier
signals, or PDA for higher temporal resolution.

## Usage

``` r
trk_pitch_rapt(
  listOfFiles,
  beginTime = 0,
  endTime = 0,
  windowShift = 10,
  minF = 60,
  maxF = 400,
  voicing_threshold = 0.6,
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

- windowShift:

  Numeric. Frame shift in milliseconds; sets output frame rate (1000 /
  windowShift Hz). Default 10.0 ms.

- minF:

  Numeric. Minimum F0 in Hz. Default 60.0 Hz.

- maxF:

  Numeric. Maximum F0 in Hz. Default 400.0 Hz.

- voicing_threshold:

  Numeric. Voicing decision threshold (0–1; higher = more conservative,
  fewer voiced frames). Default 0.6.

- toFile:

  Logical. If `TRUE`, write SSFF output files and return the count
  written invisibly. If `FALSE`, return an `AsspDataObj`. Default
  `TRUE`.

- explicitExt:

  Character. Output file extension. Default `"f0"`.

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
# Return F0 track in memory
f0 <- trk_pitch_rapt(wav, toFile = FALSE, verbose = FALSE)
track_names(f0)
#> [1] "f0"
head(as.data.frame(f0))
#>   frame_time f0
#> 1       0.00  0
#> 2       0.01  0
#> 3       0.02  0
#> 4       0.03  0
#> 5       0.04  0
#> 6       0.05  0

# Restrict the F0 search range
trk_pitch_rapt(wav, minF = 75, maxF = 300, toFile = FALSE, verbose = FALSE)
#> In-memory Assp Data Object
#> Format: SSFF (binary)
#> 404 records at 100 Hz
#> Duration: 4.040000 s
#> Number of tracks: 1 
#>   f0 (1 fields)
# }
```
