# Detect glottal closure instants using REAPER (pitch marks)

Extracts glottal closure instants (GCIs) as a binary indicator track
using the REAPER EpochTracker (Talkin 2019) from SPTK. GCIs are mapped
from irregular epoch times onto a regular grid at `windowShift`
intervals. Use
[`trk_pitch_reaper`](https://humlab-speech.github.io/superassp/reference/trk_pitch_reaper.md)
instead when F0 is also needed (avoids re-running REAPER).

## Usage

``` r
trk_pitchmark_reaper(
  listOfFiles,
  beginTime = 0,
  endTime = 0,
  windowShift = 10,
  minF = 40,
  maxF = 500,
  voicing_threshold = 0.9,
  toFile = TRUE,
  explicitExt = "rpm",
  outputDirectory = NULL,
  verbose = TRUE
)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- windowShift:

  Numeric. Frame shift in milliseconds for the output indicator grid.
  Default 10.0 ms.

- minF:

  Numeric. Minimum F0 in Hz for the internal pitch estimator. Lower
  values allow lower-pitched voices but may increase false positives.
  Default 40.0 Hz.

- maxF:

  Numeric. Maximum F0 in Hz for the internal pitch estimator. Default
  500.0 Hz.

- voicing_threshold:

  Numeric. Voicing decision threshold (0–1; higher = more conservative).
  Default 0.9.

- toFile:

  Logical. If `TRUE`, write SSFF output files and return the count
  written invisibly. If `FALSE`, return an `AsspDataObj`. Default
  `TRUE`.

- explicitExt:

  Character. Output file extension. Default `"rpm"`.

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

- `pm`:

  INT16, binary pitch mark indicator (0 = no GCI, 1 = GCI), n_frames
  × 1. Frame spacing = windowShift ms.

Additional attributes on the returned object: `epoch_times` (raw GCI
times in seconds), `n_epochs`, `polarity`. Frame rate:
`1000 / windowShift` Hz (default 100 Hz). If `toFile = TRUE`: integer
count of files written, returned invisibly.

## Details

Raw epoch times at irregular intervals are available via
`attr(result, "epoch_times")` when `toFile = FALSE`.

## References

Talkin D (2019). “REAPER: Robust epoch and pitch estimator.”
<https://github.com/google/REAPER>.

## See also

[`trk_pitch_reaper`](https://humlab-speech.github.io/superassp/reference/trk_pitch_reaper.md)
for F0 extraction (also extracts epochs as attributes)
[`trk_pitchmark_estk`](https://humlab-speech.github.io/superassp/reference/trk_pitchmark_estk.md)
for ESTK-based pitch mark detection

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage - extract pitch marks
trk_pitchmark_reaper("speech.wav")

# Get pitch marks without writing to file
result <- trk_pitchmark_reaper("speech.wav", toFile = FALSE)
pm_track <- result$pm  # Binary indicator (0 or 1)

# Access raw epoch times (irregular intervals)
epoch_times <- attr(result, "epoch_times")  # Times in seconds
n_epochs <- attr(result, "n_epochs")        # Number of epochs

# Adjust F0 range for low-pitched voice
trk_pitchmark_reaper("bass_voice.wav", minF = 50, maxF = 300)

# Process specific time window
trk_pitchmark_reaper("long_recording.wav",
              beginTime = 1.0,
              endTime = 5.0)

# Batch processing
files <- list.files(pattern = "\\.wav$", full.names = TRUE)
n_success <- trk_pitchmark_reaper(files,
                            outputDirectory = "results/",
                            verbose = TRUE)
message("Processed ", n_success, " files")
} # }
```
