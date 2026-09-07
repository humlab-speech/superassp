# Track fundamental frequency using the ESTk PDA algorithm

Extracts F0 using the Edinburgh Speech Tools super-resolution Pitch
Detection Algorithm (PDA), a normalized cross-correlation tracker with
DP-based peak tracking optimized for speech with variable pitch. Offers
sub-frame resolution and explicit voiced/unvoiced transition thresholds.

## Usage

``` r
trk_pitch_pda(listOfFiles, beginTime = 0, endTime = 0, windowShift = 5, windowSize = 10, minF = 40, maxF = 400, decimation = 4, noise_floor = 120, min_v2uv_coef_thresh = 0.75, v2uv_coef_thresh_ratio = 0.85, uv2v_coef_thresh = 0.88, anti_doubling_thresh = 0.77, peak_tracking = FALSE, toFile = FALSE, explicitExt = "pda", outputDirectory = NULL, verbose = TRUE)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- windowShift:

  Numeric. Frame shift in milliseconds; sets output frame rate (1000 /
  windowShift Hz). Default 5.0 ms.

- windowSize:

  Numeric. Analysis window length in milliseconds. Default 10.0 ms.

- minF:

  Numeric. Minimum F0 in Hz. Default 40.0 Hz.

- maxF:

  Numeric. Maximum F0 in Hz. Default 400.0 Hz.

- decimation:

  Integer. Correlation decimation factor (higher = faster but coarser).
  Default 4.

- noise_floor:

  Numeric. Silence threshold; frames below this RMS are treated as
  unvoiced. Default 120.

- min_v2uv_coef_thresh:

  Numeric. Minimum correlation for voiced-to-unvoiced transition (0–1).
  Default 0.75.

- v2uv_coef_thresh_ratio:

  Numeric. Correlation ratio threshold for voiced-to-unvoiced transition
  (0–1). Default 0.85.

- uv2v_coef_thresh:

  Numeric. Correlation threshold for unvoiced-to-voiced transition
  (0–1). Default 0.88.

- anti_doubling_thresh:

  Numeric. Threshold to suppress F0 octave doublings (0–1). Default
  0.77.

- peak_tracking:

  Logical. If `TRUE`, enable peak tracking for smoother F0 contours.
  Default `FALSE`.

- toFile:

  Logical. If `TRUE`, write SSFF output files and return the paths
  written invisibly. If `FALSE`, return an `AsspDataObj`. Default
  `FALSE`.

- explicitExt:

  Character. Output file extension. Default `"pda"`.

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

- `F0`:

  REAL32, fundamental frequency in Hz, n_frames × 1. Zero indicates
  unvoiced frames.

Frame rate: `1000 / windowShift` Hz (default 200 Hz). If
`toFile = TRUE`: output file path(s), returned invisibly.

## References

(Medan et al. 1991)

(Bagshaw et al. 1993)

## Examples

``` r
if (FALSE) { # \dontrun{
# Extract F0 from audio file
f0_data <- trk_pitch_pda("recording.wav", toFile = FALSE)

# Process with custom parameters
trk_pitch_pda("speech.mp3", minF = 75, maxF = 300, windowShift = 10)

# Enable peak tracking for smoother contours
trk_pitch_pda("recording.wav", peak_tracking = TRUE, toFile = FALSE)

# Process video file (extracts audio)
trk_pitch_pda("interview.mp4", toFile = FALSE)
} # }
```
