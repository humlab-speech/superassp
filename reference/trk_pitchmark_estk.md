# Detect glottal closure instants in laryngograph signals using ESTk pitchmark

Finds pitchmark times in laryngograph (EGG) or speech waveforms using
the Edinburgh Speech Tools pitchmark algorithm (zero-crossing detection
on the filtered, differentiated signal). Prefer
`protoscribe::draft_pitchmark()` for new code; this function is retained
for backwards compatibility only.

## Usage

``` r
trk_pitchmark_estk(
  listOfFiles,
  beginTime = 0,
  endTime = 0,
  lx_low_frequency = 400,
  lx_low_order = 19,
  lx_high_frequency = 40,
  lx_high_order = 19,
  df_low_frequency = 1000,
  df_low_order = 19,
  median_order = 19,
  fill = FALSE,
  min_period = 0.003,
  max_period = 0.02,
  def_period = 0.01,
  invert = FALSE,
  to_f0 = FALSE,
  toFile = TRUE,
  explicitExt = NULL,
  outputDirectory = NULL,
  verbose = TRUE,
  parallel = NULL,
  n_cores = NULL,
  use_cpp = TRUE
)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- lx_low_frequency:

  Numeric. Low-pass cutoff in Hz for initial denoising filter. Default
  400 Hz.

- lx_low_order:

  Integer. Order of the initial low-pass FIR filter. Default 19.

- lx_high_frequency:

  Numeric. High-pass cutoff in Hz to remove low-frequency swell. Default
  40 Hz.

- lx_high_order:

  Integer. Order of the high-pass FIR filter. Default 19.

- df_low_frequency:

  Numeric. Low-pass cutoff in Hz applied to the differentiated signal.
  Default 1000 Hz.

- df_low_order:

  Integer. Order of the differentiated-signal low-pass filter. Set to 0
  to disable. Default 19.

- median_order:

  Integer. Order of the median smoother on the differentiated signal.
  Set to 0 to disable. Default 19.

- fill:

  Logical. If `TRUE`, post-process pitchmarks: remove marks closer than
  `min_period`, interpolate gaps larger than `max_period`. Default
  `FALSE`.

- min_period:

  Numeric. Minimum pitch period in seconds (used when `fill = TRUE`).
  Default 0.003 s (\\\approx\\333 Hz max F0).

- max_period:

  Numeric. Maximum pitch period in seconds (used when `fill = TRUE`).
  Default 0.02 s (\\\approx\\50 Hz min F0).

- def_period:

  Numeric. Default pitch period for interpolated marks (used when
  `fill = TRUE`). Default 0.01 s (100 Hz).

- invert:

  Logical. Invert signal polarity before processing (use for upside-down
  EGG recordings). Default `FALSE`.

- to_f0:

  Logical. If `TRUE`, return F0 values derived from pitchmark intervals
  instead of raw pitchmark times. Default `FALSE`.

- toFile:

  Logical. If `TRUE`, write output files and return the count written
  invisibly. If `FALSE`, return results as R objects. Default `TRUE`.

- explicitExt:

  Character. Output file extension. Default `"pm"` (or `"f0"` when
  `to_f0 = TRUE`).

- parallel:

  Logical. Use parallel processing for multiple files. `NULL` (default)
  enables automatically for 2+ files.

- n_cores:

  Integer. Number of cores for parallel processing. `NULL` (default)
  uses `detectCores() - 1`.

- use_cpp:

  Logical. Use C++ implementation (default `TRUE`). Setting `FALSE`
  falls back to the ESTK binary (slower, requires temp files).

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

If `toFile = TRUE`: integer count of files written, returned invisibly.
If `toFile = FALSE` and `to_f0 = FALSE`: a data frame (single file) or
list of data frames with pitchmark times in seconds. If `toFile = FALSE`
and `to_f0 = TRUE`: a data frame (or list) with F0 values derived from
inter-pitchmark intervals.

## Note

**DEPRECATED**. Use `protoscribe::draft_pitchmark()` instead, which
follows the `draft_` naming convention and integrates with reindeer
workflows.

Pitchmarking accuracy is highest with lossless formats (WAV, FLAC).
Lossy formats (MP3, AAC) may degrade detection.

## References

(Black et al. 2020)

(Macon and Taylor 1997)

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic pitchmarking of laryngograph file
trk_pitchmark_estk("recording.egg")

# Pitchmark with filling (interpolate unvoiced regions)
trk_pitchmark_estk("recording.wav", fill = TRUE,
               min = 0.003, max = 0.02, def = 0.01)

# Extract F0 from pitchmarks
f0_data <- trk_pitchmark_estk("recording.egg", to_f0 = TRUE, toFile = FALSE)

# Process video file (extracts audio automatically)
trk_pitchmark_estk("interview.mp4")

# Batch processing with parallel execution
files <- c("rec1.wav", "rec2.wav", "rec3.wav")
trk_pitchmark_estk(files, parallel = TRUE, n_cores = 4)

# Inverted laryngograph signal
trk_pitchmark_estk("recording.egg", invert = TRUE)
} # }
```
