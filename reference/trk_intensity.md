# Sound intensity contour

Computes a time-series intensity (dB SPL) contour using Praat's
intensity algorithm via pladdrr. Window length is derived from
`minimal_f0_frequency` to ensure at least one pitch period per frame.

## Usage

``` r
trk_intensity(
  listOfFiles,
  beginTime = 0,
  endTime = 0,
  time_step = 0,
  minimal_f0_frequency = 50,
  subtract_mean = TRUE,
  windowShape = "Gaussian1",
  relativeWidth = 1,
  toFile = TRUE,
  explicitExt = "int",
  outputDirectory = NULL,
  verbose = TRUE
)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- time_step:

  Numeric. Frame shift in seconds; sets output frame rate (1 / time_step
  Hz). Set to 0 for Praat's automatic choice. Default 0.

- minimal_f0_frequency:

  Numeric. Minimum expected pitch in Hz; determines the effective
  analysis window length (3 / minimal_f0_frequency). Default 50 Hz.

- subtract_mean:

  Logical. Subtract mean intensity to correct DC offset. Default `TRUE`.

- windowShape:

  Character. Window shape applied to the extracted audio segment.
  Default `"Gaussian1"`.

- relativeWidth:

  Numeric. Relative width of the window. Default 1.0.

- toFile:

  Logical. If `TRUE`, write SSFF output files and return the count
  written (invisibly). If `FALSE`, return an `AsspDataObj`. Default
  `TRUE`.

- explicitExt:

  Character. Output file extension. Default `"int"`.

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

- `intensity`:

  REAL32, dB SPL, n_frames x 1. Sound pressure level contour. 0 encodes
  undefined frames.

Frame rate: `1 / time_step` Hz (Praat automatic when `time_step = 0`).
If `toFile = TRUE`: integer count of files written, returned invisibly.

## Details

Uses Praat's intensity-from-sound algorithm. Window length is set
automatically to 3 / `minimal_f0_frequency`. Increase
`minimal_f0_frequency` for a shorter, more time-resolved window;
decrease for better low-frequency coverage.

## Examples

``` r
if (FALSE) { # \dontrun{
# Analyze intensity for a single file
intensity <- trk_intensity("speech.wav", toFile = FALSE)

# Batch process multiple files
files <- c("speech1.wav", "speech2.wav")
trk_intensity(files, toFile = TRUE, outputDirectory = "output/")

# With time windowing
intensity <- trk_intensity("speech.wav",
                            beginTime = 1.0,
                            endTime = 3.0,
                            toFile = FALSE)
} # }
```
