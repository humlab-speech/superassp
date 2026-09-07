# Sound intensity contour

Computes a time-series intensity (dB SPL) contour using Praat's
intensity algorithm via pladdrr. Window length is derived from
`minimal_f0_frequency` to ensure at least one pitch period per frame.

## Usage

``` r
trk_intensity(listOfFiles, ...)
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
