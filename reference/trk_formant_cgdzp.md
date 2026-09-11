# Track formants using Chirp Group Delay Zero-Phase (CGDZP) analysis

Estimates formant frequencies from the chirp group delay of a zero-phase
processed spectrum with exponential envelope warping (COVAREP
implementation). Works well on nasals and fricative transitions where
LPC-based trackers (e.g.,
[`trk_formant_burg()`](https://humlab-speech.github.io/superassp/reference/trk_formant_burg.md))
struggle because it avoids all-pole model assumptions.

## Usage

``` r
trk_formant_cgdzp(
  listOfFiles,
  beginTime = 0,
  endTime = 0,
  frameSize = 30,
  frameShift = 10,
  toFile = FALSE,
  explicitExt = "cgf",
  outputDirectory = NULL,
  verbose = TRUE
)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- frameSize:

  Numeric. Analysis window length in milliseconds. Default 30 ms.

- frameShift:

  Numeric. Frame shift in milliseconds; sets output frame rate (1000 /
  frameShift Hz). Default 10 ms.

- toFile:

  Logical. If `TRUE`, write SSFF output files and return the paths
  written invisibly. If `FALSE`, return an `AsspDataObj`. Default
  `FALSE`.

- explicitExt:

  Character. Output file extension. Default `"cgf"`.

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

If `toFile = FALSE`: an `AsspDataObj` with tracks:

- `F1`:

  SHORT (integer), first formant frequency in Hz, n_frames × 1. Zero for
  unreliable or unvoiced frames.

- `F2`:

  SHORT, second formant frequency in Hz, n_frames × 1.

- `F3`:

  SHORT, third formant frequency in Hz, n_frames × 1.

- `F4`:

  SHORT, fourth formant frequency in Hz, n_frames × 1.

- `F5`:

  SHORT, fifth formant frequency in Hz, n_frames × 1.

Frame rate: `1000 / frameShift` Hz (default 100 Hz). If `toFile = TRUE`:
character vector of output file paths, returned invisibly.

## Details

Temporal continuity is enforced with a 250 Hz maximum frame-to-frame
jump threshold. The exponential envelope warping factor is adjusted
automatically to yield exactly five formant candidates per frame. Frames
with fewer than five reliable peaks are filled from adjacent frames;
remaining zeros indicate no reliable estimate.

## Examples

``` r
if (FALSE) { # \dontrun{
# Single file, return object
formants <- trk_formant_cgdzp("speech.wav", toFile = FALSE)

# Batch process multiple files
files <- c("file1.wav", "file2.wav")
trk_formant_cgdzp(files, toFile = TRUE, outputDirectory = "output/")

# Custom frame parameters
formants <- trk_formant_cgdzp("speech.wav",
                             frameSize = 40,
                             frameShift = 5,
                             toFile = FALSE)
} # }
```
