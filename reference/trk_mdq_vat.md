# Track Maxima Dispersion Quotient (MDQ) for breathy/tense voice discrimination

Returns MDQ (Kane and Gobl 2013) per frame at 100 Hz — higher values
indicate more dispersed wavelet maxima (breathier voice).

## Usage

``` r
trk_mdq_vat(listOfFiles, beginTime = 0, endTime = 0, toFile = FALSE, explicitExt = "mdq", outputDirectory = NULL, verbose = TRUE)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- beginTime:

  Start time for the extracted portion in seconds. Default: NULL
  (beginning of signal). Note: uses `beginTime`/`endTime` (seconds)
  matching DSP function conventions, unlike
  [`read_audio()`](https://humlab-speech.github.io/superassp/reference/read_audio.md)
  which uses `begin`/`end`.

- endTime:

  The end time of the section of the sound files that should be analysed
  (in seconds). Use 0 for end of file.

- toFile:

  Logical. If `TRUE`, write SSFF output files and return the count
  written. If `FALSE`, return an `AsspDataObj` (single file only).
  Default `TRUE`.

- explicitExt:

  By default, a character "d" will be prepended to the file name suffix
  when writing the output to file. The user can also specify an explicit
  extension which will be used instead.

- outputDirectory:

  The directory where the slice file should be stored. If not defiled
  (NULL), the sparse slice file will placed in the same folder as the
  media file.

- verbose:

  Logical. Show a progress bar (sequential path) or a progress-aware
  parallel apply (`pbapply`/`pbmcapply`, if installed).

## Value

If `toFile = FALSE`: an `AsspDataObj` with one track:

- `mdq`:

  REAL32, MDQ value, n_frames × 1. Higher values indicate more dispersed
  wavelet maxima (breathier voice).

Frame rate: 100 Hz (10 ms hop).

## Details

MDQ is computed per-GCI using the bit-faithful Daless wavelet bank from
voiceanalysis, then resampled to a fixed 100 Hz frame grid (10 ms hop)
so it can sit alongside other `trk_*` tracks.

## References

(Kane and Gobl 2013)

## Examples

``` r
if (FALSE) { # \dontrun{
trk_mdq_vat(
  system.file("samples", "sustained", "a1.wav", package = "superassp"),
  toFile = FALSE
)
} # }
```
