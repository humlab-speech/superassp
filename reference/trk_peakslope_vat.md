# Peak slope via voiceanalysis Daless wavelet bank

Computes the peak slope acoustic parameter (Kane and Gobl 2011) using
the bit-faithful Daless wavelet bank in voiceanalysis. Faster and more
numerically stable than the pure-R
[`trk_peakslope`](https://humlab-speech.github.io/superassp/reference/trk_peakslope.md)
(which approximates Daless via db4 wavelets) for the same algorithm.

## Usage

``` r
trk_peakslope_vat(listOfFiles, beginTime = 0, endTime = 0, toFile = FALSE, explicitExt = "psv", outputDirectory = NULL, verbose = TRUE)
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

- `peak_slope`:

  REAL32, peak-slope coefficient per frame, n_frames × 1.

Frame rate: 100 Hz (10 ms hop).

## References

(Kane and Gobl 2011)

## See also

[`trk_peakslope`](https://humlab-speech.github.io/superassp/reference/trk_peakslope.md)

## Examples

``` r
if (FALSE) { # \dontrun{
trk_peakslope_vat(
  system.file("samples", "sustained", "a1.wav", package = "superassp"),
  toFile = FALSE
)
} # }
```
