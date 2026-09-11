# Detect glottal closure instants (GCIs) using SE-VQ via voiceanalysis

Detects GCIs using the improved SEDREAMS / SE-VQ algorithm of Kane &
Gobl (2013) (Kane and Gobl 2013) as implemented by the voiceanalysis
package (bit-faithful Rcpp port of the MATLAB Voice Analysis Toolkit).
Supports both fixed-F0 and the variable-F0 variant optimised for highly
expressive speech.

## Usage

``` r
trk_gci_vat(
  listOfFiles,
  beginTime = 0,
  endTime = 0,
  var_f0 = FALSE,
  f0_min = 20,
  f0_max = 500,
  use_creak = FALSE,
  toFile = TRUE,
  explicitExt = "gciv",
  outputDirectory = NULL,
  verbose = TRUE
)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- var_f0:

  Logical. If `TRUE`, use the variable-F0 SE-VQ variant (recommended for
  expressive speech). Default `FALSE`.

- f0_min:

  Minimum F0 in Hz (default 20).

- f0_max:

  Maximum F0 in Hz (default 500).

- use_creak:

  Logical. If `TRUE`, run voiceanalysis's creak detector and feed its
  decisions into the SE-VQ creaky post-processing step. Default `FALSE`.

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

If `toFile = FALSE`: an `AsspDataObj` with tracks:

- `gci_sample`:

  INT32, GCI sample indices (1-based) in the resampled signal, n_gci ×
  1.

- `residual`:

  REAL32, LP residual signal (normalised), n_samples × 1.

The `gci_sample` track is a sparse event list; the `residual` track is
per-sample at the audio rate. If `toFile = TRUE`: invisibly returns the
count of files written.

## References

(Kane and Gobl 2013)

## See also

[`trk_covarep_vq_gci`](https://humlab-speech.github.io/superassp/reference/trk_covarep_vq_gci.md)

## Examples

``` r
if (FALSE) { # \dontrun{
trk_gci_vat(
  system.file("samples", "sustained", "a1.wav", package = "superassp"),
  toFile = FALSE
)
} # }
```
