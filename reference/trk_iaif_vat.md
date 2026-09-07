# Estimate glottal flow via IAIF using voiceanalysis

Iterative Adaptive Inverse Filtering (Alku 1992) , driven by the
bit-faithful Rcpp port in the voiceanalysis package. GCIs are detected
internally using SE-VQ.

## Usage

``` r
trk_iaif_vat(listOfFiles, beginTime = 0, endTime = 0, p = NULL, toFile = TRUE, explicitExt = "glv", outputDirectory = NULL, verbose = TRUE, ...)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- ...:

  Additional arguments (currently unused).

- p:

  LPC prediction order. `NULL` (default) sets `round(fs/1000)+2`.

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

- `glottal_flow`:

  REAL64, integrated glottal flow, n_samples × 1.

- `glottal_derivative`:

  REAL64, glottal flow derivative, n_samples × 1.

Frame rate equals the audio sample rate. Schema matches
[`trk_covarep_iaif()`](https://humlab-speech.github.io/superassp/reference/trk_covarep_iaif.md)
for swap-in convenience. If `toFile = TRUE`: invisibly returns the count
of files written.

## Details

Alternative to
[`trk_covarep_iaif`](https://humlab-speech.github.io/superassp/reference/trk_covarep_iaif.md),
which uses superassp's own native C++ IAIF implementation. The two
implementations differ in parameterisation (`trk_iaif_vat` uses the
MATLAB-VAT defaults and the SE-VQ GCI grid; `trk_covarep_iaif` uses
COVAREP defaults and a contiguous frame grid).

## References

(Alku 1992) (Kane and Gobl 2013)

## See also

[`trk_covarep_iaif`](https://humlab-speech.github.io/superassp/reference/trk_covarep_iaif.md),
[`trk_gci_vat`](https://humlab-speech.github.io/superassp/reference/trk_gci_vat.md)

## Examples

``` r
if (FALSE) { # \dontrun{
trk_iaif_vat(
  system.file("samples", "sustained", "a1.wav", package = "superassp"),
  toFile = FALSE
)
} # }
```
