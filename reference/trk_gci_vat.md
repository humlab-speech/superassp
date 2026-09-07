# Detect glottal closure instants (GCIs) using SE-VQ via voiceanalysis

Detects GCIs using the improved SEDREAMS / SE-VQ algorithm of Kane &
Gobl (2013) (Kane and Gobl 2013) as implemented by the voiceanalysis
package (bit-faithful Rcpp port of the MATLAB Voice Analysis Toolkit).
Supports both fixed-F0 and the variable-F0 variant optimised for highly
expressive speech.

## Usage

``` r
trk_gci_vat(listOfFiles, ...)
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
