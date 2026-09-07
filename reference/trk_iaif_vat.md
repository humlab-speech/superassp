# Estimate glottal flow via IAIF using voiceanalysis

Iterative Adaptive Inverse Filtering (Alku 1992) , driven by the
bit-faithful Rcpp port in the voiceanalysis package. GCIs are detected
internally using SE-VQ.

## Usage

``` r
trk_iaif_vat(listOfFiles, ...)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- ...:

  Additional arguments (currently unused).

- p:

  LPC prediction order. `NULL` (default) sets `round(fs/1000)+2`.

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
