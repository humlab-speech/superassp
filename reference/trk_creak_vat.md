# Detect creaky voice using the Kane-Drugman VAT creak detector

Returns per-frame creaky phonation probability and binary label at 100
Hz using the Kane-Drugman ANN from the MATLAB Voice Analysis Toolkit
(Kane et al. 2013) . Prefer this over
[`trk_covarep_creak`](https://humlab-speech.github.io/superassp/reference/trk_covarep_creak.md)
when MATLAB-VAT parity matters.

## Usage

``` r
trk_creak_vat(listOfFiles, ...)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- threshold:

  Numeric decision threshold for binarising the posterior (default 0.3,
  matches MATLAB).

## Value

If `toFile = FALSE`: an `AsspDataObj` with tracks:

- `creak_pp`:

  REAL32, posterior probability of creaky phonation, n_frames × 1.

- `creak_bin`:

  REAL32, binary creak decision (0 / 1), n_frames × 1.

Frame rate: 100 Hz (10 ms hop). Schema mirrors
[`trk_covarep_creak()`](https://humlab-speech.github.io/superassp/reference/trk_covarep_creak.md)
for swap-in convenience.

## Details

Bit-faithful Rcpp port of the Kane-Drugman + Ishi 36-feature pipeline.
Uses the full 12-base × (static + delta + delta-delta) feature stack and
the logistic-output ANN that MATLAB's `patternnet` uses. The partial
R-side reimplementation in
[`trk_covarep_creak`](https://humlab-speech.github.io/superassp/reference/trk_covarep_creak.md)
may differ in borderline frames.

## References

(Kane et al. 2013)

## See also

[`trk_covarep_creak`](https://humlab-speech.github.io/superassp/reference/trk_covarep_creak.md)

## Examples

``` r
if (FALSE) { # \dontrun{
trk_creak_vat(
  system.file("samples", "sustained", "a1.wav", package = "superassp"),
  toFile = FALSE
)
} # }
```
