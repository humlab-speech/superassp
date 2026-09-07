# Detect creaky voice using the Kane-Drugman VAT creak detector

Returns per-frame creaky phonation probability and binary label at 100
Hz using the Kane-Drugman ANN from the MATLAB Voice Analysis Toolkit
(Kane et al. 2013) . Prefer this over
[`trk_covarep_creak`](https://humlab-speech.github.io/superassp/reference/trk_covarep_creak.md)
when MATLAB-VAT parity matters.

## Usage

``` r
trk_creak_vat(listOfFiles, beginTime = 0, endTime = 0, threshold = 0.3, toFile = FALSE, explicitExt = "crv", outputDirectory = NULL, verbose = TRUE)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- threshold:

  Numeric decision threshold for binarising the posterior (default 0.3,
  matches MATLAB).

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
