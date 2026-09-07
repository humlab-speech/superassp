# Detect creaky voice (vocal fry) per frame

Estimates the probability of creaky phonation in each analysis frame
using a shallow neural network trained on LPC residual and spectral
features from the COVAREP toolkit. Returns both a continuous posterior
and a binary label. Prefer this over HNR-based measures when
distinguishing creak from breathiness.

## Usage

``` r
trk_covarep_creak(listOfFiles, beginTime = 0, endTime = 0, toFile = FALSE, explicitExt = "crk", outputDirectory = NULL, verbose = TRUE)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- toFile:

  Logical. If `TRUE`, write SSFF output files and return the paths
  written invisibly. If `FALSE`, return an `AsspDataObj`. Default
  `FALSE`.

- explicitExt:

  Character. Output file extension. Default `"crk"`.

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

- `creak_pp`:

  FLOAT, posterior probability of creaky phonation, 0–1, n_frames × 1.
  Values above 0.5 indicate creak.

- `creak_bin`:

  FLOAT, binary creakiness label (0 = non-creaky, 1 = creaky), n_frames
  × 1. Threshold applied at 0.3.

Frame rate: 100 Hz (10 ms hop). If `toFile = TRUE`: character vector of
output file paths, returned invisibly.

## Details

Features are 12 static predictors (H2H1, LPC residual peak, ZCR, frame
energy, power SD, spectral band energies) plus first and second
derivatives (36-D total), normalized and passed through a two-layer ANN
with tanh/sigmoid activations. A 3-frame median filter is applied to the
posterior before thresholding.

## References

There are no references for Rd macro `\insertAllCites` on this help
page.

## Examples

``` r
if (FALSE) { # \dontrun{
trk_covarep_creak(
  system.file("samples", "sustained", "a1.wav", package = "superassp"),
  toFile = FALSE
)
} # }
```
