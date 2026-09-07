# Track GCI-anchored voice quality measures as a time series

Computes five glottal voice quality parameters (NAQ, QOQ, H1H2, HRF,
PSP) at each glottal closure instant and interpolates them onto a
regular 10 ms frame grid. GCIs are detected automatically via SEDREAMS
when not supplied. Use this function for time-varying voice quality
trajectories; for scalar summaries see
[`lst_covarep_vq()`](https://humlab-speech.github.io/superassp/reference/lst_covarep_vq.md).

## Usage

``` r
trk_covarep_vq_gci(listOfFiles, ...)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- gci_times:

  Optional numeric vector (or list of vectors for multiple files) of GCI
  times in seconds. If `NULL` (default), GCIs are computed internally
  via SEDREAMS.

- toFile:

  Logical. If `TRUE`, write SSFF output files and return the paths
  written invisibly. If `FALSE`, return an `AsspDataObj`. Default
  `FALSE`.

- explicitExt:

  Character. Output file extension. Default `"vqg"`.

## Value

If `toFile = FALSE`: an `AsspDataObj` with tracks:

- `naq`:

  FLOAT, Normalized Amplitude Quotient, 0–1, n_frames × 1. Low = creaky,
  high = breathy (typical range 0.5–0.9).

- `qoq`:

  FLOAT, Quasi-Open Quotient, 0–1, n_frames × 1. Fraction of pitch
  period in open phase (typical range 0.3–0.7).

- `h1h2`:

  FLOAT, H1–H2 spectral difference in dB, n_frames × 1. Positive =
  brighter/breathier, negative = darker/creakier.

- `hrf`:

  FLOAT, Harmonic Richness Factor, dimensionless, n_frames × 1. Higher
  values indicate more periodic phonation (typical range 0.5–0.95).

- `psp`:

  FLOAT, Parabolic Spectral Peak, dimensionless, n_frames × 1. Higher
  values indicate smoother spectral envelope (typical range 0.5–0.95).

Frame rate: 100 Hz (10 ms grid). If `toFile = TRUE`: character vector of
output file paths, returned invisibly.

## Details

GCI detection uses SEDREAMS on the LPC residual. Glottal flow is
obtained via IAIF. Parameter values are linearly interpolated to the 10
ms grid; frames with no voiced GCIs nearby are set to zero.

## References

There are no references for Rd macro `\insertAllCites` on this help
page.

## Examples

``` r
if (FALSE) { # \dontrun{
trk_covarep_vq_gci(
  system.file("samples", "sustained", "a1.wav", package = "superassp"),
  toFile = FALSE
)
} # }
```
