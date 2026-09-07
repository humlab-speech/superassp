# Estimate spectral envelope using the True Envelope (Teager energy) method

Extracts a smooth spectral envelope and compact cepstral representation
per frame using cepstral liftering. The True Envelope method separates
the slow-varying spectral shape from fine harmonic structure, making it
useful for voice quality and vocal tract characterization independent of
F0.

## Usage

``` r
trk_covarep_env_te(listOfFiles, ...)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- frameSize:

  Numeric. Analysis window length in milliseconds. Default 30 ms.

- frameShift:

  Numeric. Frame shift in milliseconds; sets output frame rate (1000 /
  frameShift Hz). Default 5 ms.

- cep_order:

  Integer. Cepstral order for envelope reconstruction. Controls
  smoothness of the recovered envelope: lower values yield smoother
  envelopes (safe range 12–40). Default 24.

- toFile:

  Logical. If `TRUE`, write SSFF output files and return the paths
  written invisibly. If `FALSE`, return an `AsspDataObj`. Default
  `FALSE`.

- explicitExt:

  Character. Output file extension. Default `"ete"`.

## Value

If `toFile = FALSE`: an `AsspDataObj` with tracks:

- `env_te`:

  FLOAT, log10 spectral envelope, n_frames × 513 (DFT size 1024,
  half-spectrum). Units: log10(amplitude).

- `env_cc`:

  FLOAT, cepstral coefficients c0…c\_{cep_order}, n_frames ×
  (cep_order + 1). Dimensionless.

Frame rate: `1000 / frameShift` Hz (default 200 Hz). If `toFile = TRUE`:
character vector of output file paths, returned invisibly.

## Details

Envelope is estimated by cepstral liftering: log FFT magnitude is
inverse-DFT'd, the top `cep_order + 1` cepstral coefficients are kept,
and the result is forward-DFT'd to recover a smooth log spectrum.

## References

There are no references for Rd macro `\insertAllCites` on this help
page.

## Examples

``` r
if (FALSE) { # \dontrun{
trk_covarep_env_te(
  system.file("samples", "sustained", "a1.wav", package = "superassp"),
  toFile = FALSE
)
} # }
```
