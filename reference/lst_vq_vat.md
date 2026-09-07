# Per-GCI voice-quality summary via voiceanalysis

Extracts the canonical voice-source parameters (Kane and Gobl 2013) from
an audio file using voiceanalysis: NAQ (Normalised Amplitude Quotient),
QOQ (Quasi-Open Quotient), H1H2 spectral tilt, and HRF (Harmonic
Richness Factor). Internally runs SE-VQ for GCIs, IAIF for the glottal
flow derivative, then `vat_voice_quality()`.

## Usage

``` r
lst_vq_vat(listOfFiles, ...)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths.

- beginTime, endTime:

  Analysis window (seconds). Defaults: full file.

- toFile:

  Logical. If TRUE, persist results as JSTF and return paths. Default
  FALSE.

- explicitExt:

  File extension. Default "vqv".

- outputDirectory:

  Output directory. Default NULL = next to input.

- verbose:

  Logical. Default TRUE.

## Value

If `toFile = FALSE` and `length(listOfFiles) == 1`: a named list with
vectors `gci_time` (seconds), `NAQ`, `QOQ`, `H1H2`, `HRF`. For multiple
files, a list of such lists. If `toFile = TRUE`: invisible vector of
output paths.

## Details

Returns one row per GCI. To get continuous-track equivalents, interp to
a 10 ms grid via
[`stats::approx`](https://rdrr.io/r/stats/approxfun.html).

## References

(Kane and Gobl 2013)

## See also

[`lst_covarep_vq`](https://humlab-speech.github.io/superassp/reference/lst_covarep_vq.md)
