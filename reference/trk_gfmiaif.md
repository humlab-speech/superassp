# Decompose speech into vocal tract, glottis, and lip radiation LP filters (GFM-IAIF)

Returns per-frame LP coefficient tracks for the vocal tract, glottis,
and lip radiation filters using GFM-IAIF source-filter separation
(Perrotin and d'Alessandro 2019) . Prefer this over
[`trk_covarep_iaif()`](https://humlab-speech.github.io/superassp/reference/trk_covarep_iaif.md)
when you need frame-level LP coefficients rather than a sample-domain
glottal waveform.

## Usage

``` r
trk_gfmiaif(listOfFiles, ...)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- centerTime:

  Logical. If `TRUE`, timestamps refer to window centres; if `FALSE`, to
  window starts. Default `FALSE`.

- windowShift:

  Numeric. Frame shift in milliseconds; sets output frame rate (1000 /
  windowShift Hz). Default 10.0 ms.

- windowSize:

  Numeric. Analysis window length in milliseconds. Must be large enough
  for `nv + 1` samples at the audio sample rate. Default 32.0 ms.

- nv:

  Integer. Vocal tract LPC order (1–100). Higher values capture more
  spectral detail but increase compute. Default 48.

- ng:

  Integer. Glottis LPC order. The value 3 is recommended by the original
  authors; departing from it may degrade estimates. Default 3.

- d:

  Numeric. Leaky integration coefficient for lip radiation (0.9–0.999).
  Default 0.99.

- window:

  Character. Analysis window type: `"HANN"`, `"HAMMING"`, or
  `"BLACKMAN"`. Default `"HANN"`.

- explicitExt:

  Character. Output file extension. Default `"gfm"`.

- toFile:

  Logical. If `TRUE`, write SSFF output files and return the count
  written. If `FALSE`, return an `AsspDataObj` (single file only).
  Default `TRUE`.

## Value

If `toFile = FALSE`: an `AsspDataObj` with tracks:

- `av_0` … `av_{nv}`:

  REAL64, vocal tract LP coefficients, n_frames × (nv + 1).
  Dimensionless.

- `ag_0` … `ag_{ng}`:

  REAL64, glottis LP coefficients, n_frames × (ng + 1). Dimensionless.

- `al_0`, `al_1`:

  REAL64, lip radiation LP coefficients, n_frames × 2. Dimensionless.

Frame rate: `1000 / windowShift` Hz (default 100 Hz). If
`toFile = TRUE`: integer count of files written.

## Details

GFM-IAIF extends classical IAIF with a wide-band glottis model that
captures both glottal formant and spectral tilt. It estimates three LP
filters per frame (vocal tract, glottis, lip radiation) via alternating
inverse filtering steps.

`ng = 3` is strongly recommended. When `toFile = FALSE`, only
single-file input is permitted. The `ng` warning is issued automatically
for other values.

## References

Perrotin O, d'Alessandro C (2019). “Glottal flow model estimation by
modified Iterative Adaptive Inverse Filtering.” In *Proceedings of
Interspeech 2019*, 234–238.
[doi:10.21437/Interspeech.2019-1499](https://doi.org/10.21437/Interspeech.2019-1499)
.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage
result <- trk_gfmiaif("speech.wav", toFile = FALSE)

# Lower vocal tract order for faster processing
result <- trk_gfmiaif("speech.wav", nv = 24, windowShift = 5.0, toFile = FALSE)

# Batch processing to SSFF files
files <- c("file1.wav", "file2.wav", "file3.mp3")
trk_gfmiaif(files, toFile = TRUE)
} # }
```
