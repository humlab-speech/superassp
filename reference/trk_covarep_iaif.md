# Extract glottal flow waveform using Iterative Adaptive Inverse Filtering (IAIF)

Separates the glottal source from the vocal tract contribution using the
IAIF algorithm (Alku 1992) . Returns sample-rate glottal flow and its
derivative (MFDR), which are used for voice quality analysis and
dysphonia assessment. Prefer this over GFM-IAIF when a sample-domain
output at the native audio rate is needed.

## Usage

``` r
trk_covarep_iaif(listOfFiles, ...)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- ...:

  Additional arguments (currently unused).

- order_vt:

  Integer. Vocal tract LPC order. `NULL` (default) sets it automatically
  to `2 * round(fs/2000) + 4` (typically 12–20).

- order_gl:

  Integer. Glottal source LPC order. `NULL` (default) sets it
  automatically to `2 * round(fs/4000)` (typically 2–4).

- leaky_coef:

  Numeric. Leaky integration coefficient for lip radiation compensation
  (0 \< leaky_coef \< 1). Default 0.99.

- hpfilt:

  Logical. Apply high-pass filter before processing. Default `TRUE`.

- toFile:

  Logical. If `TRUE`, write SSFF output files and return file paths. If
  `FALSE`, return an `AsspDataObj`. Default `TRUE`.

- explicitExt:

  Character. Output file extension. Default `"glf"`.

## Value

If `toFile = FALSE`: an `AsspDataObj` with tracks:

- `glottal_flow`:

  REAL64, estimated glottal flow waveform, dimensionless, n_samples × 1.
  Same length as the input audio.

- `glottal_derivative`:

  REAL64, glottal flow derivative (MFDR proxy), dimensionless, n_samples
  × 1.

Frame rate equals the audio sample rate. If `toFile = TRUE`: character
vector of output file paths.

## Details

Default LPC orders follow the original COVAREP/MATLAB implementation.
`order_gl = 3` is recommended; departing from this may yield unreliable
source estimates.

## References

(Alku 1992)

## See also

[`lst_covarep_vq`](https://humlab-speech.github.io/superassp/reference/lst_covarep_vq.md)
for voice quality parameters,
[`trk_gfmiaif`](https://humlab-speech.github.io/superassp/reference/trk_gfmiaif.md)
for GFM-IAIF (LP-coefficient output)

## Examples

``` r
if (FALSE) { # \dontrun{
# Single file - extract glottal flow
glottal <- trk_covarep_iaif("vowel.wav", toFile = FALSE)
plot(glottal$glottal_flow, type = "l", ylab = "Glottal Flow")
plot(glottal$glottal_derivative, type = "l", ylab = "Glottal Derivative")

# Batch processing
files <- c("a.wav", "e.wav", "i.wav")
result <- trk_covarep_iaif(files, toFile = TRUE)

# Custom filter orders
glottal_custom <- trk_covarep_iaif("audio.wav",
                                   order_vt = 16,
                                   order_gl = 3,
                                   toFile = FALSE)
} # }
```
