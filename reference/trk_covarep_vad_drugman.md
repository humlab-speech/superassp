# Detect voiced frames using Drugman's multi-branch VAD

Estimates per-frame speech activity as posterior probabilities using
three independent ANN classifiers (MFCC-based, Sadjadi pitch-related,
and CPP/SRH features) combined by geometric mean (Drugman et al. 2012) .
Suitable for pre-filtering frames before pitch or voice quality
analysis.

## Usage

``` r
trk_covarep_vad_drugman(listOfFiles, ...)
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

  Character. Output file extension. Default `"cvd"`.

## Value

If `toFile = FALSE`: an `AsspDataObj` with tracks:

- `vad_final`:

  FLOAT, ensemble voicing posterior (geometric mean of three branches),
  0–1, n_frames × 1. Recommended for downstream use.

- `vad_mfcc`:

  FLOAT, MFCC-branch posterior, 0–1, n_frames × 1.

- `vad_sadjadi`:

  FLOAT, Sadjadi pitch-feature posterior, 0–1, n_frames × 1.

- `vad_new`:

  FLOAT, CPP/SRH-feature posterior, 0–1, n_frames × 1.

Frame rate: 200 Hz (5 ms hop, interpolated from 10 ms internal hop). If
`toFile = TRUE`: character vector of output file paths, returned
invisibly.

## Details

Audio is resampled to 16 kHz internally. A threshold of 0.5 on
`vad_final` gives a reasonable binary voiced/unvoiced decision; 0.7 is
more conservative. Each branch applies an 11-frame median filter to the
posterior before combination.

## References

Drugman T, Soria-Olivas E, Perez-Córdoba JL, Alwan A (2012). “A
comparative study of different feature sets for acoustic voice quality
assessment.” *IEEE Transactions on Audio, Speech, and Language
Processing*, **20**(6), 1690–1703.
[doi:10.1109/TASL.2012.2188377](https://doi.org/10.1109/TASL.2012.2188377)
. Multi-branch voice activity detection using MFCC and Sadjadi features.

## Examples

``` r
if (FALSE) { # \dontrun{
trk_covarep_vad_drugman(
  system.file("samples", "sustained", "a1.wav", package = "superassp"),
  toFile = FALSE
)
} # }
```
