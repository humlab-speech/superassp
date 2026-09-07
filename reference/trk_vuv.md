# Voiced/unvoiced segmentation via two-pass adaptive pitch detection

Classifies each frame as voiced (V) or unvoiced (U) using a two-pass,
speaker-adaptive pitch detection strategy (Al-Tamimi & Khattab 2015,
2018). Output is either a Praat TextGrid with a VUV interval tier or an
SSFF binary `voicing` track (0 = unvoiced, 1 = voiced).

## Usage

``` r
trk_vuv(listOfFiles, ...)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- timeStep:

  Numeric. Frame shift for pitch analysis in seconds. Also sets the SSFF
  output frame rate (1 / timeStep Hz) when `outputFormat = "ssff"`.
  Default 0.005 s (200 Hz).

- initialMinPitch:

  Numeric. Lower F0 bound for the first-pass pitch estimate in Hz.
  Default 50 Hz.

- initialMaxPitch:

  Numeric. Upper F0 bound for the first-pass pitch estimate in Hz.
  Default 800 Hz.

- voicingThreshold:

  Numeric. Minimum pitch strength for a frame to be voiced (0–1).
  Default 0.45.

- vuvMaxPeriod:

  Numeric. Maximum glottal period (s) for TextGrid VUV conversion.
  Default 0.02 s (minimum 50 Hz).

- minPeriod:

  Numeric. Minimum glottal period (s) accepted when computing mean
  period from the PointProcess. Default 0.0001 s.

- maxPeriod:

  Numeric. Maximum glottal period (s) accepted when computing mean
  period from the PointProcess. Default 0.02 s.

- maxPeriodFactor:

  Numeric. Maximum ratio between consecutive periods still treated as
  periodic. Default 1.3.

- windowShape:

  Character. Window shape applied to audio before loading. Default
  `"Gaussian1"`.

- relativeWidth:

  Numeric. Relative width of the window. Default 1.0.

- outputFormat:

  Character. Output format: `"textgrid"` (Praat TextGrid with interval
  tier) or `"ssff"` (binary AsspDataObj track). Default `"textgrid"`.

- toFile:

  Logical. If `TRUE`, write output files and return paths (invisibly).
  If `FALSE`, return the in-memory object. Default `TRUE`.

- explicitExt:

  Character. Output file extension. Defaults to `"TextGrid"` when
  `outputFormat = "textgrid"` and `"vuv"` when `"ssff"`.

## Value

Depends on `outputFormat` and `toFile`:

- `outputFormat = "textgrid"`, `toFile = FALSE`:

  A pladdrr TextGrid object (or list) with one interval tier containing
  V/U labels.

- `outputFormat = "textgrid"`, `toFile = TRUE`:

  Character vector of output TextGrid paths, returned invisibly.

- `outputFormat = "ssff"`, `toFile = FALSE`:

  An `AsspDataObj` with track `voicing` (INT16, binary 0/1, n_frames x
  1). Frame rate: `1 / timeStep` Hz (default 200 Hz).

- `outputFormat = "ssff"`, `toFile = TRUE`:

  Character vector of output SSFF file paths, returned invisibly.

## Details

Pass 1: pitch estimated across `initialMinPitch`–`initialMaxPitch` Hz on
a 0–500 Hz bandpass-filtered signal. Pass 2: adaptive bounds set to Q1 x
0.75 and Q3 x 1.5 of the voiced frames from pass 1. The PointProcess
derived from the refined pitch drives TextGrid VUV interval creation.

## References

(Al-Tamimi and Khattab 2015)

(Al-Tamimi and Khattab 2018)

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate VUV TextGrid
trk_vuv("speech.wav", toFile = TRUE)
# Creates speech.TextGrid with VUV tier

# Get TextGrid object in memory
tg <- trk_vuv("speech.wav", toFile = FALSE)

# Generate binary SSFF track instead
result <- trk_vuv("speech.wav", outputFormat = "ssff", toFile = FALSE)
plot(result$voicing, type = "l")

# Custom pitch range
trk_vuv("speech.wav", initialMinPitch = 75, initialMaxPitch = 500)
} # }
```
