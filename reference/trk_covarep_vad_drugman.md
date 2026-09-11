# Detect voiced frames using Drugman's multi-branch VAD

Estimates per-frame speech activity as posterior probabilities using
three independent ANN classifiers (MFCC-based, Sadjadi pitch-related,
and CPP/SRH features) combined by geometric mean (Drugman et al. 2016) .
Suitable for pre-filtering frames before pitch or voice quality
analysis.

## Usage

``` r
trk_covarep_vad_drugman(
  listOfFiles,
  beginTime = 0,
  endTime = 0,
  toFile = FALSE,
  explicitExt = "cvd",
  outputDirectory = NULL,
  verbose = TRUE
)
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

Drugman T, Stylianou Y, Kida Y, Akamine M (2016). “Voice Activity
Detection: Merging Source and Filter-based Information.” *IEEE Signal
Processing Letters*, **23**(2), 252–256.
[doi:10.1109/LSP.2015.2495219](https://doi.org/10.1109/LSP.2015.2495219)
. Multi-branch voice activity detection merging MFCC (filter) and
Sadjadi/SRH (source) features.

## Examples

``` r
if (FALSE) { # \dontrun{
trk_covarep_vad_drugman(
  system.file("samples", "sustained", "a1.wav", package = "superassp"),
  toFile = FALSE
)
} # }
```
