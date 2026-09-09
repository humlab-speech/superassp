# Track pitch and voiced speech using the TANDEM-STRAIGHT algorithm

Estimates F0 and per-frame voicing probability using a gammatone
filterbank combined with neural network-based pitch tracking (Hu & Wang
2010), which simultaneously segregates voiced speech from noise. TANDEM
is robust to noise and reverberation and can track multiple simultaneous
pitch sources.

## Usage

``` r
trk_tandem(listOfFiles, minF = 50, maxF = 500, target_sample_rate = 20000, return_mask = FALSE, toFile = FALSE, explicitExt = "tnd", outputDirectory = NULL, verbose = TRUE, ...)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; audio is resampled to `target_sample_rate` Hz internally.

- ...:

  Additional arguments (currently unused).

- minF:

  Numeric. Minimum F0 in Hz. Default 50 Hz.

- maxF:

  Numeric. Maximum F0 in Hz. Default 500 Hz.

- target_sample_rate:

  Numeric. Internal processing sample rate in Hz. TANDEM requires 20000
  Hz. Default 20000.

- return_mask:

  Logical. Return time-frequency voiced mask (currently unused). Default
  `FALSE`.

- toFile:

  Logical. If `TRUE`, write SSFF output files and return the paths
  written. If `FALSE`, return an `AsspDataObj`. Default `FALSE`.

- explicitExt:

  Character. Output file extension. Default `"tnd"`.

- outputDirectory:

  Character. Directory for output files. `NULL` (default) writes
  alongside the input file.

- verbose:

  Logical. Print per-file progress. Default `TRUE`.

## Value

If `toFile = FALSE`: an `AsspDataObj` with tracks:

- `pitch`:

  REAL64, fundamental frequency in Hz, n_frames × 1. Zero indicates
  unvoiced frames.

- `voicing_prob`:

  REAL64, voicing probability, 0–1, n_frames × 1.

Frame rate: 100 Hz (fixed 10 ms hop). If `toFile = TRUE`: character
vector of output file paths.

## Note

The core processing is currently a placeholder; full TANDEM C++
integration is under development. Results reflect the algorithm
framework but may not match the published TANDEM-STRAIGHT output.

## References

Hu G, Wang D (2010). “A tandem algorithm for pitch estimation and voiced
speech segregation.” *IEEE Transactions on Audio, Speech, and Language
Processing*, **18**(8), 2067–2079.
[doi:10.1109/TASL.2010.2041110](https://doi.org/10.1109/TASL.2010.2041110)
.

Hu K, Wang D (2011). “Unvoiced speech segregation from nonspeech
interference via CASA and spectral subtraction.” *IEEE Transactions on
Audio, Speech, and Language Processing*, **19**(6), 1600–1609.
[doi:10.1109/TASL.2010.2093893](https://doi.org/10.1109/TASL.2010.2093893)
.

## See also

[`trk_pitch_rapt`](https://humlab-speech.github.io/superassp/reference/trk_pitch_rapt.md),
[`trk_pitch_swipe`](https://humlab-speech.github.io/superassp/reference/trk_pitch_swipe.md),
[`trk_pitch_yin`](https://humlab-speech.github.io/superassp/reference/trk_pitch_yin.md)
for other pitch tracking methods

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic pitch tracking
result <- trk_tandem("speech.wav")
plot(result$pitch, type = "l", main = "TANDEM Pitch Track")

# With noisy speech
result <- trk_tandem("noisy_speech.wav", minF = 80, maxF = 400)

# Batch processing
files <- c("speaker1.wav", "speaker2.wav", "speaker3.wav")
results <- trk_tandem(files, verbose = TRUE)

# Save to files
trk_tandem("speech.wav", toFile = TRUE, outputDirectory = "output/")
} # }
```
