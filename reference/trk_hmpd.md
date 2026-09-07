# Extract harmonic model phase distortion features (HMPD): AE, PDM, PDD

Returns amplitude envelope (AE), phase deviation mean (PDM), and phase
deviation deviation (PDD) per frame using the HMPD method (Degottex et
al. 2014) . Suitable for voice quality characterization and
vocoder-based resynthesis.

## Usage

``` r
trk_hmpd(listOfFiles, ...)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- f0s:

  Optional numeric matrix with 2 columns (time in seconds, F0 in Hz). If
  `NULL` (default), a constant 100 Hz F0 is assumed with a warning;
  supply an accurate F0 track for reliable PDM/PDD values.

- f0min:

  Numeric. Minimum F0 in Hz for sinusoidal analysis. Default 60 Hz.

- f0max:

  Numeric. Maximum F0 in Hz for sinusoidal analysis. Default 440 Hz.

- toFile:

  Logical. If `TRUE`, write SSFF output files and return the paths
  written invisibly. If `FALSE`, return an `AsspDataObj`. Default
  `FALSE`.

- explicitExt:

  Character. Output file extension. Default `"hpd"`.

## Value

If `toFile = FALSE`: an `AsspDataObj` with tracks:

- `ae`:

  FLOAT, Mel-frequency-warped cepstral amplitude envelope, n_frames × 25
  (coefficients 0–24). Dimensionless (log-amplitude units).

- `pdm`:

  FLOAT, phase deviation mean on a log-harmonic scale, n_frames × 25.
  Units: radians.

- `pdd`:

  FLOAT, phase deviation deviation (aperiodicity), n_frames × 13. Larger
  values indicate more noise-like phase scatter.

Frame rate: 100 Hz (10 ms grid, resampled from variable-rate sinusoidal
analysis). If `toFile = TRUE`: character vector of output file paths,
returned invisibly.

## Details

HMPD captures both spectral shape (AE) and aperiodicity-related phase
scatter (PDM, PDD). Phase deviations are computed relative to a
deterministic sinusoidal model at each harmonic; the mean and standard
deviation across harmonics give PDM and PDD.

When `f0s` is `NULL` the function falls back to a constant 100 Hz F0 and
emits a warning. AE is computed by Mel-cepstrum (method 2 in COVAREP).
PDM uses a log-harmonic scale with 8 harmonics in the linear region
(Bezier interpolation above). PDD is derived from circular phase
variance.

## References

Degottex G, Kane J, Drugman T, Raitio T, Scherer S (2014). “COVAREP: a
collaborative voice analysis repository for speech technologies.” In
*Proceedings of the 2014 IEEE International Conference on Acoustics,
Speech and Signal Processing (ICASSP)*, 960–964.
[doi:10.1109/ICASSP.2014.6853739](https://doi.org/10.1109/ICASSP.2014.6853739)
. Open-source toolkit for voice analysis including HMPD, GCI detection,
formants, pitch, and voice quality measures,
<https://github.com/covarep/covarep>.

## Examples

``` r
if (FALSE) { # \dontrun{
# Single file, return object
hmpd <- trk_hmpd("speech.wav", toFile = FALSE)

# Batch process multiple files
files <- c("file1.wav", "file2.wav")
hmpd_results <- trk_hmpd(files, toFile = TRUE, outputDirectory = "output/")
} # }
```
