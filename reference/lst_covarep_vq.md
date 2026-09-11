# Extract voice quality parameters (NAQ, QOQ, H1-H2, HRF, PSP) per utterance

Computes five IAIF-derived glottal voice quality measures (NAQ, QOQ,
H1-H2, HRF, PSP) summarised per utterance. For frame-by-frame
trajectories see
[`trk_covarep_vq_gci`](https://humlab-speech.github.io/superassp/reference/trk_covarep_vq_gci.md).

## Usage

``` r
lst_covarep_vq(
  listOfFiles,
  beginTime = 0,
  endTime = 0,
  f0 = NULL,
  gci = NULL,
  gci_in_samples = FALSE,
  verbose = TRUE,
  toFile = FALSE,
  explicitExt = "cvq",
  outputDirectory = NULL
)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- f0:

  Optional F0 estimate. Can be:

  - Scalar: Single F0 value for entire signal

  - Vector: F0 contour (median of voiced frames will be used)

  - NULL: H1-H2 will not be computed (returns NA)

- gci:

  Optional glottal closure instants. Can be:

  - Numeric vector of sample indices

  - Numeric vector of times in seconds (will be converted to samples)

  - NULL: NAQ and QOQ will not be computed (return NA)

- gci_in_samples:

  Logical; if TRUE, gci is in sample indices; if FALSE, gci is in
  seconds (default: FALSE)

- toFile:

  Logical. If TRUE, write results to JSTF file. Default FALSE.

- explicitExt:

  Character. File extension for output. Default "cvq".

- beginTime:

  Start time for the extracted portion in seconds. Default: NULL
  (beginning of signal). Note: uses `beginTime`/`endTime` (seconds)
  matching DSP function conventions, unlike
  [`read_audio()`](https://humlab-speech.github.io/superassp/reference/read_audio.md)
  which uses `begin`/`end`.

- endTime:

  The end time of the section of the sound files that should be analysed
  (in seconds). Use 0 for end of file.

- verbose:

  Logical. Show a progress bar (sequential path) or a progress-aware
  parallel apply (`pbapply`/`pbmcapply`, if installed).

- outputDirectory:

  The directory where the slice file should be stored. If not defiled
  (NULL), the sparse slice file will placed in the same folder as the
  media file.

## Value

If `toFile=FALSE` (default), for single file: Named list with voice
quality parameters. For multiple files: List of named lists. If
`toFile=TRUE`, invisibly returns the path(s) to the written JSTF
file(s).

Each parameter list contains:

- `glottal_flow_max`:

  Peak glottal flow amplitude

- `glottal_flow_min`:

  Minimum glottal flow value

- `glottal_derivative_peak`:

  Maximum flow derivative (MFDR)

- `NAQ`:

  Normalized Amplitude Quotient (NA if no GCI)

- `QOQ`:

  Quasi-Open Quotient (NA if no GCI)

- `H1_H2`:

  First two harmonics difference in dB (NA if no F0)

- `HRF`:

  Harmonic Richness Factor

- `PSP`:

  Parabolic Spectral Parameter

## Details

Internally applies IAIF (Iterative Adaptive Inverse Filtering, (Alku
1992) ) to extract the glottal source, then computes GCI-anchored
measures from the glottal flow and its derivative. Implemented in native
C++ (no Python dependency).

## References

(Alku 1992)

(Kane and Gobl 2013)

## See also

[`trk_covarep_iaif`](https://humlab-speech.github.io/superassp/reference/trk_covarep_iaif.md)
for glottal waveforms,
[`trk_pitch_srh`](https://humlab-speech.github.io/superassp/reference/trk_pitch_srh.md)
for F0 estimation

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage (HRF and PSP only, no F0 or GCI)
vq <- lst_covarep_vq("vowel.wav")
print(vq$HRF)
print(vq$PSP)

# With F0 for H1-H2 computation
vq <- lst_covarep_vq("vowel.wav", f0 = 150)
print(vq$H1_H2)

# Batch processing
files <- c("a.wav", "e.wav", "i.wav", "o.wav", "u.wav")
vq_all <- lst_covarep_vq(files)
} # }
```
