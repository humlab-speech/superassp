# Decompose speech into vocal tract, glottis, and lip radiation LP filters (GFM-IAIF)

Returns per-frame LP coefficient tracks for the vocal tract, glottis,
and lip radiation filters using GFM-IAIF source-filter separation
(Perrotin and McLoughlin 2019) . Prefer this over
[`trk_covarep_iaif()`](https://humlab-speech.github.io/superassp/reference/trk_covarep_iaif.md)
when you need frame-level LP coefficients rather than a sample-domain
glottal waveform.

## Usage

``` r
trk_gfmiaif(listOfFiles, beginTime = 0, centerTime = FALSE, endTime = 0, windowShift = 10, windowSize = 32, nv = 48L, ng = 3L, d = 0.99, window = "HANN", explicitExt = "gfm", outputDirectory = NULL, toFile = TRUE, verbose = TRUE)
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

Perrotin O, McLoughlin I (2019). “A Spectral Glottal Flow Model for
Source-filter Separation of Speech.” In *ICASSP 2019 - 2019 IEEE
International Conference on Acoustics, Speech and Signal Processing
(ICASSP)*, 7160–7164.
[doi:10.1109/ICASSP.2019.8682625](https://doi.org/10.1109/ICASSP.2019.8682625)
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
