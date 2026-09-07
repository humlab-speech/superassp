# Estimate band aperiodicity using the D4C algorithm (WORLD vocoder)

Computes per-frame, per-band aperiodicity using the D4C (Death,
Destruction, Diversion and Disgrace) estimator from the WORLD vocoder
via SPTK. Aperiodicity quantifies noise-to-harmonic energy ratio across
frequency bands and is the primary input to WORLD's noise excitation
model. Use alongside
[`trk_pitch_rapt()`](https://humlab-speech.github.io/superassp/reference/trk_pitch_rapt.md)
or
[`trk_pitch_swipe()`](https://humlab-speech.github.io/superassp/reference/trk_pitch_swipe.md)
for full WORLD vocoder analysis/resynthesis.

## Usage

``` r
trk_d4c(listOfFiles, ...)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- windowShift:

  Numeric. Frame shift in milliseconds; sets output frame rate (1000 /
  windowShift Hz). Default 5.0 ms.

- minF:

  Numeric. Minimum F0 in Hz used for the internal pitch estimator.
  Default 60.0 Hz.

- maxF:

  Numeric. Maximum F0 in Hz used for the internal pitch estimator.
  Default 400.0 Hz.

- voicing_threshold:

  Numeric. Voicing threshold for the internal F0 detector (0–1; higher =
  more conservative). Default 0.85.

- threshold:

  Numeric. D4C aperiodicity clipping threshold (0–1). Default 0.85.

- toFile:

  Logical. If `TRUE`, write SSFF output files and return the count
  written invisibly. If `FALSE`, return an `AsspDataObj`. Default
  `TRUE`.

- explicitExt:

  Character. Output file extension. Default `"ap"`.

## Value

If `toFile = FALSE`: an `AsspDataObj` with track:

- `aperiodicity`:

  REAL32, band aperiodicity, 0–1 (0 = fully periodic, 1 = fully
  aperiodic), n_frames × n_bands where n_bands = floor(fs/2 / 3000) + 1.

Frame rate: `1000 / windowShift` Hz (default 200 Hz). If
`toFile = TRUE`: integer count of files written, returned invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
# Estimate aperiodicity
trk_d4c("recording.wav")

# Process with custom parameters
trk_d4c("speech.wav", minF = 80, maxF = 350, windowShift = 10)

# Process multiple files
trk_d4c(c("file1.wav", "file2.wav"))
} # }
```
