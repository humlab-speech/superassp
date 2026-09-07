# Track fundamental frequency using the Snack/ESPS dp_f0 algorithm

Extracts F0, voicing probability, RMS energy, and autocorrelation peak
using the Snack Sound Toolkit normalized cross-correlation +
dynamic-programming pitch tracker (`dp_f0`, Talkin 1995). Unlike
[`trk_pitch_rapt`](https://humlab-speech.github.io/superassp/reference/trk_pitch_rapt.md),
which returns only F0, this function exposes all four tracks for
downstream signal quality assessment.

## Usage

``` r
trk_pitch_snack(listOfFiles, ...)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- beginTime:

  Numeric. Start of analysis window in seconds. Default 0 (file start).

- endTime:

  Numeric. End of analysis window in seconds. Default 0 (file end).

- windowShift:

  Numeric. Frame shift in milliseconds; sets output frame rate (1000 /
  windowShift Hz). Default 10.0 ms.

- minF:

  Numeric. Minimum F0 in Hz. Default 50.0 Hz.

- maxF:

  Numeric. Maximum F0 in Hz. Default 550.0 Hz.

- voiceBias:

  Numeric. Bias toward the voiced hypothesis in the DP cost function
  (range approximately −0.5 to 0.5; positive = more voiced frames).
  Default 0.0.

- toFile:

  Logical. If `TRUE`, write SSFF output files and return the count
  written invisibly. If `FALSE`, return an `AsspDataObj`. Default
  `TRUE`.

- explicitExt:

  Character. Output file extension. Default `"snackpitch"`.

- outputDirectory:

  Character. Directory for output files. `NULL` (default) writes
  alongside the input file.

- verbose:

  Logical. Print per-file progress. Default `TRUE`.

## Value

If `toFile = FALSE`: an `AsspDataObj` with tracks:

- `f0`:

  REAL32, fundamental frequency in Hz, n_frames × 1. Zero indicates
  unvoiced frames.

- `voicing`:

  REAL32, voicing probability, 0–1, n_frames × 1.

- `rms`:

  REAL32, RMS energy (linear), dimensionless, n_frames × 1.

- `acpeak`:

  REAL32, autocorrelation peak magnitude, 0–1, n_frames × 1.

Frame rate: `1000 / windowShift` Hz (default 100 Hz). If
`toFile = TRUE`: integer count of files written, returned invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
# Full 4-track pitch analysis
res <- trk_pitch_snack("recording.wav", toFile = FALSE)
names(res)  # "f0" "voicing" "rms" "acpeak"

# Custom F0 range
trk_pitch_snack("speech.mp3", minF = 75, maxF = 300)
} # }
```
