# Track fundamental frequency using the YIN algorithm

Extracts F0 using the YIN algorithm (Cheveigné and Kawahara 2002) ,
which applies cumulative mean normalized difference (CMND) to the
autocorrelation function for reliable pitch detection even in moderate
noise. Returns both F0 and a per-frame voicing probability derived from
the CMND minimum. For better noise robustness, prefer the probabilistic
extension
[`trk_pitch_pyin`](https://humlab-speech.github.io/superassp/reference/trk_pitch_pyin.md).

## Usage

``` r
trk_pitch_yin(listOfFiles, ...)
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
  windowShift Hz). Default 5.0 ms.

- windowSize:

  Numeric. Analysis window length in milliseconds. Default 30.0 ms.

- minF:

  Numeric. Minimum F0 in Hz. Default 70.0 Hz.

- maxF:

  Numeric. Maximum F0 in Hz. Default 200.0 Hz.

- threshold:

  Numeric. YIN CMND threshold (0–1); lower values accept weaker pitch
  candidates (more voiced frames). Default 0.1.

- toFile:

  Logical. If `TRUE`, write SSFF output files and return the paths
  written invisibly. If `FALSE`, return an `AsspDataObj`. Default
  `FALSE`.

- explicitExt:

  Character. Output file extension. Default `"yip"`.

- outputDirectory:

  Character. Directory for output files. `NULL` (default) writes
  alongside the input file.

- verbose:

  Logical. Print per-file progress. Default `TRUE`.

## Value

If `toFile = FALSE`: an `AsspDataObj` with tracks:

- `F0`:

  REAL32, fundamental frequency in Hz, n_frames × 1. Zero indicates
  unvoiced frames.

- `prob`:

  REAL32, voicing probability, 0–1, n_frames × 1.

Frame rate: `1000 / windowShift` Hz (default 200 Hz). If
`toFile = TRUE`: output file path(s), returned invisibly.

## References

Cheveigné Ad, Kawahara H (2002). “YIN, a fundamental frequency estimator
for speech and music.” *The Journal of the Acoustical Society of
America*, **111**(4), 1917–1930. ISSN 0001-4966.
[doi:10.1121/1.1458024](https://doi.org/10.1121/1.1458024) .

## Examples

``` r
if (FALSE) { # \dontrun{
# Extract F0 from audio file
f0_data <- trk_pitch_yin("recording.wav", toFile = FALSE)

# Process with custom parameters
trk_pitch_yin("speech.mp3", minF = 75, maxF = 300, windowShift = 10)

# Process video file (extracts audio)
trk_pitch_yin("interview.mp4", toFile = FALSE)
} # }
```
