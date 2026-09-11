# Track fundamental frequency using probabilistic YIN (pYIN)

Extracts F0 and a per-frame voicing probability using a simplified C++
implementation of the probabilistic YIN algorithm (Mauch and Dixon 2014)
, which extends YIN (Cheveigné and Kawahara 2002) by evaluating multiple
pitch candidates. More robust than plain YIN in noise and lower-voiced
speech; the `prob` track can serve as a soft voicing mask.

## Usage

``` r
trk_pitch_pyin(
  listOfFiles,
  beginTime = 0,
  endTime = 0,
  windowShift = 5,
  windowSize = 30,
  minF = 70,
  maxF = 200,
  threshold = 0.1,
  toFile = FALSE,
  explicitExt = "pyp",
  outputDirectory = NULL,
  verbose = TRUE
)
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

  Numeric. YIN dip threshold (0–1); lower values accept weaker pitch
  candidates (more voiced frames). Default 0.1.

- toFile:

  Logical. If `TRUE`, write SSFF output files and return the paths
  written invisibly. If `FALSE`, return an `AsspDataObj`. Default
  `FALSE`.

- explicitExt:

  Character. Output file extension. Default `"pyp"`.

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
  
Mauch M, Dixon S (2014). “PYIN: A Fundamental Frequency Estimator using
Probabilistic Threshold Distributions.” *2014 IEEE International
Conference on Acoustics, Speech and Signal Processing (ICASSP)*,
659–663.
[doi:10.1109/icassp.2014.6853678](https://doi.org/10.1109/icassp.2014.6853678)
.

## Examples

``` r
if (FALSE) { # \dontrun{
# Extract F0 from audio file
f0_data <- trk_pitch_pyin("recording.wav", toFile = FALSE)

# Process with custom parameters
trk_pitch_pyin("speech.mp3", minF = 75, maxF = 300, windowShift = 10)

# Process video file (extracts audio)
trk_pitch_pyin("interview.mp4", toFile = FALSE)
} # }
```
