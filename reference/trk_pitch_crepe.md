# Track fundamental frequency and periodicity using CREPE (ONNX)

Estimates F0 by applying a deep convolutional neural network directly to
the time-domain waveform (CREPE; (Kim et al. 2018) ). Two model sizes
are available: `"tiny"` (~1.9 MB, fast) and `"full"` (~85 MB, more
accurate). No Python or PyTorch required — inference uses the bundled
ONNX Runtime.

## Usage

``` r
trk_pitch_crepe(
  listOfFiles,
  beginTime = 0,
  endTime = 0,
  windowShift = 10,
  windowSize = 15,
  minF = 50,
  maxF = 550,
  voicing.threshold = 0.21,
  silence.threshold = -60,
  model = c("tiny", "full"),
  decoder = c("viterbi", "argmax"),
  batch_size = 512L,
  explicitExt = "crp",
  outputDirectory = NULL,
  toFile = TRUE,
  verbose = TRUE
)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted.

- beginTime:

  Numeric. Start of analysis window in seconds. Default 0 (file start).

- endTime:

  Numeric. End of analysis window in seconds. Default 0 (file end).

- windowShift:

  Numeric. Frame shift in milliseconds; sets output frame rate
  (`1000 / windowShift` Hz). Default 10 ms (100 Hz).

- windowSize:

  Numeric. Smoothing filter window size in milliseconds, applied to both
  median (periodicity) and mean (F0) post-processing filters. Default 15
  ms.

- minF:

  Numeric. Minimum F0 in Hz. Frames with estimated F0 below this value
  are set to 0 (unvoiced). Default 50 Hz.

- maxF:

  Numeric. Maximum F0 in Hz. Frames above this value are set to 0.
  Default 550 Hz.

- voicing.threshold:

  Numeric. Periodicity threshold for voicing decisions. Frames with
  periodicity below this value are set to unvoiced (F0 = 0). Default
  0.21.

- silence.threshold:

  Numeric. A-weighted dB threshold relative to global maximum. Frames
  below this are treated as silent (periodicity = 0). Default -60 dB.

- model:

  Character. Model variant: `"tiny"` (fast) or `"full"` (more accurate).
  Default `"tiny"`.

- decoder:

  Character. Decoding method: `"viterbi"` (recommended, smoother) or
  `"argmax"` (faster but noisier). Default `"viterbi"`.

- batch_size:

  Integer. Number of frames per ONNX inference batch. Default 512.

- explicitExt:

  Character. Output file extension. Default `"crp"`.

- outputDirectory:

  Character. Directory for output files. `NULL` (default) writes
  alongside the input file.

- toFile:

  Logical. If `TRUE`, write SSFF output files and return the count
  written (invisibly). If `FALSE`, return an `AsspDataObj`. Default
  `TRUE`.

- verbose:

  Logical. Print per-file progress. Default `TRUE`.

## Value

If `toFile = FALSE`: an `AsspDataObj` with tracks:

- `f0`:

  REAL32, Hz, *n_frames* × 1. Fundamental frequency; 0 in
  unvoiced/silent frames.

- `periodicity`:

  REAL32, 0–1, *n_frames* × 1. Model confidence; values below
  `voicing.threshold` are treated as unvoiced.

Frame rate: `1000 / windowShift` Hz (default 100 Hz, 10 ms hop). If
`toFile = TRUE`: integer count of files written, returned invisibly.

## Details

ONNX Runtime is installed automatically on first use (~30 MB, cached in
the R user directory) and persists across R sessions and package
reinstalls. The model file (`tiny`: ~2 MB, `full`: ~85 MB) is downloaded
from the [torchcrepe-onnx Hugging Face Hub
repo](https://huggingface.co/FredrikKarlssonSpeech/torchcrepe-onnx) on
first use (requires the huggingfaceR package and a network connection)
and cached in the R user directory; subsequent calls read the cached
copy with no network access.

Post-processing (matching torchcrepe): median filter on periodicity →
A-weighted silence detection → voicing threshold → NaN-aware mean filter
on F0 → F0 range clamping.

## References

Kim JW, Salamon J, Li P, Bello JP (2018). “Crepe: A Convolutional
Representation for Pitch Estimation.” *2018 IEEE International
Conference on Acoustics, Speech and Signal Processing (ICASSP)*, **00**,
161–165.
[doi:10.1109/icassp.2018.8461329](https://doi.org/10.1109/icassp.2018.8461329)
.

## Examples

``` r
if (FALSE) { # \dontrun{
trk_pitch_crepe(
  system.file("samples", "sustained", "a1.wav", package = "superassp"),
  toFile = FALSE
)
} # }
```
