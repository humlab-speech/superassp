# Track fundamental frequency using SwiftF0 (ONNX)

Estimates F0 by applying a convolutional neural network directly to a
short-time spectrogram computed inside the ONNX graph (SwiftF0;
(Nieradzik 2025) ). SwiftF0 targets real-time speed (~130 ms for 5 s of
audio on CPU) at accuracy competitive with CREPE. No Python or PyTorch
required — inference uses ONNX Runtime.

## Usage

``` r
trk_pitch_swiftf0(listOfFiles, ...)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted.

- beginTime:

  Numeric. Start of analysis window in seconds. Default 0 (file start).

- endTime:

  Numeric. End of analysis window in seconds. Default 0 (file end).

- minF:

  Numeric. Minimum F0 in Hz to treat as voiced. Default 75 Hz (speech).
  Must be \>= 46.875 Hz (model minimum; G1). For music, use 46.875.

- maxF:

  Numeric. Maximum F0 in Hz to treat as voiced. Default 400 Hz (speech).
  Must be \<= 2093.75 Hz (model maximum; C7). For music, use 2093.75.

- confidence_threshold:

  Numeric (0–1). Frames with model confidence below this value are
  marked unvoiced. Default 0.9.

- toFile:

  Logical. If `TRUE`, write SSFF output files and return the count
  written (invisibly). If `FALSE`, return an `AsspDataObj`. Default
  `TRUE`.

- explicitExt:

  Character. Output file extension. Default `"sf0"`.

- outputDirectory:

  Character. Directory for output files. `NULL` (default) writes
  alongside the input file.

- verbose:

  Logical. Print per-file progress. Default `TRUE`.

## Value

If `toFile = FALSE`: an `AsspDataObj` with tracks:

- `f0`:

  REAL32, Hz, *n\\frames* × 1. Fundamental frequency; 0 in unvoiced
  frames.

- `confidence`:

  REAL32, 0–1, *n\\frames* × 1. Model voicing confidence.

Frame rate: fixed 62.5 Hz (16 ms hop; not configurable — fixed by the
ONNX graph's internal STFT). If `toFile = TRUE`: integer count of files
written, returned invisibly.

## Details

ONNX Runtime is installed automatically on first use (~30 MB, cached in
the R user directory). The model file (~400 KB) is downloaded from the
[swift-f0-onnx Hugging Face Hub
repo](https://huggingface.co/FredrikKarlssonSpeech/swift-f0-onnx) on
first use (requires the huggingfaceR package and a network connection)
and cached in the R user directory; subsequent calls read the cached
copy with no network access.

Pre-processing (fixed by model training): resample to 16 kHz mono,
normalise to full-scale float32 in `[-1, 1]`. Framing, STFT (1024-sample
window, 256-sample hop), and pitch/confidence estimation all happen
inside the ONNX graph — no manual windowing is needed. Post-processing
(applied outside the graph, matching the reference implementation): a
frame is voiced when `confidence > confidence_threshold` AND its pitch
estimate falls within `[minF, maxF]`; unvoiced frames get `f0 = 0`.

## References

Nieradzik L (2025). “SwiftF0: Fast and Accurate Monophonic Pitch
Detection.” Swift-F0 is a fast and accurate F0 detector using CNN-based
approach on STFT spectrograms, achieving 132ms processing time for 5
seconds of audio on CPU, 2508.18440, <https://arxiv.org/abs/2508.18440>.

## Examples

``` r
if (FALSE) { # \dontrun{
trk_pitch_swiftf0(
  system.file("samples", "sustained", "a1.wav", package = "superassp"),
  toFile = FALSE
)
} # }
```
