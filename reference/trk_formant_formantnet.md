# Track formant frequencies and bandwidths using FormantNet (ONNX)

Predicts formant frequencies and bandwidths from log-scale smoothed
spectral envelopes using a bidirectional LSTM (FormantNet; (Sakamoto et
al. 2021) ). The model was trained on 16 kHz speech; all audio is
resampled automatically. No Python or TensorFlow required — inference
uses ONNX Runtime.

## Usage

``` r
trk_formant_formantnet(listOfFiles, beginTime = 0, endTime = 0, numFormants = 3L, windowShift = 5, toFile = TRUE, explicitExt = "fnf", outputDirectory = NULL, verbose = TRUE)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted.

- beginTime:

  Numeric. Start of analysis window in seconds. Default 0 (file start).

- endTime:

  Numeric. End of analysis window in seconds. Default 0 (file end).

- numFormants:

  Integer (1–6). Number of formants to return. The model always predicts
  6 formants internally; this selects the lowest `numFormants` after
  sorting by mean frequency. Default 3.

- windowShift:

  Numeric. Frame shift in milliseconds; sets output frame rate
  (`1000 / windowShift` Hz). Default 5.0 ms (200 Hz). Must be strictly
  less than 32 ms (the 512-sample analysis window at 16 kHz). Values
  other than the training default (5 ms) may slightly reduce accuracy.

- toFile:

  Logical. If `TRUE`, write SSFF output files and return the count
  written (invisibly). If `FALSE`, return an `AsspDataObj`. Default
  `TRUE`.

- explicitExt:

  Character. Output file extension. Default `"fnf"`.

- outputDirectory:

  Character. Directory for output files. `NULL` (default) writes
  alongside the input file.

- verbose:

  Logical. Print per-file progress. Default `TRUE`.

## Value

If `toFile = FALSE`: an `AsspDataObj` with tracks:

- `fm`:

  REAL32, Hz, *n\\frames* × `numFormants`. Formant frequencies; column 1
  = F1, column 2 = F2, etc.

- `bw`:

  REAL32, Hz, *n\\frames* × `numFormants`. Formant bandwidths
  corresponding to each frequency column.

Frame rate: `1000 / windowShift` Hz (default 200 Hz, 5 ms hop). If
`toFile = TRUE`: integer count of files written, returned invisibly.

## Details

ONNX Runtime is installed automatically on first use (~30 MB, cached in
the R user directory). The model file (6.4 MB) is downloaded from the
[FormantNet Hugging Face Hub
repo](https://huggingface.co/FredrikKarlssonSpeech/FormantNet) on first
use (requires the huggingfaceR package and a network connection) and
cached in the R user directory; subsequent calls read the cached copy
with no network access.

Pre-processing (fixed by model training): resample to 16 kHz →
pre-emphasis (0.98) → 512-sample Hann-windowed STFT → 6-pass binomial
spectral envelope smoothing → dB conversion → 257-bin truncation →
global mean/SD normalisation. Post-processing: rescale sigmoid output to
Hz/dB → sort formants by mean frequency → 10-pass time smoothing.

## References

Sakamoto Y, Bunnell HT, Hosner PA (2021). “Neural Formant Tracking.” In
*Proceedings of the 9th International Conference on Speech Prosody (PaPE
2021)*. FormantNet: LSTM-based formant tracker using log-scale spectral
envelopes.

## Examples

``` r
if (FALSE) { # \dontrun{
trk_formant_formantnet(
  system.file("samples", "sustained", "a1.wav", package = "superassp"),
  toFile = FALSE
)
} # }
```
