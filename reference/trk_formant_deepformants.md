# Track formant frequencies using DeepFormants (ONNX)

Predicts formant frequencies from composite LPC/cepstral features using
a stacked LSTM (DeepFormants; (Dissen and Keshet 2016) ). The model was
trained on 16 kHz speech and predicts up to four formants without
bandwidth estimates. All audio is resampled automatically; no Python
required — inference uses the bundled ONNX Runtime.

## Usage

``` r
trk_formant_deepformants(
  listOfFiles,
  beginTime = 0,
  endTime = 0,
  numFormants = 3L,
  windowShift = 10,
  toFile = TRUE,
  explicitExt = "dff",
  outputDirectory = NULL,
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

- numFormants:

  Integer (1–4). Number of formants to return. The model always predicts
  4 formants; this selects the first `numFormants`. Default 3.

- windowShift:

  Numeric. Frame shift in milliseconds; sets output frame rate
  (`1000 / windowShift` Hz). Default 10.0 ms (100 Hz). Frame length is
  fixed at 30 ms; `windowShift` controls overlap. Values other than the
  training default (10 ms) may reduce accuracy.

- toFile:

  Logical. If `TRUE`, write SSFF output files and return the count
  written (invisibly). If `FALSE`, return an `AsspDataObj`. Default
  `TRUE`.

- explicitExt:

  Character. Output file extension. Default `"dff"`.

- outputDirectory:

  Character. Directory for output files. `NULL` (default) writes
  alongside the input file.

- verbose:

  Logical. Print per-file progress. Default `TRUE`.

## Value

If `toFile = FALSE`: an `AsspDataObj` with track:

- `fm`:

  REAL32, Hz, *n_frames* × `numFormants`. Formant frequencies; column 1
  = F1, column 2 = F2, etc. No bandwidth track is produced.

Frame rate: `1000 / windowShift` Hz (default 100 Hz, 10 ms hop). If
`toFile = TRUE`: integer count of files written, returned invisibly.

## Details

ONNX Runtime is installed automatically on first use (~30 MB, cached in
the R user directory). The model file (~10 MB) is downloaded from the
[DeepFormants Hugging Face Hub
repo](https://huggingface.co/FredrikKarlssonSpeech/DeepFormants) on
first use (requires the huggingfaceR package and a network connection)
and cached in the R user directory; subsequent calls read the cached
copy with no network access.

Pre-processing (fixed by model training): resample to 16 kHz (int16
scale, no normalisation) → 480-sample (30 ms) frames → 350-dim feature
vector per frame: 50 specPS coefficients (averaged periodogram → log →
DCT-II) plus 300 arspec coefficients (Levinson-Durbin LPC orders 8–17 →
AR spectrum → same transform, 30 per order). Post-processing: raw output
× 1000 = Hz.

## References

Dissen Y, Keshet J (2016). “DeepFormants: Automatic formant tracking and
estimation using deep networks.” Software package for formant tracking
using LSTM networks with LPC/cepstral features,
<https://github.com/MLSpeech/DeepFormants>.

## Examples

``` r
if (FALSE) { # \dontrun{
trk_formant_deepformants(
  system.file("samples", "sustained", "a1.wav", package = "superassp"),
  toFile = FALSE
)
} # }
```
