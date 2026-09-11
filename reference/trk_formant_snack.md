# Track formants and bandwidths using the Snack/ESPS LPC tracker

Extracts formant frequencies and bandwidths using the Snack Sound
Toolkit dynamic-programming LPC formant tracker (Talkin / AT&T / KTH,
`jkFormant.c`). This tracker applies normalized cross-correlation for
voicing detection and a DP smoother for temporal continuity. It is a
strong general-purpose alternative to Burg-method trackers, particularly
for lower sample-rate or telephone-bandwidth speech.

## Usage

``` r
trk_formant_snack(
  listOfFiles,
  beginTime = 0,
  endTime = 0,
  numFormants = 4,
  lpcOrder = 12,
  windowLength = 0.049,
  windowShift = 10,
  preEmphasis = 0.7,
  dsFreq = 10000,
  nomF1 = -10,
  lpcType = 0,
  windowType = 2,
  toFile = TRUE,
  explicitExt = "snackfmt",
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

- numFormants:

  Integer. Number of formants to track (1–7). Default 4.

- lpcOrder:

  Integer. LPC order; must satisfy `lpcOrder >= numFormants * 2 + 4`.
  Default 12.

- windowLength:

  Numeric. Analysis window duration in seconds. Default 0.049 s.

- windowShift:

  Numeric. Frame shift in milliseconds; sets output frame rate (1000 /
  windowShift Hz). Default 10.0 ms.

- preEmphasis:

  Numeric. Pre-emphasis filter coefficient (0–1). Default 0.7.

- dsFreq:

  Numeric. Target downsample frequency in Hz; limits formant search to
  0–dsFreq/2 Hz. Default 10000 Hz.

- nomF1:

  Numeric. Nominal F1 for DP cost function. `-10` (default) uses
  built-in defaults derived from `dsFreq`.

- lpcType:

  Integer. LPC method: 0 = autocorrelation, 1 = stabilized covariance, 2
  = covariance. Default 0.

- windowType:

  Integer. Window type: 0 = rectangular, 1 = Hamming, 2 = cos^4, 3 =
  Hanning. Default 2.

- toFile:

  Logical. If `TRUE`, write SSFF output files and return the count
  written invisibly. If `FALSE`, return an `AsspDataObj`. Default
  `TRUE`.

- explicitExt:

  Character. Output file extension. Default `"snackfmt"`.

- outputDirectory:

  Character. Directory for output files. `NULL` (default) writes
  alongside the input file.

- verbose:

  Logical. Print per-file progress. Default `TRUE`.

## Value

If `toFile = FALSE`: an `AsspDataObj` with tracks:

- `fm`:

  REAL32, formant frequencies in Hz, n_frames × numFormants. Columns
  correspond to F1, F2, …, F{numFormants}.

- `bw`:

  REAL32, formant bandwidths in Hz, n_frames × numFormants.

Frame rate: `1000 / windowShift` Hz (default 100 Hz). If
`toFile = TRUE`: integer count of files written, returned invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
# 4-formant tracking
res <- trk_formant_snack("recording.wav", toFile = FALSE)
names(res)  # "fm" "bw"

# Custom LPC order and formant count
trk_formant_snack("speech.mp3", numFormants = 5, lpcOrder = 14)
} # }
```
