# Track formants using Time-Varying Weighted Linear Prediction (TVWLP)

Estimates formant frequencies and bandwidths using GCI-anchored
quasi-closed phase (QCP) weighted, time-varying linear prediction at 8
kHz. TVWLP reduces glottal source contamination compared to standard
LPC, yielding more stable F1 estimates on voiced speech. Prefer this
over Burg-method trackers when F1 tracking in the closed phase is the
priority.

## Usage

``` r
trk_formant_tvwlp(listOfFiles, beginTime = 0, endTime = 0, windowShift = 10, npeaks = 3L, p = 8L, q = 3L, preemp = 0.97, lptype = "tvwlp_l2", toFile = TRUE, explicitExt = "tvf", outputDirectory = NULL, verbose = TRUE)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; audio is resampled to 8 kHz internally.

- beginTime:

  Numeric. Start of analysis window in seconds. Default 0 (file start).

- endTime:

  Numeric. End of analysis window in seconds. Default 0 (file end).

- windowShift:

  Numeric. Frame shift in milliseconds; sets output frame rate (1000 /
  windowShift Hz). Default 10.0 ms.

- npeaks:

  Integer. Number of formants to extract from LP roots. Default 3.

- p:

  Integer. LPC order; must exceed 2 × npeaks for reliable root
  extraction. Default 8.

- q:

  Integer. Polynomial degree for the time-varying LP coefficients
  (controls within-segment variation). Default 3.

- preemp:

  Numeric. Pre-emphasis coefficient (0–1). Default 0.97.

- lptype:

  Character. LP solver method: `"tvwlp_l2"` (GCI-weighted L2, default),
  `"tvlp_l2"` (unweighted L2), `"tvwlp_l1"` (weighted L1, requires
  quantreg), `"tvlp_l1"` (unweighted L1).

- toFile:

  Logical. If `TRUE`, write SSFF output files and return the count
  written invisibly. If `FALSE`, return an `AsspDataObj`. Default
  `TRUE`.

- explicitExt:

  Character. Output file extension. Default `"tvf"`.

- outputDirectory:

  Character. Directory for output files. `NULL` (default) writes
  alongside the input file.

- verbose:

  Logical. Print per-file progress. Default `TRUE`.

## Value

If `toFile = FALSE`: an `AsspDataObj` with tracks:

- `fm`:

  REAL32, formant frequencies in Hz, n_frames × npeaks. Columns
  correspond to F1, F2, …, F{npeaks}.

- `bw`:

  REAL32, formant bandwidths in Hz, n_frames × npeaks.

Frame rate: `1000 / windowShift` Hz (default 100 Hz). If
`toFile = TRUE`: integer count of files written, returned invisibly.

## Details

The weighted methods (`tvwlp_l2`, `tvwlp_l1`) internally estimate pitch
via SRH and GCIs via SEDREAMS at 8 kHz. The unweighted variants skip
pitch estimation and are faster but less accurate on voiced speech. A
5-frame median filter is applied to each formant track before output.

## References

(El-Jaroudi and Makhoul 1991)

## Examples

``` r
if (FALSE) { # \dontrun{
# Default 3-formant TVWLP tracking
res <- trk_formant_tvwlp("recording.wav", toFile = FALSE)
dim(res$fm)  # n_frames x 3

# Unweighted LP (no pitch/GCI estimation)
trk_formant_tvwlp("speech.mp3", lptype = "tvlp_l2")
} # }
```
