# Dysphonia Severity Index (DSI) Analysis (pladdrr)

Compute the Dysphonia Severity Index (DSI) using pladdrr. DSI is a
multiparametric voice quality assessment combining maximum phonation
time, softest intensity, highest fundamental frequency, and jitter.

## Usage

``` r
lst_dsi(
  softDF,
  highpitchDF,
  maxprolongedDF,
  stableDF = NULL,
  use.calibration = FALSE,
  db.calibration = 10,
  speaker.name = NULL,
  speaker.ID = NULL,
  speaker.dob = NULL,
  session.datetime = NULL,
  pdf.path = NULL,
  simple.output = FALSE,
  overwrite.pdfs = FALSE,
  praat_path = NULL,
  toFile = FALSE,
  return_jstf = FALSE,
  explicitExt = "dsi",
  outputDirectory = NULL
)
```

## Arguments

- softDF:

  Data frame with soft phonation samples. Columns: absolute_file_path,
  start, end

- highpitchDF:

  Data frame with high pitch samples. Columns: absolute_file_path,
  start, end

- maxprolongedDF:

  Data frame with maximally prolonged vowel samples. Columns:
  absolute_file_path, start, end

- stableDF:

  Data frame with stable vowel samples for jitter (optional, defaults to
  maxprolongedDF)

- use.calibration:

  Logical. Apply calibration to intensity measurements. Default FALSE

- db.calibration:

  Numeric. Calibration factor in dB. Default 10

- speaker.name:

  Character. Speaker name (optional, for identification)

- speaker.ID:

  Character. Speaker ID (optional, for file naming)

- speaker.dob:

  Character. Speaker date of birth (optional)

- session.datetime:

  Character. Session datetime (optional)

- pdf.path:

  Character. Path for PDF output (not implemented in pladdrr version)

- simple.output:

  Logical. Simplified output (not used in pladdrr version)

- overwrite.pdfs:

  Logical. Overwrite PDFs (not used in pladdrr version)

- praat_path:

  Character. Praat path (not used in pladdrr version)

- toFile:

  Logical. Write to JSTF file? Default FALSE

- explicitExt:

  Character. Output file extension. Default "dsi"

- outputDirectory:

  Character. Output directory. NULL = first file's directory. Default
  NULL

## Value

If toFile=FALSE, list with 6 elements:

- ID:

  Speaker ID

- Maximum_phonation_time:

  MPT in seconds

- Softest_intensity_of_voiced_speech:

  Minimum intensity in dB

- Maximum_fundamental_frequency:

  Maximum F0 in Hz

- Jitter_ppq5:

  5-point period perturbation quotient (%)

- Dysphonia_Severity_Index:

  DSI composite score

If toFile=TRUE, invisibly returns output file path.

## Details

DSI Formula (Wuyts et al., 2000): DSI = 1.127 + 0.164*MPT -
0.038*I-low + 0.0053*F0-high - 5.30*Jitter

## References

(Wuyts et al. 2000)

## Examples

``` r
if (FALSE) { # \dontrun{
# Example DSI analysis
soft <- data.frame(
  absolute_file_path = "soft1.wav",
  start = 0,
  end = 2
)

high <- data.frame(
  absolute_file_path = "high1.wav",
  start = 0,
  end = 1.5
)

prolonged <- data.frame(
  absolute_file_path = c("vowel1.wav", "vowel2.wav"),
  start = c(0, 0),
  end = c(5, 6)
)

result <- lst_dsi(
  softDF = soft,
  highpitchDF = high,
  maxprolongedDF = prolonged,
  speaker.ID = "P001"
)

print(result)
} # }
```
