# Extract Dysprosody Prosodic Features

Extracts 193 prosodic features using pladdrr-based dysprosody pipeline
including F0 analysis (MOMEL/INTSINT), intensity, formants, and spectral
tilt measures with Iseli-Alwan harmonic correction.

## Usage

``` r
lst_dysprosody(listOfFiles, ...)
```

## Arguments

- listOfFiles:

  Vector of file paths to audio files

- beginTime:

  Start time in seconds (default: 0.0)

- endTime:

  End time in seconds (default: 0.0 = full duration)

- minF:

  Minimum F0 in Hz (default: 60)

- maxF:

  Maximum F0 in Hz (default: 750)

- windowShift:

  Window shift in milliseconds (default: 10)

- toFile:

  Write output to .dyp files (default: FALSE)

- explicitExt:

  Output file extension (default: "dyp")

- outputDirectory:

  Output directory (default: NULL = same as input)

- verbose:

  Show progress (default: TRUE)

## Value

If toFile=FALSE: list of data frames (one per file) with 193 features.
If toFile=TRUE: invisibly returns vector of output file paths.

## Details

**Requires**: pladdrr \>= 4.8.23

**Features Extracted** (193 total):

- F0 analysis: MOMEL quadratic spline modeling, INTSINT tonal coding

- Intensity measures: statistics of intensity contour at INTSINT targets

- Spectral tilt: 7 measures including Iseli-Alwan harmonic correction
  (L2L1, L2cL1c, L1cLF3c, L1LF3, SLF, C1, Spectral Balance)

- Time-series statistics: mean, std, variation, IQR, max, min for all
  measures and their first differences

- Global measures: duration, pitch key, pitch range, INTSINT
  concentration

**Performance**: ~10-12 seconds per file (optimized via batch queries)

**Output Format**: When toFile=TRUE, writes JSON Track Format (.dyp)
files with all 193 features organized by time slice.

## References

(Villarubia and others 2025)

## Examples

``` r
if (FALSE) { # \dontrun{
# Extract dysprosody features from audio
result <- lst_dysprosody("speech.wav", toFile = FALSE)

# Access specific features
cat("Duration:", result$Duration, "seconds\n")
cat("Pitch mean:", result$momelpitchtmean, "Hz\n")
cat("Pitch key:", result$PitchKey, "Hz\n")
cat("Pitch range:", result$PitchRange, "octaves\n")

# Batch processing with file output
files <- c("speech1.wav", "speech2.wav")
lst_dysprosody(files, toFile = TRUE, outputDirectory = "output")
} # }
```
