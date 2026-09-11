# Signal polarity detection (RESKEW algorithm)

Detects microphone/recording polarity from LP residual skewness. The
RESKEW algorithm compares residual skewness characteristics with and
without high-pass filtering to determine signal polarity.

## Usage

``` r
lst_polarity(
  listOfFiles,
  beginTime = 0,
  endTime = 0,
  toFile = FALSE,
  return_jstf = FALSE,
  verbose = TRUE
)
```

## Arguments

- listOfFiles:

  Vector of file paths (WAV, MP3, MP4, etc.) to analyze

- beginTime:

  Start time in seconds (0 for beginning of file)

- endTime:

  End time in seconds (0 for end of file)

- toFile:

  Write output to file (default: FALSE, not supported for scalar output)

- verbose:

  Show progress messages (default: TRUE)

- return_jstf:

  Logical. Return JsonTrackObj instead of data.frame? Default FALSE.
  When both toFile and return_jstf are TRUE, the file is written AND the
  object returned.

## Value

Data frame with columns:

- `file`: Input file basename

- `polarity`: Polarity sign (+1 or -1)

## Details

**Algorithm** (RESKEW, Drugman et al.):

1.  High-pass filter signal at 490 Hz cutoff (removes low-frequency
    noise/drift)

2.  Compute LP residual with order = fs/1000 + 2 samples, 25ms frames,
    5ms shift

3.  Compare residual skewness: unfiltered vs. high-pass filtered

4.  Polarity = sign(skew_filtered - skew_unfiltered)

**Interpretation**:

- `+1`: Normal polarity (positive peaks are vocal pulses)

- `-1`: Inverted polarity (flip signal before further processing)

**Note**: If skewness values are very close, result may be unstable.
Recommend confidence threshold: \|result\| \> 0.1 for reliable polarity
detection.

## Examples

``` r
if (FALSE) { # \dontrun{
# Single file
pol <- lst_polarity("speech.wav")
cat("Polarity:", pol$polarity, "\n")

# Batch process
files <- c("file1.wav", "file2.wav")
polarities <- lst_polarity(files)
} # }
```
