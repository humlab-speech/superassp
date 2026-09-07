# Acoustic Voice Quality Index (AVQI) using pladdrr

Computes the Acoustic Voice Quality Index (AVQI) from continuous speech
and sustained vowel recordings using pladdrr's Praat bindings. Supports
both AVQI v2.03 (Maryn et al. 2010) and v3.01 (Barsties & Maryn 2015).

## Usage

``` r
lst_avqi(
  svDF,
  csDF,
  version = "v2.03",
  min.sv = 1000,
  toFile = FALSE,
  return_jstf = FALSE,
  explicitExt = "avqi",
  outputDirectory = NULL,
  verbose = TRUE
)
```

## Arguments

- svDF:

  Data frame with sustained vowel samples. Must contain columns:
  listOfFiles, start, end (in milliseconds)

- csDF:

  Data frame with continuous speech samples. Must contain columns:
  listOfFiles, start, end (in milliseconds)

- version:

  Character. AVQI version: "v2.03" (default) or "v3.01"

- min.sv:

  Minimum sustained vowel duration in milliseconds (default: 1000)

- toFile:

  Logical. If TRUE, write results to JSTF file. Default FALSE.

- explicitExt:

  Character. File extension for output. Default "avqi".

- outputDirectory:

  Character. Output directory path. Default NULL (use input directory).

- verbose:

  Logical. Print progress messages (default TRUE)

## Value

If `toFile=FALSE` (default), a list with AVQI measurements. If
`toFile=TRUE`, invisibly returns the path to the written JSTF file.

The list contains:

- version:

  AVQI version used

- avqi:

  Acoustic Voice Quality Index (0-10 scale)

- cpps:

  Cepstral Peak Prominence Smoothed (dB)

- hnr:

  Harmonics-to-Noise Ratio (dB)

- shimmer_local:

  Local shimmer (percent)

- shimmer_db:

  Local shimmer in dB

- slope:

  LTAS slope (dB)

- tilt:

  LTAS tilt (dB)

## Details

The AVQI combines 6 acoustic measures: CPPS, HNR, shimmer (local and
dB), LTAS slope, and LTAS tilt into a single voice quality index (0-10
scale).

## References

(Maryn et al. 2010)

(Barsties and Maryn 2015)

## Examples

``` r
if (FALSE) { # \dontrun{
# Define sustained vowel samples
sv <- data.frame(
  listOfFiles = c("sv1.wav", "sv2.wav"),
  start = c(100, 150),  # milliseconds
  end = c(2500, 2800)
)

# Define continuous speech samples
cs <- data.frame(
  listOfFiles = c("cs1.wav", "cs2.wav"),
  start = c(80, 120),
  end = c(3500, 4000)
)

# Compute AVQI (v2.03)
result <- lst_avqi(sv, cs)
print(result$avqi)

# Compute AVQI v3.01
result_v3 <- lst_avqi(sv, cs, version = "v3.01")

# Write to JSTF file
lst_avqi(sv, cs, toFile = TRUE)
track <- read_track("sv1.avqi")
df <- as.data.frame(track)
} # }
```
