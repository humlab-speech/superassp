# Pharyngeal Voice Quality Analysis

Extract pharyngealization/voice quality measures from labeled vowel
intervals. Analyzes spectral measures (H1-H2, H1-A1, H1-A2, H1-A3, etc.)
at vowel onset and midpoint with Iseli & Alwan (2004) normalization.

## Usage

``` r
lst_pharyngeal(listOfFiles, textgridPath = NULL, intervalTier = 3, intervalNumber = 1, beginTime = NULL, endTime = NULL, minPitchInitial = 50, maxPitchInitial = 800, toFile = FALSE, return_jstf = FALSE, explicitExt = "pha", outputDirectory = NULL, verbose = TRUE)
```

## Arguments

- listOfFiles:

  Character vector of file paths to audio files (WAV, MP3, etc.)

- textgridPath:

  Character vector of TextGrid file paths (same length as listOfFiles)

- intervalTier:

  Integer specifying which TextGrid tier contains the vowel intervals
  (default: 3)

- intervalNumber:

  Integer specifying which interval to analyze per file (default: 1)

- beginTime:

  Numeric vector of start times in seconds (optional, overrides TextGrid
  if provided)

- endTime:

  Numeric vector of end times in seconds (optional, overrides TextGrid
  if provided)

- minPitchInitial:

  Minimum pitch for initial pitch detection in Hz (default: 50)

- maxPitchInitial:

  Maximum pitch for initial pitch detection in Hz (default: 800)

- toFile:

  Logical - write results to JSTF file (default: FALSE)

- explicitExt:

  File extension for output files (default: "pha")

- outputDirectory:

  Directory for output files (default: NULL = same as input)

- verbose:

  Logical - show progress messages (default: TRUE)

- return_jstf:

  Logical. Return JsonTrackObj instead of data.frame? Default FALSE.
  When both toFile and return_jstf are TRUE, the file is written AND the
  object returned.

## Value

Data frame with 68 columns per file:

- file:

  File path

- start_time, mid_time, end_time, duration_ms:

  Timing information

- f0_start, f0_mid:

  F0 at onset and midpoint (Hz)

- f1/f2/f3_start, f1/f2/f3_mid:

  Formant frequencies at onset/mid (Hz)

- bw1/bw2/bw3_start, bw1/bw2/bw3_mid:

  Formant bandwidths (Hz)

- bw1/bw2/bw3_start_norm, bw1/bw2/bw3_mid_norm:

  Normalized bandwidths (Hz)

- intensity_start, intensity_mid:

  Intensity at onset/mid (dB)

- h1/h2_onset, h1/h2_mid:

  Harmonic amplitudes (dB)

- h1_hz/h2_hz_onset, h1_hz/h2_hz_mid:

  Harmonic frequencies (Hz)

- h1/h2_onset_norm, h1/h2_mid_norm:

  Normalized harmonic amplitudes (dB)

- a1/a2/a3_onset, a1/a2/a3_mid:

  Formant peak amplitudes (dB)

- a1_hz/a2_hz/a3_hz_onset, a1_hz/a2_hz/a3_hz_mid:

  Formant peak frequencies (Hz)

- a3_onset_norm, a3_mid_norm:

  Normalized A3 amplitudes (dB)

- h1_minus_h2/h1_minus_a1/h1_minus_a2/h1_minus_a3_onset:

  Raw differences at onset (dB)

- h1_minus_h2/h1_minus_a1/h1_minus_a2/h1_minus_a3_onset_norm:

  Normalized differences at onset (dB)

- a1_minus_a2/a1_minus_a3/a2_minus_a3_onset:

  Additional raw differences at onset (dB)

- a1_minus_a3_onset_norm, a2_minus_a3_onset_norm:

  Additional normalized differences at onset (dB)

- ...:

  Same structure repeated for midpoint (\*\_mid)

## Details

**Algorithm** (based on scriptPharyFullV4.praat):

1.  Two-pass adaptive pitch detection (speaker-specific F0 range)

2.  Find intensity maxima at onset/midpoint

3.  Extract formants F1/F2/F3 with bandwidths

4.  Extract 40ms Kaiser2 window at onset

5.  Create spectrum with pre-emphasis

6.  Find H1, H2 harmonics near F0, 2×F0

7.  Find A1, A2, A3 formant peaks near F1/F2/F3

8.  Apply Iseli & Alwan (2004) normalization

9.  Calculate differences (raw and normalized)

10. If duration \> 120ms, repeat for midpoint

**Key Measures**:

- **H1-H2**: Open quotient indicator (higher = breathier)

- **H1-A1**: Voice quality indicator (higher = breathier)

- **H1-A2**: Spectral tilt measure

- **H1-A3**: High-frequency energy indicator

- **Normalized (\*\_norm)**: Corrected for formant influence (Iseli &
  Alwan 2004)

**Typical Values**:

- H1-H2: -10 to +10 dB

- H1-A1: -5 to +5 dB

- H1-A3: Often negative (A3 \> H1)

**Performance**: ~24ms per vowel (15.7x faster than v4.8.14)

## Dependencies

Requires `pladdrr` package (\>= 4.8.16)

## References

(Iseli and Alwan 2004)

## See also

[`lst_vq`](https://humlab-speech.github.io/superassp/reference/lst_vq.md),
[`trk_praatsauce`](https://humlab-speech.github.io/superassp/reference/trk_praatsauce.md),
[`lst_voice_report`](https://humlab-speech.github.io/superassp/reference/lst_voice_report.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Analyze labeled vowels from TextGrid
results <- lst_pharyngeal(
  listOfFiles = "speech.wav",
  textgridPath = "speech.TextGrid",
  intervalTier = 3,
  intervalNumber = 1
)

# Analyze specific time ranges (no TextGrid)
results <- lst_pharyngeal(
  listOfFiles = c("vowel1.wav", "vowel2.wav"),
  beginTime = c(0.5, 1.0),
  endTime = c(0.7, 1.3)
)

# Write to JSTF files
lst_pharyngeal(
  listOfFiles = "speech.wav",
  textgridPath = "speech.TextGrid",
  toFile = TRUE
)  # Creates speech.pha

# Access key measures
cat(sprintf("H1-H2* (onset): %.2f dB\\n", results$h1_minus_h2_onset_norm[1]))
cat(sprintf("H1-A1* (onset): %.2f dB\\n", results$h1_minus_a1_onset_norm[1]))
} # }
```
