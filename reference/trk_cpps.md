# Cepstral Peak Prominence Smoothed (CPPS)

Extracts time-series Cepstral Peak Prominence Smoothed (CPPS) via
Praat's PowerCepstrogram. CPPS quantifies voice periodicity and is a
robust correlate of breathiness and dysphonia. Prefer this over
instantaneous CPP when temporal smoothing is desired.

## Usage

``` r
trk_cpps(listOfFiles, ...)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- minF:

  Numeric. Lower quefrency bound for cepstral peak search, in Hz (as
  reciprocal of quefrency). Sets the minimum F0 detectable. Default 60
  Hz.

- maxF:

  Numeric. Upper quefrency bound for cepstral peak search, in Hz. Sets
  the maximum F0 detectable. Default 333 Hz.

- timeStep:

  Numeric. Frame shift for the PowerCepstrogram in seconds. Sets output
  frame rate (1 / timeStep Hz). Default 0.002 s (500 Hz).

- maximumFrequency:

  Numeric. Highest frequency included in the cepstrum in Hz. Default
  5000 Hz. Set to 0 for Nyquist.

- preEmphFrom:

  Numeric. Pre-emphasis onset frequency in Hz. Default 50 Hz.

- windowShape:

  Character. Window shape applied to each analysis frame. Default
  `"Hanning"`.

- relativeWidth:

  Numeric. Relative width of the analysis window. Default 1.0.

- subtractTilt:

  Logical. If `TRUE`, subtract the fitted spectral tilt trend before
  measuring peak prominence (gives CPPS rather than CPP). Default
  `TRUE`.

- timeAveragingWindow:

  Numeric. Duration of the smoothing window along the time axis in
  seconds. Default 0.02 s.

- quefrencyAveragingWindow:

  Numeric. Width of the smoothing window along the quefrency axis in
  seconds. Default 0.0005 s.

- interpolation:

  Character. Peak interpolation method: one of `"none"`, `"parabolic"`,
  `"cubic"`, `"sinc70"`, `"sinc700"`. Default `"parabolic"`.

- trendLineQuefrencyMin:

  Numeric. Minimum quefrency (s) for trend line fitting. Default 0.001
  s.

- trendLineQuefrencyMax:

  Numeric. Maximum quefrency (s) for trend line fitting. Default 0.05 s.

- trendType:

  Character. Shape of the fitted trend: `"straight"` or
  `"exponential decay"`. Default `"exponential decay"`.

- fitMethod:

  Character. Regression method: `"robust"`, `"least squares"`, or
  `"robust slow"`. Default `"robust"`.

- toFile:

  Logical. If `TRUE`, write SSFF output files and return the paths
  written (invisibly). If `FALSE`, return an `AsspDataObj`. Default
  `TRUE`.

- explicitExt:

  Character. Output file extension. Default `"cps"`.

## Value

If `toFile = FALSE`: an `AsspDataObj` with track:

- `cpp`:

  REAL32, dB, n_frames x 1. Cepstral Peak Prominence Smoothed. Higher
  values indicate more periodic (healthier) phonation.

Frame rate: `1 / timeStep` Hz (default 500 Hz). If `toFile = TRUE`:
character vector of output file paths, returned invisibly.

## Details

CPPS is computed via Praat's PowerCepstrogram. Each frame's cepstral
peak prominence is measured relative to a fitted trend line (removing
spectral tilt), then smoothed over `timeAveragingWindow` and
`quefrencyAveragingWindow`. Typical values: 15–25 dB for normal voice;
below 10 dB for breathy or dysphonic voice.

## References

(Hillenbrand et al. 1994)

(Hillenbrand and Houde 1996)

(Heman-Ackah 2003)

## Examples

``` r
if (FALSE) { # \dontrun{
# Extract CPPS from audio file
result <- trk_cpps("speech.wav", toFile = FALSE)

# Plot CPPS over time
plot(result$cpp, type = "l", main = "CPPS", ylab = "CPP (dB)", xlab = "Frame")

# Custom pitch range for female speaker
result <- trk_cpps("speech.wav", minF = 100, maxF = 400, toFile = FALSE)

# Batch process multiple files
trk_cpps(c("f1.wav", "f2.wav", "f3.wav"), toFile = TRUE)
} # }
```
