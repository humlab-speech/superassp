#!/usr/bin/env Rscript
#
# Demonstration of ISO 226:2023 Phon (Loudness Level) Conversions
#
# This script demonstrates the new ISO 226:2023 phon bivariate conversion functions
# for converting between sound pressure level (dB), frequency (Hz), and phon (loudness level).

# Load the package (use devtools::load_all() for development)
if (requireNamespace("devtools", quietly = TRUE) &&
    file.exists("DESCRIPTION")) {
  devtools::load_all(".", quiet = TRUE)
} else {
  library(superassp)
}

cat("=== ISO 226:2023 Phon (Loudness Level) Conversions ===\n\n")

# ---- Introduction ----
cat("## What is Phon?\n\n")
cat("Phon is a unit of loudness level based on equal-loudness contours.\n")
cat("Definition: N phon = sounds as loud as N dB at 1000 Hz reference frequency.\n")
cat("Key insight: Equal SPL at different frequencies doesn't mean equal loudness!\n\n")

# ---- Example 1: Reference Frequency (1000 Hz) ----
cat("## Example 1: Reference Frequency (1000 Hz)\n\n")
cat("At 1000 Hz, phon equals dB by definition:\n")

test_dbs <- c(20, 40, 60, 80)
for (db in test_dbs) {
  phon <- db_and_hz_to_phon(db, 1000)
  cat(sprintf("  %3d dB at 1000 Hz = %5.1f phon\n", db, phon))
}
cat("\n")

# ---- Example 2: Low Frequency Example (100 Hz) ----
cat("## Example 2: Low Frequency Requires More SPL for Same Loudness\n\n")
cat("Human hearing is less sensitive at low frequencies.\n")
cat("To sound as loud as 60 dB at 1000 Hz, we need:\n\n")

target_phon <- 60
freq_100hz <- phon_and_hz_to_db(target_phon, 100)
freq_1000hz <- phon_and_hz_to_db(target_phon, 1000)

cat(sprintf("  60 phon at  100 Hz requires: %5.1f dB SPL\n", freq_100hz))
cat(sprintf("  60 phon at 1000 Hz requires: %5.1f dB SPL\n", freq_1000hz))
cat(sprintf("  Difference: %5.1f dB more at low frequency!\n\n", freq_100hz - freq_1000hz))

# ---- Example 3: High Frequency Sensitivity (4000 Hz) ----
cat("## Example 3: Peak Sensitivity Around 3-4 kHz\n\n")
cat("Human hearing is MOST sensitive around 3-4 kHz.\n")
cat("To sound as loud as 40 dB at 1000 Hz:\n\n")

target_phon <- 40
freq_1000hz <- phon_and_hz_to_db(target_phon, 1000)
freq_4000hz <- phon_and_hz_to_db(target_phon, 4000)

cat(sprintf("  40 phon at 1000 Hz requires: %5.1f dB SPL\n", freq_1000hz))
cat(sprintf("  40 phon at 4000 Hz requires: %5.1f dB SPL\n", freq_4000hz))
cat(sprintf("  Difference: %5.1f dB LESS at 4 kHz (more sensitive!)\n\n", freq_4000hz - freq_1000hz))

# ---- Example 4: Equal-Loudness Contour ----
cat("## Example 4: Equal-Loudness Contour (40 phon)\n\n")
cat("This shows the SPL required at different frequencies for constant loudness:\n\n")

target_phon <- 40
freqs <- c(20, 50, 100, 200, 500, 1000, 2000, 4000, 8000, 12500)

cat("Frequency (Hz) | SPL (dB) for 40 phon\n")
cat("---------------|---------------------\n")
for (freq in freqs) {
  spl <- phon_and_hz_to_db(target_phon, freq)
  cat(sprintf("%14.0f | %8.1f dB\n", freq, spl))
}
cat("\n")

# ---- Example 5: Vectorization ----
cat("## Example 5: Vectorized Computation\n\n")
cat("Both functions support vectorized input for efficient computation:\n\n")

# Generate frequency range
freqs <- seq(100, 8000, by = 100)
phon_level <- 60

# Compute SPL for all frequencies at once
spls <- phon_and_hz_to_db(phon_level, freqs)

cat(sprintf("Computed %d SPL values for %d phon in one call:\n", length(spls), phon_level))
cat(sprintf("  Frequency range: %.0f - %.0f Hz\n", min(freqs), max(freqs)))
cat(sprintf("  SPL range: %.1f - %.1f dB\n", min(spls), max(spls)))
cat(sprintf("  Minimum SPL at: %.0f Hz (%.1f dB) - peak sensitivity\n",
            freqs[which.min(spls)], min(spls)))
cat(sprintf("  Maximum SPL at: %.0f Hz (%.1f dB) - lowest sensitivity\n\n",
            freqs[which.max(spls)], max(spls)))

# ---- Example 6: Round-Trip Conversion ----
cat("## Example 6: Round-Trip Conversion Accuracy\n\n")
cat("Testing conversion accuracy: dB → phon → dB\n\n")

test_cases <- data.frame(
  freq = c(100, 500, 1000, 2000, 4000),
  spl = c(70, 60, 50, 55, 45)
)

cat("Freq (Hz) | Original SPL | Phon | Recovered SPL | Error (dB)\n")
cat("----------|--------------|------|---------------|------------\n")

for (i in 1:nrow(test_cases)) {
  freq <- test_cases$freq[i]
  spl_orig <- test_cases$spl[i]

  phon <- db_and_hz_to_phon(spl_orig, freq)
  spl_recovered <- phon_and_hz_to_db(phon, freq)
  error <- abs(spl_recovered - spl_orig)

  cat(sprintf("%9.0f | %12.1f | %4.1f | %13.1f | %10.4f\n",
              freq, spl_orig, phon, spl_recovered, error))
}
cat("\n")

# ---- Example 7: Hearing Threshold Region ----
cat("## Example 7: Near Hearing Threshold (< 20 phon)\n\n")
cat("ISO 226 is normative for 20-90 phon. Below 20 phon is informative only.\n\n")

suppressWarnings({
  # At threshold, different frequencies require very different SPLs
  threshold_phon <- 10
  freqs_threshold <- c(100, 1000, 4000)

  cat(sprintf("At %d phon (near hearing threshold):\n", threshold_phon))
  for (freq in freqs_threshold) {
    spl <- phon_and_hz_to_db(threshold_phon, freq)
    cat(sprintf("  %4.0f Hz requires: %5.1f dB SPL\n", freq, spl))
  }
})
cat("\n")

# ---- Example 8: Practical Application ----
cat("## Example 8: Practical Application - Sound System Equalization\n\n")
cat("When designing a sound system, equal SPL ≠ equal loudness!\n")
cat("To achieve perceptually flat response at 70 phon:\n\n")

target_phon <- 70
eq_freqs <- c(63, 125, 250, 500, 1000, 2000, 4000, 8000)

cat("Frequency Band (Hz) | Required SPL (dB) | Adjustment from 1kHz\n")
cat("--------------------|-------------------|---------------------\n")

ref_spl <- phon_and_hz_to_db(target_phon, 1000)

for (freq in eq_freqs) {
  spl <- phon_and_hz_to_db(target_phon, freq)
  adjustment <- spl - ref_spl
  cat(sprintf("%19.0f | %17.1f | %+8.1f dB\n", freq, spl, adjustment))
}
cat("\n")

# ---- Example 9: Comparison of Multiple Loudness Levels ----
cat("## Example 9: Equal-Loudness Contours at Different Levels\n\n")
cat("The shape of equal-loudness contours changes with level!\n\n")

phon_levels <- c(20, 40, 60, 80)
freq <- 100

cat("At 100 Hz:\n")
for (phon in phon_levels) {
  spl <- phon_and_hz_to_db(phon, freq)
  cat(sprintf("  %2d phon requires: %5.1f dB SPL\n", phon, spl))
}
cat("\n")

# ---- Example 10: Interpolation for Non-Standard Frequencies ----
cat("## Example 10: Interpolation for Non-Standard Frequencies\n\n")
cat("ISO 226 provides 29 standard 1/3-octave frequencies.\n")
cat("Other frequencies are interpolated (log-frequency linear interpolation):\n\n")

# Standard frequency
spl_standard <- phon_and_hz_to_db(50, 100)
cat(sprintf("  100 Hz (standard):      %5.1f dB for 50 phon\n", spl_standard))

# Interpolated frequency
spl_interp <- phon_and_hz_to_db(50, 110)
cat(sprintf("  110 Hz (interpolated):  %5.1f dB for 50 phon\n", spl_interp))

# Next standard frequency
spl_next <- phon_and_hz_to_db(50, 125)
cat(sprintf("  125 Hz (standard):      %5.1f dB for 50 phon\n\n", spl_next))

cat("Interpolated value falls between the two standard frequencies, as expected.\n\n")

# ---- Summary ----
cat("## Summary\n\n")
cat("Key takeaways:\n")
cat("  - Phon is a bivariate measure: requires both frequency (Hz) AND SPL (dB)\n")
cat("  - At 1000 Hz reference: phon = dB by definition\n")
cat("  - Low frequencies require more SPL for same loudness\n")
cat("  - Peak sensitivity around 3-4 kHz (requires LESS SPL)\n")
cat("  - Equal-loudness contours are non-linear and level-dependent\n")
cat("  - Essential for audio engineering, psychoacoustics, and hearing research\n\n")

cat("Functions:\n")
cat("  - db_and_hz_to_phon(spl_db, freq_hz): Convert (SPL, frequency) → phon\n")
cat("  - phon_and_hz_to_db(phon, freq_hz):   Convert (phon, frequency) → SPL\n\n")

cat("Valid ranges:\n")
cat("  - Frequency: 20-12500 Hz\n")
cat("  - Phon: 20-90 phon (20-4000 Hz), 20-80 phon (5000-12500 Hz)\n")
cat("  - Below 20 phon: informative only (near hearing threshold)\n\n")

cat("Reference: ISO 226:2023 - Acoustics - Normal equal-loudness-level contours\n\n")

cat("=== Demonstration Complete ===\n")
