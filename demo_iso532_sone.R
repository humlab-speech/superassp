#!/usr/bin/env Rscript
#
# Demonstration of ISO 532 Sone (Loudness) Conversions
#
# This script demonstrates the ISO 532 sone conversions for measuring
# perceived loudness on a linear scale, building on ISO 226:2023 phon conversions.

# Load the package (use devtools::load_all() for development)
if (requireNamespace("devtools", quietly = TRUE) &&
    file.exists("DESCRIPTION")) {
  devtools::load_all(".", quiet = TRUE)
} else {
  library(superassp)
}

cat("=== ISO 532 Sone (Loudness) Conversions ===\n\n")

# ---- Introduction ----
cat("## What is Sone?\n\n")
cat("Sone is a linear unit of perceived loudness proposed by Stanley Smith Stevens (1936).\n")
cat("Definition: 1 sone = 40 phons = 40 dB SPL at 1 kHz reference frequency.\n")
cat("Key insight: Sones provide a LINEAR scale where 2 sones is TWICE as loud as 1 sone!\n\n")

cat("Relationship to Phon:\n")
cat("  - Phon is LOGARITHMIC (like dB)\n")
cat("  - Sone is LINEAR (matches perceived loudness directly)\n")
cat("  - Each 10 phon increase ≈ doubles loudness in sones\n\n")

# ---- Example 1: Reference Value ----
cat("## Example 1: Reference Value (1 sone = 40 phon)\n\n")

sone <- phon_to_sone(40)
phon <- sone_to_phon(1)

cat(sprintf("  40 phon = %.2f sone (by definition)\n", sone))
cat(sprintf("  1 sone = %.1f phon (by definition)\n\n", phon))

# ---- Example 2: Doubling Property ----
cat("## Example 2: Doubling Property\n\n")
cat("The key advantage of sones: doubling the value = doubling perceived loudness!\n\n")

phons <- c(40, 50, 60, 70)
sones <- phon_to_sone(phons)

cat("Phon to Sone:\n")
for (i in seq_along(phons)) {
  ratio <- if (i > 1) sones[i] / sones[i-1] else NA
  if (is.na(ratio)) {
    cat(sprintf("  %2d phon = %5.2f sone\n", phons[i], sones[i]))
  } else {
    cat(sprintf("  %2d phon = %5.2f sone (×%.2f from previous)\n",
                phons[i], sones[i], ratio))
  }
}
cat("\n")

cat("Sone to Phon:\n")
sone_values <- c(1, 2, 4, 8)
phon_values <- sone_to_phon(sone_values)

for (i in seq_along(sone_values)) {
  diff <- if (i > 1) phon_values[i] - phon_values[i-1] else NA
  if (is.na(diff)) {
    cat(sprintf("  %.0f sone = %5.1f phon\n", sone_values[i], phon_values[i]))
  } else {
    cat(sprintf("  %.0f sone = %5.1f phon (+%.1f phon from previous)\n",
                sone_values[i], phon_values[i], diff))
  }
}
cat("\n")

# ---- Example 3: Linear Loudness Perception ----
cat("## Example 3: Linear Loudness Perception\n\n")
cat("Unlike phon (logarithmic), sones match human perception of 'twice as loud':\n\n")

cat("Loudness Level | Perceived Loudness\n")
cat("---------------|-------------------\n")
loudness_table <- data.frame(
  sone = c(0.5, 1, 2, 4, 8),
  description = c("Half as loud", "Reference", "Twice as loud",
                  "Four times as loud", "Eight times as loud")
)

for (i in 1:nrow(loudness_table)) {
  phon_val <- sone_to_phon(loudness_table$sone[i])
  cat(sprintf("%6.1f sone (%4.1f phon) | %s\n",
              loudness_table$sone[i], phon_val, loudness_table$description[i]))
}
cat("\n")

# ---- Example 4: Zwicker vs Moore-Glasberg Methods ----
cat("## Example 4: Zwicker vs Moore-Glasberg Methods\n\n")
cat("ISO 532 defines two methods for phon-sone conversion:\n")
cat("  - Zwicker (ISO 532-1): Piecewise mathematical formulas\n")
cat("  - Moore-Glasberg (ISO 532-2): Lookup table with interpolation\n\n")

test_phons <- c(20, 40, 60, 80, 100)

cat("Phon | Zwicker Sone | Moore-Glasberg Sone | Difference\n")
cat("-----|--------------|---------------------|------------\n")

for (phon_val in test_phons) {
  sone_z <- phon_to_sone(phon_val, method = "zwicker")
  sone_mg <- phon_to_sone(phon_val, method = "moore-glasberg")
  diff_pct <- abs(sone_z - sone_mg) / sone_z * 100

  cat(sprintf("%4.0f | %12.2f | %19.2f | %6.1f%%\n",
              phon_val, sone_z, sone_mg, diff_pct))
}
cat("\n")

# ---- Example 5: Round-Trip Conversion Accuracy ----
cat("## Example 5: Round-Trip Conversion Accuracy\n\n")
cat("Testing: phon → sone → phon (Zwicker method)\n\n")

cat("Original Phon | Sone | Recovered Phon | Error\n")
cat("--------------|------|----------------|-------\n")

test_phons_rt <- c(10, 20, 30, 40, 50, 60, 70, 80, 90)

for (phon_orig in test_phons_rt) {
  sone <- phon_to_sone(phon_orig)
  phon_recovered <- sone_to_phon(sone)
  error <- abs(phon_recovered - phon_orig)

  cat(sprintf("%13.0f | %4.2f | %14.1f | %5.3f\n",
              phon_orig, sone, phon_recovered, error))
}
cat("\n")

# ---- Example 6: Stevens' Power Law ----
cat("## Example 6: Stevens' Power Law\n\n")
cat("Loudness follows a power law: Loudness ∝ Intensity^0.3\n")
cat("This means 10× intensity (10 dB) gives 10^0.3 ≈ 2× loudness\n\n")

intensities_db <- c(0, 10, 20, 30, 40)
phons_ref <- 40 + intensities_db  # At 1 kHz, phon = dB
sones_ref <- phon_to_sone(phons_ref)

cat("Intensity Increase | dB at 1kHz | Phon | Sone | Loudness Ratio\n")
cat("-------------------|------------|------|------|----------------\n")

for (i in seq_along(intensities_db)) {
  if (i == 1) {
    cat(sprintf("%18s | %10.0f | %4.0f | %4.2f | (reference)\n",
                "1×", phons_ref[i], phons_ref[i], sones_ref[i]))
  } else {
    intensity_ratio <- 10^(intensities_db[i] / 10)
    loudness_ratio <- sones_ref[i] / sones_ref[1]
    cat(sprintf("%17.0f× | %10.0f | %4.0f | %4.2f | %.2f×\n",
                intensity_ratio, phons_ref[i], phons_ref[i],
                sones_ref[i], loudness_ratio))
  }
}
cat("\n")

# ---- Example 7: Combined dB/Hz to Sone Conversion ----
cat("## Example 7: Direct dB + Hz → Sone Conversion\n\n")
cat("Convenience functions combine ISO 226 (phon) and ISO 532 (sone) conversions:\n\n")

# At 1 kHz reference
sone_1khz <- db_and_hz_to_sone(40, 1000)
cat(sprintf("40 dB at 1000 Hz = %.2f sone\n", sone_1khz))

sone_1khz_50 <- db_and_hz_to_sone(50, 1000)
cat(sprintf("50 dB at 1000 Hz = %.2f sone (2× louder)\n\n", sone_1khz_50))

# Frequency affects perceived loudness
cat("Same SPL at different frequencies gives different loudness:\n\n")

freqs <- c(100, 500, 1000, 2000, 4000, 8000)
spl <- 60

cat("Frequency (Hz) | 60 dB SPL → Sone | Relative Loudness\n")
cat("---------------|------------------|------------------\n")

sones_at_60db <- sapply(freqs, function(f) db_and_hz_to_sone(spl, f))
ref_sone <- db_and_hz_to_sone(spl, 1000)

for (i in seq_along(freqs)) {
  rel_loudness <- sones_at_60db[i] / ref_sone
  cat(sprintf("%14.0f | %16.3f | %.2f× (vs 1 kHz)\n",
              freqs[i], sones_at_60db[i], rel_loudness))
}
cat("\n")

# ---- Example 8: Inverse Conversion (Sone + Hz → dB) ----
cat("## Example 8: Sone + Hz → dB Conversion\n\n")
cat("To achieve constant perceived loudness across frequencies:\n\n")

target_sone <- 2.0  # Target: 2 sones (twice reference loudness)

cat(sprintf("Target loudness: %.1f sone\n", target_sone))
cat("Required SPL at different frequencies:\n\n")

cat("Frequency (Hz) | Required SPL (dB) | Adjustment from 1kHz\n")
cat("---------------|-------------------|---------------------\n")

spl_at_1khz <- sone_and_hz_to_db(target_sone, 1000)

for (freq in freqs) {
  spl_required <- sone_and_hz_to_db(target_sone, freq)
  adjustment <- spl_required - spl_at_1khz

  cat(sprintf("%14.0f | %17.1f | %+8.1f dB\n",
              freq, spl_required, adjustment))
}
cat("\n")

# ---- Example 9: Practical Application - Audio Mixing ----
cat("## Example 9: Practical Application - Audio Mixing\n\n")
cat("When mixing audio, adjustments in dB don't directly map to perceived loudness.\n")
cat("Use sones for perceptually-balanced mixing:\n\n")

# Mix scenario: balance vocals and instruments
vocal_db <- 50
instrument_db <- 48

vocal_sone <- db_and_hz_to_sone(vocal_db, 1000)
instrument_sone <- db_and_hz_to_sone(instrument_db, 500)

cat(sprintf("Vocals:      %2d dB at 1000 Hz = %.3f sone\n", vocal_db, vocal_sone))
cat(sprintf("Instrument:  %2d dB at  500 Hz = %.3f sone\n", instrument_db, instrument_sone))
cat(sprintf("Loudness ratio: %.2f:1 (vocals:instrument)\n\n", vocal_sone / instrument_sone))

# Adjust for equal loudness
target_ratio <- 1.0  # Equal perceived loudness
new_instrument_sone <- vocal_sone / target_ratio
new_instrument_db <- sone_and_hz_to_db(new_instrument_sone, 500)

cat(sprintf("To achieve 1:1 perceived loudness balance:\n"))
cat(sprintf("  Adjust instrument from %d dB to %.1f dB (+%.1f dB)\n\n",
            instrument_db, new_instrument_db, new_instrument_db - instrument_db))

# ---- Example 10: Low Loudness Levels (< 40 phon) ----
cat("## Example 10: Low Loudness Levels (< 40 phon)\n\n")
cat("Below 40 phon (1 sone), different formulas apply:\n\n")

low_phons <- c(0, 10, 20, 30, 40)
low_sones <- phon_to_sone(low_phons)

cat("Phon | Sone | Notes\n")
cat("-----|------|-------\n")

for (i in seq_along(low_phons)) {
  notes <- if (low_phons[i] < 20) {
    "Near hearing threshold"
  } else if (low_phons[i] < 40) {
    "Below reference (uses S = (L/40)^(1/0.35))"
  } else {
    "Reference level"
  }

  cat(sprintf("%4.0f | %4.3f | %s\n", low_phons[i], low_sones[i], notes))
}
cat("\n")

# ---- Example 11: Comparison Across Scales ----
cat("## Example 11: Complete Psychoacoustic Scale Comparison\n\n")
cat("Example: 200 Hz tone at 60 dB SPL\n\n")

freq <- 200
spl <- 60

# Convert through scales
phon <- db_and_hz_to_phon(spl, freq)
sone <- phon_to_sone(phon)
bark <- hz_to_bark(freq)
erb <- hz_to_erb(freq)
mel <- hz_to_mel(freq)
semitone <- hz_to_semitone(freq)

cat("Physical:\n")
cat(sprintf("  Frequency: %.0f Hz\n", freq))
cat(sprintf("  SPL: %.0f dB\n\n", spl))

cat("Psychoacoustic (Pitch):\n")
cat(sprintf("  Bark: %.2f\n", bark))
cat(sprintf("  ERB: %.2f\n", erb))
cat(sprintf("  Mel: %.1f\n", mel))
cat(sprintf("  Semitone: %.1f (MIDI note %.0f)\n\n", semitone, semitone))

cat("Psychoacoustic (Loudness):\n")
cat(sprintf("  Phon: %.1f (loudness level)\n", phon))
cat(sprintf("  Sone: %.2f (perceived loudness)\n\n", sone))

cat("Interpretation:\n")
cat(sprintf("  - This tone is perceived as %.2f× the loudness of 40 dB at 1 kHz\n", sone))
cat(sprintf("  - It has the same loudness level as %.1f dB at 1 kHz\n", phon))
cat(sprintf("  - Its pitch is in critical band %.1f (Bark scale)\n\n", bark))

# ---- Summary ----
cat("## Summary\n\n")
cat("Key advantages of sones over phon:\n")
cat("  - LINEAR scale matches perceived loudness (2 sones = 2× louder)\n")
cat("  - Easier for perceptual calculations (mixing, balancing)\n")
cat("  - Direct representation of Stevens' power law\n\n")

cat("Key relationships:\n")
cat("  - 1 sone = 40 phons = 40 dB at 1 kHz (by definition)\n")
cat("  - Each 10 phon ≈ doubles loudness in sones\n")
cat("  - 10× intensity (10 dB) ≈ 2× loudness (Stevens' power law)\n\n")

cat("Functions:\n")
cat("  - phon_to_sone(phon): Convert loudness level to loudness\n")
cat("  - sone_to_phon(sone): Convert loudness to loudness level\n")
cat("  - db_and_hz_to_sone(spl_db, freq_hz): Direct (dB, Hz) → sone\n")
cat("  - sone_and_hz_to_db(sone, freq_hz): Direct (sone, Hz) → dB\n\n")

cat("Methods:\n")
cat("  - Zwicker (ISO 532-1): Piecewise formulas (default)\n")
cat("  - Moore-Glasberg (ISO 532-2): Lookup table with interpolation\n\n")

cat("References:\n")
cat("  - ISO 532-1:2017 - Zwicker method\n")
cat("  - ISO 532-2:2017 - Moore-Glasberg method\n")
cat("  - ISO 226:2023 - Equal-loudness contours (for phon calculation)\n")
cat("  - Stevens (1936) - Original sone scale definition\n\n")

cat("=== Demonstration Complete ===\n")
