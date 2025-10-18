#!/usr/bin/env Rscript
# Benchmark for Snack-compatible pitch and formant functions
#
# This script benchmarks trk_snackp() and trk_snackf() to establish
# their performance characteristics and compare with other methods.

library(superassp)
library(bench)
library(dplyr)

cat("=== Snack Functions Benchmark ===\n\n")

# Get test file
test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")

if (!file.exists(test_wav)) {
  stop("Test file not found")
}

# Get file info
info <- av::av_media_info(test_wav)
duration <- info$duration
cat(sprintf("Test file: %s\n", basename(test_wav)))
cat(sprintf("Duration: %.2f seconds\n", duration))
cat(sprintf("Sample rate: %d Hz\n\n", info$audio$sample_rate))

# =============================================================================
# Benchmark snack_pitch
# =============================================================================

cat("Benchmarking trk_snackp()...\n")

pitch_bench <- bench::mark(
  snack_pitch = trk_snackp(test_wav, toFile = FALSE, verbose = FALSE),
  min_iterations = 10,
  max_iterations = 50,
  check = FALSE
)

cat("\nsnack_pitch Results:\n")
summary_cols <- intersect(names(pitch_bench), c("expression", "min", "median", "max", "mean", "mem_alloc"))
print(pitch_bench[, summary_cols])

# =============================================================================
# Benchmark snack_formant
# =============================================================================

cat("\n\nBenchmarking trk_snackf()...\n")

formant_bench <- bench::mark(
  snack_formant_3f = trk_snackf(test_wav, numFormants = 3, toFile = FALSE, verbose = FALSE),
  snack_formant_4f = trk_snackf(test_wav, numFormants = 4, toFile = FALSE, verbose = FALSE),
  snack_formant_5f = trk_snackf(test_wav, numFormants = 5, toFile = FALSE, verbose = FALSE),
  min_iterations = 10,
  max_iterations = 30,
  check = FALSE
)

cat("\nsnack_formant Results:\n")
summary_cols <- intersect(names(formant_bench), c("expression", "min", "median", "max", "mean", "mem_alloc"))
print(formant_bench[, summary_cols])

# =============================================================================
# Comparison with other methods (if available)
# =============================================================================

cat("\n\n=== Comparison with Other Methods ===\n\n")

# Compare pitch tracking methods
cat("Pitch Tracking Comparison:\n")
pitch_comparison <- bench::mark(
  snack_pitch = trk_snackp(test_wav, toFile = FALSE, verbose = FALSE),
  rapt = trk_rapt(test_wav, toFile = FALSE, verbose = FALSE),
  swipe = trk_swipe(test_wav, toFile = FALSE, verbose = FALSE),
  fo_ksvfo = fo(test_wav, toFile = FALSE, verbose = FALSE),
  min_iterations = 8,
  max_iterations = 20,
  check = FALSE
)

cat("\n")
summary_cols <- intersect(names(pitch_comparison), c("expression", "min", "median", "mean", "mem_alloc"))
print(pitch_comparison[, summary_cols])

# Calculate relative speed
snack_median_idx <- which(pitch_comparison$expression == "trk_snackp")
if (length(snack_median_idx) > 0) {
  snack_median <- as.numeric(pitch_comparison$median[snack_median_idx])
  pitch_comparison$relative_speed <- snack_median / as.numeric(pitch_comparison$median)
  
  cat("\n\nRelative to snack_pitch:\n")
  print(pitch_comparison[, c("expression", "relative_speed")])
}

# Compare formant tracking methods
cat("\n\nFormant Tracking Comparison:\n")
formant_comparison <- bench::mark(
  snack_formant = trk_snackf(test_wav, numFormants = 4, toFile = FALSE, verbose = FALSE),
  forest = trk_forest(test_wav, toFile = FALSE, verbose = FALSE),
  min_iterations = 10,
  max_iterations = 30,
  check = FALSE
)

cat("\n")
summary_cols <- intersect(names(formant_comparison), c("expression", "min", "median", "mean", "mem_alloc"))
print(formant_comparison[, summary_cols])

# Calculate relative speed
snack_f_median_idx <- which(formant_comparison$expression == "trk_snackf")
if (length(snack_f_median_idx) > 0) {
  snack_f_median <- as.numeric(formant_comparison$median[snack_f_median_idx])
  formant_comparison$relative_speed <- snack_f_median / as.numeric(formant_comparison$median)
  
  cat("\n\nRelative to snack_formant:\n")
  print(formant_comparison[, c("expression", "relative_speed")])
}

# =============================================================================
# Summary
# =============================================================================

cat("\n\n=== Summary ===\n\n")

cat(sprintf("snack_pitch median time: %.1f ms\n", 
            as.numeric(pitch_bench$median) * 1000))

# Get the correct formant time
formant_4f_idx <- which(formant_bench$expression == "snack_formant_4f")
if (length(formant_4f_idx) > 0) {
  cat(sprintf("snack_formant median time: %.1f ms\n", 
              as.numeric(formant_bench$median[formant_4f_idx]) * 1000))
}

cat("\nSpeed categories:\n")
cat("  - Fastest (< 100 ms): C/C++ implementations (RAPT, SWIPE, fo, forest)\n")
cat("  - Reference (~1500-2500 ms): Snack, Kaldi - Compatibility layers\n")

cat("\n\nUse cases for Snack functions:\n")
cat("  ✓ Replication studies citing Snack\n")
cat("  ✓ Method comparison with Snack-based measurements\n")
cat("  ✓ Historical data processing\n")
cat("  ✓ Benchmark/baseline comparisons\n")

cat("\nRecommendations:\n")
cat("  - General use: trk_rapt()/trk_swipe() for pitch, trk_forest() for formants\n")
cat("  - Snack compatibility: trk_snackp()/trk_snackf()\n")
cat("  - Maximum speed: fo() for pitch\n")

cat("\n=== Benchmark Complete ===\n")
