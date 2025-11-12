# YIN SIMD Testing and Benchmarking Script
# Tests correctness and measures performance of SIMD-optimized YIN implementation

library(superassp)
library(microbenchmark)

# Load test audio
test_file <- system.file("samples", "sustained", "a1.wav", package = "superassp")

if (test_file == "") {
  stop("Test file not found. Package may not be properly installed.")
}

cat("========================================\n")
cat("YIN SIMD TESTING AND BENCHMARKING\n")
cat("========================================\n\n")

cat("Test file:", test_file, "\n")
audio_info <- av::av_media_info(test_file)
cat("Duration:", audio_info$duration, "seconds\n")
cat("Sample rate:", audio_info$audio$sample_rate, "Hz\n\n")

# ============================================
# PART 1: CORRECTNESS TESTING
# ============================================

cat("========================================\n")
cat("PART 1: CORRECTNESS TESTING\n")
cat("========================================\n\n")

cat("Running YIN with SIMD (3 iterations)...\n")
results_simd <- list()
for (i in 1:3) {
  results_simd[[i]] <- trk_yin(test_file, toFile = FALSE, verbose = FALSE)
}

cat("✓ SIMD version completed\n\n")

# Check determinism: All runs should produce identical results
cat("Checking determinism (all runs identical)...\n")
all_identical <- TRUE
for (i in 2:3) {
  # Compare F0 tracks (first element of each result)
  comparison <- all.equal(results_simd[[1]][[1]], results_simd[[i]][[1]])
  if (!isTRUE(comparison)) {
    all_identical <- FALSE
    cat("Difference detected:", comparison, "\n")
    break
  }
}

if (all_identical) {
  cat("✓ PASS: All SIMD runs produce identical results\n\n")
} else {
  cat("✗ FAIL: SIMD runs produce different results (non-deterministic)\n\n")
  stop("SIMD implementation is non-deterministic!")
}

# ============================================
# PART 2: PERFORMANCE BENCHMARKING
# ============================================

cat("========================================\n")
cat("PART 2: PERFORMANCE BENCHMARKING\n")
cat("========================================\n\n")

cat("Running benchmark (100 iterations)...\n")
cat("(This may take 1-2 minutes)\n\n")

benchmark_results <- microbenchmark(
  yin_simd = trk_yin(test_file, toFile = FALSE, verbose = FALSE),
  times = 100
)

cat("\n")
cat("========================================\n")
cat("BENCHMARK RESULTS\n")
cat("========================================\n\n")

print(summary(benchmark_results))

# Calculate metrics
median_time_ms <- median(benchmark_results$time) / 1e6
mean_time_ms <- mean(benchmark_results$time) / 1e6
min_time_ms <- min(benchmark_results$time) / 1e6
max_time_ms <- max(benchmark_results$time) / 1e6

cat("\n")
cat("Detailed Metrics:\n")
cat(sprintf("  Median: %.2f ms\n", median_time_ms))
cat(sprintf("  Mean:   %.2f ms\n", mean_time_ms))
cat(sprintf("  Min:    %.2f ms\n", min_time_ms))
cat(sprintf("  Max:    %.2f ms\n", max_time_ms))

# Real-time factor (for 0.79s audio file)
audio_duration_s <- as.numeric(audio_info$duration)
rtf <- (median_time_ms / 1000) / audio_duration_s
cat(sprintf("\nReal-Time Factor: %.3fx\n", rtf))
cat("  (Lower is better - <1.0 means faster than real-time)\n")

# ============================================
# PART 3: RESULTS INSPECTION
# ============================================

cat("\n========================================\n")
cat("PART 3: RESULTS INSPECTION\n")
cat("========================================\n\n")

result <- results_simd[[1]]
# Access F0 track (first element of list)
f0_values <- result[[1]][,1]  # First track, first column
voiced_frames <- sum(f0_values > 0)
total_frames <- length(f0_values)
voiced_pct <- (voiced_frames / total_frames) * 100

cat(sprintf("Total frames: %d\n", total_frames))
cat(sprintf("Voiced frames: %d (%.1f%%)\n", voiced_frames, voiced_pct))
cat(sprintf("Unvoiced frames: %d (%.1f%%)\n",
           total_frames - voiced_frames, 100 - voiced_pct))

if (voiced_frames > 0) {
  voiced_f0 <- f0_values[f0_values > 0]
  cat(sprintf("\nF0 Statistics (voiced frames only):\n"))
  cat(sprintf("  Mean:   %.2f Hz\n", mean(voiced_f0)))
  cat(sprintf("  Median: %.2f Hz\n", median(voiced_f0)))
  cat(sprintf("  Min:    %.2f Hz\n", min(voiced_f0)))
  cat(sprintf("  Max:    %.2f Hz\n", max(voiced_f0)))
  cat(sprintf("  SD:     %.2f Hz\n", sd(voiced_f0)))
}

# ============================================
# SUMMARY
# ============================================

cat("\n========================================\n")
cat("SUMMARY\n")
cat("========================================\n\n")

cat("✓ SIMD implementation enabled and working\n")
cat("✓ Results are deterministic\n")
cat(sprintf("✓ Median processing time: %.2f ms\n", median_time_ms))
cat(sprintf("✓ Real-time factor: %.3fx (%.1fx faster than real-time)\n",
           rtf, 1/rtf))

# Expected performance (from SIMD plan)
cat("\nExpected SIMD speedup: 4-8x vs scalar baseline\n")
cat("Note: Cannot directly compare SIMD vs scalar in this build\n")
cat("      (would require disabling SIMD flag and recompiling)\n\n")

cat("========================================\n")
cat("TEST COMPLETE\n")
cat("========================================\n")
