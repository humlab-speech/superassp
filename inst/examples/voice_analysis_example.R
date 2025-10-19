#!/usr/bin/env Rscript
#
# Example script demonstrating the Voice Analysis Toolbox (lst_vat) function
# This script shows how to compute 132 dysphonia measures from sustained vowel recordings
#

library(superassp)

# =============================================================================
# Setup
# =============================================================================

cat("Voice Analysis Toolbox - Example Usage\n")
cat("=======================================\n\n")

# Check if voice_analysis module is available
if (!voice_analysis_available()) {
  cat("Installing voice_analysis module...\n")
  cat("This may take a few minutes on first run.\n\n")

  # Install with automatic method selection
  install_voice_analysis()

  cat("\nInstallation complete!\n\n")
}

# Display system information
cat("System Information:\n")
info <- voice_analysis_info()
cat(sprintf("  Platform: %s (%s)\n", info$platform, info$machine))
cat(sprintf("  CPU cores: %d physical, %d logical\n",
            info$cpu_count_physical, info$cpu_count))
cat(sprintf("  Cython extensions: %s\n",
            ifelse(info$cython_available, "available", "not available")))
cat(sprintf("  Numba support: %s\n",
            ifelse(info$numba_available, "available", "not available")))
cat(sprintf("  Recommended workers: %d\n\n", info$recommended_workers))

# =============================================================================
# Example 1: Basic Usage - Single File
# =============================================================================

cat("Example 1: Basic Usage\n")
cat("----------------------\n")

# Get a test file
test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")

if (test_wav == "" || !file.exists(test_wav)) {
  cat("Warning: Test file not found. Skipping examples.\n")
  cat("Please provide your own sustained vowel recording.\n")
  quit(save = "no", status = 0)
}

cat(sprintf("Analyzing: %s\n", basename(test_wav)))

# Analyze the recording
result <- lst_vat(test_wav, verbose = TRUE)

if (result$success) {
  cat("\nAnalysis successful!\n")
  cat(sprintf("Computed %d measures\n", length(result$measures)))
  cat(sprintf("Sample rate: %d Hz\n\n", result$fs))

  # Display some key measures
  cat("Key Dysphonia Measures:\n")
  cat("  Nonlinear Dynamics:\n")
  cat(sprintf("    DFA:  %.6f (Detrended Fluctuation Analysis)\n", result$measures$DFA))
  cat(sprintf("    RPDE: %.6f (Recurrence Period Density Entropy)\n", result$measures$RPDE))
  cat(sprintf("    PPE:  %.6f (Pitch Period Entropy)\n\n", result$measures$PPE))

  # Display some HNR measures
  hnr_measures <- result$measures[grepl("^HNR", names(result$measures))]
  if (length(hnr_measures) > 0) {
    cat("  Harmonic-to-Noise Ratios:\n")
    for (name in names(hnr_measures)) {
      cat(sprintf("    %s: %.2f dB\n", name, hnr_measures[[name]]))
    }
    cat("\n")
  }

  # Display measure categories
  cat("Feature Categories:\n")
  cat(sprintf("  Jitter measures:   %d\n",
              sum(grepl("jitter", names(result$measures), ignore.case = TRUE))))
  cat(sprintf("  Shimmer measures:  %d\n",
              sum(grepl("shimmer", names(result$measures), ignore.case = TRUE))))
  cat(sprintf("  HNR/NHR measures:  %d\n",
              sum(grepl("HNR|NHR", names(result$measures)))))
  cat(sprintf("  GNE measures:      %d\n",
              sum(grepl("GNE", names(result$measures)))))
  cat(sprintf("  MFCC measures:     %d\n",
              sum(grepl("MFCC", names(result$measures), ignore.case = TRUE))))
  cat(sprintf("  Wavelet measures:  %d\n",
              sum(grepl("wavelet", names(result$measures), ignore.case = TRUE))))

} else {
  cat(sprintf("\nAnalysis failed: %s\n", result$error))
}

cat("\n")

# =============================================================================
# Example 2: Custom F0 Range
# =============================================================================

cat("Example 2: Custom F0 Range (Male Voice)\n")
cat("----------------------------------------\n")

# Analyze with male voice F0 range (75-300 Hz)
result_male <- lst_vat(test_wav,
                       f0_min = 75,
                       f0_max = 300,
                       verbose = FALSE)

if (result_male$success) {
  cat(sprintf("DFA (male range): %.6f\n", result_male$measures$DFA))
}

# Analyze with female voice F0 range (100-500 Hz)
result_female <- lst_vat(test_wav,
                        f0_min = 100,
                        f0_max = 500,
                        verbose = FALSE)

if (result_female$success) {
  cat(sprintf("DFA (female range): %.6f\n", result_female$measures$DFA))
}

cat("\n")

# =============================================================================
# Example 3: Time Windowing
# =============================================================================

cat("Example 3: Time Windowing\n")
cat("-------------------------\n")

# Get file duration
file_info <- av::av_media_info(test_wav)
duration <- file_info$duration

cat(sprintf("File duration: %.2f seconds\n", duration))

if (duration > 1.0) {
  # Analyze middle portion
  mid_start <- duration * 0.25
  mid_end <- duration * 0.75

  cat(sprintf("Analyzing %.2f to %.2f seconds\n", mid_start, mid_end))

  result_windowed <- lst_vat(test_wav,
                            beginTime = mid_start,
                            endTime = mid_end,
                            verbose = FALSE)

  if (result_windowed$success) {
    cat(sprintf("DFA (windowed): %.6f\n", result_windowed$measures$DFA))
  }
}

cat("\n")

# =============================================================================
# Example 4: Parallel Processing
# =============================================================================

cat("Example 4: Parallel Processing\n")
cat("-------------------------------\n")

# Use recommended number of workers for best performance
n_workers <- min(info$recommended_workers, 4)  # Cap at 4 for demo

cat(sprintf("Using %d workers for parallel processing\n", n_workers))

# Time the analysis
start_time <- Sys.time()
result_parallel <- lst_vat(test_wav,
                          n_cores = n_workers,
                          verbose = FALSE)
end_time <- Sys.time()

elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))

if (result_parallel$success) {
  cat(sprintf("Analysis completed in %.2f seconds\n", elapsed))
  cat(sprintf("DFA (parallel): %.6f\n", result_parallel$measures$DFA))
}

cat("\n")

# =============================================================================
# Example 5: F0 Contour Extraction
# =============================================================================

cat("Example 5: F0 Contour Extraction\n")
cat("---------------------------------\n")

# Extract F0 contour along with measures
result_f0 <- lst_vat(test_wav,
                    return_f0 = TRUE,
                    verbose = FALSE)

if (result_f0$success && !is.null(result_f0$f0)) {
  f0_values <- unlist(result_f0$f0)
  voiced_f0 <- f0_values[f0_values > 0]

  cat(sprintf("F0 frames: %d total, %d voiced\n",
              length(f0_values), length(voiced_f0)))

  if (length(voiced_f0) > 0) {
    cat(sprintf("F0 statistics:\n"))
    cat(sprintf("  Mean:   %.2f Hz\n", mean(voiced_f0)))
    cat(sprintf("  Median: %.2f Hz\n", median(voiced_f0)))
    cat(sprintf("  SD:     %.2f Hz\n", sd(voiced_f0)))
    cat(sprintf("  Range:  %.2f - %.2f Hz\n", min(voiced_f0), max(voiced_f0)))
  }

  # Create a simple plot if in interactive mode
  if (interactive()) {
    plot(f0_values, type = "l",
         main = "F0 Contour",
         xlab = "Frame",
         ylab = "F0 (Hz)",
         col = "blue")
    abline(h = mean(voiced_f0), col = "red", lty = 2)
  }
}

cat("\n")

# =============================================================================
# Example 6: Batch Processing (if multiple files available)
# =============================================================================

cat("Example 6: Batch Processing\n")
cat("----------------------------\n")

# Get all available test files
test_files <- list.files(
  system.file("samples", "sustained", package = "superassp"),
  pattern = "\\.wav$",
  full.names = TRUE
)

if (length(test_files) > 1) {
  # Limit to first 3 for demo
  test_files <- head(test_files, 3)

  cat(sprintf("Processing %d files...\n", length(test_files)))

  results <- lst_vat(test_files, verbose = FALSE)

  cat("\nResults:\n")
  for (i in seq_along(results)) {
    if (results[[i]]$success) {
      cat(sprintf("  %s: DFA = %.6f, RPDE = %.6f\n",
                  basename(results[[i]]$file),
                  results[[i]]$measures$DFA,
                  results[[i]]$measures$RPDE))
    } else {
      cat(sprintf("  %s: FAILED - %s\n",
                  basename(test_files[i]),
                  results[[i]]$error))
    }
  }

  # Extract specific measure from all files
  dfa_values <- sapply(results, function(r) {
    if (r$success) r$measures$DFA else NA
  })

  cat(sprintf("\nDFA across files: %.6f ± %.6f (mean ± SD)\n",
              mean(dfa_values, na.rm = TRUE),
              sd(dfa_values, na.rm = TRUE)))

} else {
  cat("Only one test file available. Skipping batch processing example.\n")
}

cat("\n")

# =============================================================================
# Summary
# =============================================================================

cat("Examples Complete!\n")
cat("==================\n\n")

cat("For more information:\n")
cat("  ?lst_vat              - Function documentation\n")
cat("  ?install_voice_analysis - Installation help\n")
cat("  voice_analysis_info() - System capabilities\n\n")

cat("Feature categories available:\n")
cat("  - Jitter (22-25 measures)\n")
cat("  - Shimmer (22-25 measures)\n")
cat("  - HNR/NHR (4 measures)\n")
cat("  - DFA, RPDE, PPE (3 measures)\n")
cat("  - GNE (6 measures)\n")
cat("  - Glottal Quotient (3 measures)\n")
cat("  - VFER (7 measures)\n")
cat("  - MFCCs (84 measures)\n")
cat("  - Wavelet (~50 measures)\n")
cat("  - EMD (6 measures)\n\n")

cat("Citation:\n")
cat("  Tsanas, A., Little, M., McSharry, P., & Ramig, L. (2011).\n")
cat("  Nonlinear speech analysis algorithms mapped to a standard metric\n")
cat("  achieve clinically useful quantification of average Parkinson's\n")
cat("  disease symptom severity. Journal of the Royal Society Interface,\n")
cat("  8(59), 842-855.\n")
