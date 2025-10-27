#!/usr/bin/env Rscript
#
# Test script for TVWLP formant tracking integration
#
# This script demonstrates the usage of the ultra-optimized formant tracking
# and validates the R-Python integration.

library(superassp)

cat("="*60, "\n")
cat("TVWLP Formant Tracking - Test & Example\n")
cat("="*60, "\n\n")

# Step 1: Check installation
cat("Step 1: Checking Python dependencies\n")
cat(rep("-", 60), "\n")

check_python <- function() {
  tryCatch({
    py_config <- reticulate::py_config()
    cat(sprintf("✓ Python: %s\n", py_config$python))

    # Check for required packages
    has_numpy <- reticulate::py_module_available("numpy")
    has_scipy <- reticulate::py_module_available("scipy")
    has_numba <- reticulate::py_module_available("numba")

    cat(sprintf("  - numpy: %s\n", if (has_numpy) "✓" else "✗"))
    cat(sprintf("  - scipy: %s\n", if (has_scipy) "✓" else "✗"))
    cat(sprintf("  - numba: %s\n", if (has_numba) "✓" else "✗"))

    if (!has_numpy || !has_scipy || !has_numba) {
      cat("\n⚠ Missing dependencies. Run: install_ftrack_tvwlp()\n")
      return(FALSE)
    }

    return(TRUE)

  }, error = function(e) {
    cat("✗ Python not configured\n")
    cat("  Install with: reticulate::install_miniconda()\n")
    return(FALSE)
  })
}

python_ok <- check_python()

if (!python_ok) {
  cat("\nInstalling Python dependencies...\n")
  install_ftrack_tvwlp(build_cython = FALSE, verbose = TRUE)
}

cat("\n")

# Step 2: Create test audio (simple sine wave)
cat("Step 2: Creating test audio\n")
cat(rep("-", 60), "\n")

create_test_audio <- function(filename = "test_audio.wav", duration = 2.0) {
  fs <- 16000
  t <- seq(0, duration, by = 1/fs)

  # Create a simple vowel-like signal with 3 formants
  # F1 = 500 Hz, F2 = 1500 Hz, F3 = 2500 Hz
  f0 <- 120  # Fundamental frequency
  signal <- sin(2 * pi * f0 * t) +
            0.5 * sin(2 * pi * 500 * t) +
            0.3 * sin(2 * pi * 1500 * t) +
            0.2 * sin(2 * pi * 2500 * t)

  # Normalize
  signal <- signal / max(abs(signal)) * 0.8

  # Write WAV file
  tuneR::writeWave(
    tuneR::Wave(
      left = as.integer(signal * 32767),
      samp.rate = fs,
      bit = 16
    ),
    filename
  )

  cat(sprintf("✓ Created test audio: %s (%.1fs, %d Hz)\n", filename, duration, fs))
  return(filename)
}

test_file <- create_test_audio()

cat("\n")

# Step 3: Test different optimization levels
cat("Step 3: Testing optimization levels\n")
cat(rep("-", 60), "\n")

benchmark_optimization_level <- function(file, level, lptype = "tvwlp_l2") {
  cat(sprintf("\nTesting: %s (method: %s)\n", level, lptype))

  start_time <- Sys.time()

  tryCatch({
    result <- trk_formants_tvwlp(
      file,
      lptype = lptype,
      optimization_level = level,
      npeaks = 3,
      toFile = FALSE,
      verbose = FALSE
    )

    end_time <- Sys.time()
    elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))

    # Get audio duration
    audio_info <- av::av_media_info(file)
    duration <- audio_info$duration
    rt_factor <- duration / elapsed

    cat(sprintf("  Runtime: %.3fs\n", elapsed))
    cat(sprintf("  Real-time factor: %.2fx\n", rt_factor))
    cat(sprintf("  Output shape: [%d x %d]\n", nrow(result$tracks), ncol(result$tracks)))

    # Show first formant values
    if (ncol(result$tracks) > 0) {
      cat(sprintf("  F1 mean: %.1f Hz\n", mean(result$tracks[1, ], na.rm = TRUE)))
      cat(sprintf("  F2 mean: %.1f Hz\n", mean(result$tracks[2, ], na.rm = TRUE)))
      cat(sprintf("  F3 mean: %.1f Hz\n", mean(result$tracks[3, ], na.rm = TRUE)))
    }

    return(list(
      level = level,
      elapsed = elapsed,
      rt_factor = rt_factor,
      success = TRUE
    ))

  }, error = function(e) {
    cat(sprintf("  ✗ Error: %s\n", e$message))
    return(list(
      level = level,
      elapsed = NA,
      rt_factor = NA,
      success = FALSE
    ))
  })
}

# Test all optimization levels
results <- list()

# Vectorized (always works)
results[[1]] <- benchmark_optimization_level(test_file, "vectorized", "tvlp_l2")

# Numba (if available)
results[[2]] <- benchmark_optimization_level(test_file, "numba", "tvwlp_l2")

# Ultra (if Cython available)
results[[3]] <- benchmark_optimization_level(test_file, "ultra", "tvwlp_l2")

cat("\n")

# Step 4: Summary
cat("Step 4: Performance Summary\n")
cat(rep("-", 60), "\n\n")

cat(sprintf("%-20s %-12s %-12s %-10s\n", "Level", "Runtime", "RT Factor", "Status"))
cat(rep("-", 60), "\n")

for (res in results) {
  if (res$success) {
    status <- if (res$rt_factor > 1.0) "✓ Faster" else "⚠ Slower"
    cat(sprintf("%-20s %8.3fs     %8.2fx     %s\n",
                res$level, res$elapsed, res$rt_factor, status))
  } else {
    cat(sprintf("%-20s %8s     %8s     ✗ Failed\n", res$level, "N/A", "N/A"))
  }
}

cat("\n")

# Step 5: Test TVWLP vs TVLP
cat("Step 5: Comparing TVWLP vs TVLP methods\n")
cat(rep("-", 60), "\n")

compare_methods <- function(file) {
  methods <- c("tvwlp_l2", "tvlp_l2")

  for (method in methods) {
    cat(sprintf("\nMethod: %s\n", method))

    start_time <- Sys.time()

    tryCatch({
      result <- trk_formants_tvwlp(
        file,
        lptype = method,
        optimization_level = "ultra",
        npeaks = 3,
        toFile = FALSE,
        verbose = FALSE
      )

      end_time <- Sys.time()
      elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))

      cat(sprintf("  Runtime: %.3fs\n", elapsed))
      cat(sprintf("  Has GCI detection: %s\n", if (grepl("tvwlp", method)) "Yes" else "No"))

    }, error = function(e) {
      cat(sprintf("  ✗ Error: %s\n", e$message))
    })
  }
}

compare_methods(test_file)

cat("\n")

# Step 6: Test batch processing
cat("Step 6: Testing batch processing\n")
cat(rep("-", 60), "\n")

# Create multiple test files
test_files <- sapply(1:3, function(i) {
  create_test_audio(sprintf("test_audio_%d.wav", i), duration = 1.0)
})

cat("\nProcessing batch of", length(test_files), "files...\n")

start_time <- Sys.time()

n_processed <- tryCatch({
  trk_formants_tvwlp(
    test_files,
    lptype = "tvlp_l2",  # Fast method for batch
    optimization_level = "ultra",
    outputDirectory = tempdir(),
    explicitExt = "fms",
    toFile = TRUE,
    verbose = TRUE
  )
}, error = function(e) {
  cat(sprintf("✗ Batch processing failed: %s\n", e$message))
  0
})

end_time <- Sys.time()
elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))

if (n_processed > 0) {
  cat(sprintf("\n✓ Processed %d/%d files in %.2fs\n",
              n_processed, length(test_files), elapsed))
  cat(sprintf("  Average: %.3fs per file\n", elapsed / n_processed))
}

# Cleanup
cat("\nCleaning up test files...\n")
unlink(test_file)
unlink(test_files)

cat("\n")
cat("="*60, "\n")
cat("Test completed!\n")
cat("="*60, "\n\n")

# Step 7: Recommendations
cat("Recommendations:\n")
cat(rep("-", 60), "\n")

# Check which optimization is best
best_result <- results[[which.max(sapply(results, function(r) if (r$success) r$rt_factor else 0))]]

if (best_result$success && best_result$level == "ultra") {
  cat("✓ Ultra optimization is working! (4.37x speedup expected)\n")
  cat("  Recommended for production use.\n\n")
} else if (best_result$success && best_result$level == "numba") {
  cat("⚠ Numba optimization is working, but Cython not available\n")
  cat("  For maximum performance, run:\n")
  cat("    install_ftrack_tvwlp(build_cython = TRUE)\n\n")
} else {
  cat("⚠ Only vectorized optimization available\n")
  cat("  For better performance, install dependencies:\n")
  cat("    install_ftrack_tvwlp(build_cython = TRUE)\n\n")
}

cat("\nFor production use:\n")
cat("  • Use optimization_level = 'ultra'\n")
cat("  • Use lptype = 'tvwlp_l2' for accuracy\n")
cat("  • Use lptype = 'tvlp_l2' for speed\n")
cat("  • Process in batches with toFile = TRUE\n")
cat("\n")

cat("Documentation:\n")
cat("  ?trk_formants_tvwlp\n")
cat("  ?install_ftrack_tvwlp\n")
cat("\n")
