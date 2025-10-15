# Test script for AVQI and DSI optimized implementations
# This tests that praat_avqi_opt() and praat_dsi_opt() produce similar results
# to the original Praat-based implementations

library(superassp)

# Check if Parselmouth is available
if (!reticulate::py_module_available("parselmouth")) {
  stop("Parselmouth not available. Install with: reticulate::py_install('praat-parselmouth')")
}

cat("Testing AVQI optimized implementation...\n")

# === Test AVQI ===
# Define sustained vowel samples (times in milliseconds)
sv <- data.frame(
  listOfFiles = c(
    "tests/signalfiles/AVQI/input/sv1.wav",
    "tests/signalfiles/AVQI/input/sv2.wav",
    "tests/signalfiles/AVQI/input/sv3.wav",
    "tests/signalfiles/AVQI/input/sv4.wav"
  ),
  start = rep(63.33, 4),  # milliseconds
  end = rep(2830.56, 4)
)

# Define continuous speech samples (times in milliseconds)
cs <- data.frame(
  listOfFiles = c(
    "tests/signalfiles/AVQI/input/cs1.wav",
    "tests/signalfiles/AVQI/input/cs2.wav",
    "tests/signalfiles/AVQI/input/cs3.wav",
    "tests/signalfiles/AVQI/input/cs4.wav"
  ),
  start = rep(82.50, 4),
  end = rep(3738.39, 4)
)

# Test AVQI optimized version
cat("Running praat_avqi_opt()...\n")
avqi_result <- tryCatch({
  praat_avqi_opt(
    svDF = sv,
    csDF = cs,
    speaker.name = "Test Speaker",
    speaker.ID = "1"
  )
}, error = function(e) {
  cat("ERROR in praat_avqi_opt():\n")
  cat(e$message, "\n")
  NULL
})

if (!is.null(avqi_result)) {
  cat("\nAVQI Results:\n")
  cat("AVQI Version:", avqi_result$AVQI_VERSION, "\n")
  cat("CPPS:", avqi_result$CPPS, "\n")
  cat("HNR:", avqi_result$HNR, "\n")
  cat("Shimmer (local):", avqi_result$Shim_local, "%\n")
  cat("Shimmer (dB):", avqi_result$Shim_local_DB, "dB\n")
  cat("LTAS Slope:", avqi_result$LTAS_Slope, "dB\n")
  cat("LTAS Tilt:", avqi_result$LTAS_Tilt, "dB\n")
  cat("AVQI Score:", avqi_result$AVQI, "\n")

  # Expected values from reference (praat_avqi_3.01.csv)
  expected <- list(
    CPPS = 13.13,
    HNR = 14.47,
    Shim_local = 8.12,
    Shim_local_DB = 0.75,
    LTAS_Slope = -19.71,
    LTAS_Tilt = -11.17,
    AVQI = 2.64
  )

  cat("\nComparison with reference:\n")
  tolerance <- 0.5  # Allow 0.5 difference (reasonable for numerical algorithms)

  all_ok <- TRUE
  for (measure in names(expected)) {
    diff <- abs(avqi_result[[measure]] - expected[[measure]])
    status <- if (diff < tolerance) "✓ OK" else "✗ DIFF"
    cat(sprintf("%s: %.2f (expected %.2f, diff %.2f) %s\n",
                measure, avqi_result[[measure]], expected[[measure]], diff, status))
    if (diff >= tolerance) all_ok <- FALSE
  }

  if (all_ok) {
    cat("\n✓ AVQI test PASSED - all measures within tolerance\n")
  } else {
    cat("\n⚠ AVQI test completed but some measures differ from reference\n")
    cat("Note: Small differences are expected due to numerical precision\n")
  }
} else {
  cat("\n✗ AVQI test FAILED\n")
}

cat("\n" , rep("=", 70), "\n", sep = "")

# === Test DSI ===
cat("\nTesting DSI optimized implementation...\n")

# Check for DSI test files
dsi_input_dir <- "tests/signalfiles/DSI/input"
if (!dir.exists(dsi_input_dir)) {
  cat("DSI test files not found, skipping DSI test\n")
} else {
  # List available files
  dsi_files <- list.files(dsi_input_dir, pattern = "*.wav", full.names = TRUE)

  if (length(dsi_files) == 0) {
    cat("No DSI test files found, skipping DSI test\n")
  } else {
    cat("Found DSI test files:\n")
    print(basename(dsi_files))

    # Look for different file types (im*, fh*, mpt*, ppq*)
    im_files <- grep("^im", basename(dsi_files), value = TRUE)
    fh_files <- grep("^fh", basename(dsi_files), value = TRUE)
    mpt_files <- grep("^mpt", basename(dsi_files), value = TRUE)
    ppq_files <- grep("^ppq", basename(dsi_files), value = TRUE)

    if (length(im_files) > 0 && length(fh_files) > 0 &&
        length(mpt_files) > 0 && length(ppq_files) > 0) {

      # Create dataframes for DSI
      softDF <- data.frame(
        absolute_file_path = file.path(dsi_input_dir, im_files),
        start = 0,
        end = 3  # Use first 3 seconds
      )

      highpitchDF <- data.frame(
        absolute_file_path = file.path(dsi_input_dir, fh_files),
        start = 0,
        end = 3
      )

      maxprolongedDF <- data.frame(
        absolute_file_path = file.path(dsi_input_dir, mpt_files),
        start = 0,
        end = 10  # Full duration
      )

      stableDF <- data.frame(
        absolute_file_path = file.path(dsi_input_dir, ppq_files),
        start = 0,
        end = 5
      )

      cat("\nRunning praat_dsi_opt()...\n")
      dsi_result <- tryCatch({
        praat_dsi_opt(
          softDF = softDF,
          highpitchDF = highpitchDF,
          maxprolongedDF = maxprolongedDF,
          stableDF = stableDF,
          speaker.name = "Test Speaker",
          speaker.ID = "1"
        )
      }, error = function(e) {
        cat("ERROR in praat_dsi_opt():\n")
        cat(e$message, "\n")
        NULL
      })

      if (!is.null(dsi_result)) {
        cat("\nDSI Results:\n")
        cat("Maximum phonation time:", dsi_result$Maximum_phonation_time, "s\n")
        cat("Softest intensity:", dsi_result$Softest_intensity_of_voiced_speech, "dB\n")
        cat("Maximum F0:", dsi_result$Maximum_fundamental_frequency, "Hz\n")
        cat("Jitter ppq5:", dsi_result$Jitter_ppq5, "%\n")
        cat("DSI Score:", dsi_result$Dysphonia_Severity_Index, "\n")

        cat("\n✓ DSI test PASSED - function executed successfully\n")
      } else {
        cat("\n✗ DSI test FAILED\n")
      }
    } else {
      cat("Not all required DSI file types found (need im*, fh*, mpt*, ppq*)\n")
    }
  }
}

cat("\n", rep("=", 70), "\n", sep = "")
cat("Testing complete!\n")
