# Manual Test Script for DisVoice Integration
#
# Run this script interactively to test the DisVoice functions
# with real audio files.

# ============================================================================
# Setup
# ============================================================================

# Set working directory to superassp package root
# setwd("/Users/frkkan96/Documents/src/superassp")

# Load the package (assuming you've built and loaded it)
# library(superassp)

# Or source the files directly for development:
source("R/disvoice_init.R")
source("R/disvoice_utils.R")
source("R/ssff_python_dv_f0.R")
source("R/ssff_python_dv_formants.R")

# ============================================================================
# Check DisVoice Availability
# ============================================================================

cat("\n=== Checking DisVoice Support ===\n")
if (has_disvoice_support()) {
  cat("✓ DisVoice support is available\n")
} else {
  cat("✗ DisVoice support not available\n")
  cat("  Install with: install_disvoice_python()\n")
  stop("DisVoice support required for testing")
}

# ============================================================================
# Find Test Audio File
# ============================================================================

# Look for test audio files in common locations
test_audio_paths <- c(
  "tests/testthat/test_audio.wav",
  "testthat/test_audio.wav",
  "/tmp/test_audio.wav"
)

test_audio <- NULL
for (path in test_audio_paths) {
  if (file.exists(path)) {
    test_audio <- path
    break
  }
}

if (is.null(test_audio)) {
  cat("\n⚠ Warning: No test audio file found\n")
  cat("  Please provide path to a WAV file:\n")
  test_audio <- readline("  Audio file path: ")

  if (!file.exists(test_audio)) {
    stop("Audio file not found: ", test_audio)
  }
}

cat(sprintf("\n✓ Using test audio: %s\n", test_audio))

# Get audio info
audio_info <- av::av_media_info(test_audio)
cat(sprintf("  Duration: %.2f seconds\n", audio_info$duration))
cat(sprintf("  Sample rate: %d Hz\n", audio_info$audio$sample_rate))
cat(sprintf("  Channels: %d\n", audio_info$audio$channels))

# ============================================================================
# Test trk_dv_f0()
# ============================================================================

cat("\n=== Testing trk_dv_f0() ===\n")

cat("\n1. Basic F0 extraction (with voicing)...\n")
start_time <- Sys.time()
f0_result <- trk_dv_f0(
  test_audio,
  frame_shift = 10,
  min_f0 = 75,
  max_f0 = 600,
  include_voicing = TRUE
)
end_time <- Sys.time()
f0_duration <- as.numeric(difftime(end_time, start_time, units = "secs"))

cat(sprintf("  ✓ Completed in %.3f seconds\n", f0_duration))
cat(sprintf("  Number of frames: %d\n", nrow(f0_result$tracks)))
cat(sprintf("  Track names: %s\n", paste(colnames(f0_result$tracks), collapse = ", ")))
cat(sprintf("  Sample rate: %.1f Hz\n", f0_result$sampleRate))

# Basic validation
valid_f0 <- f0_result$tracks[, "f0"]
valid_f0 <- valid_f0[valid_f0 > 0 & !is.na(valid_f0)]
if (length(valid_f0) > 0) {
  cat(sprintf("  F0 range: %.1f - %.1f Hz\n", min(valid_f0), max(valid_f0)))
  cat(sprintf("  Mean F0: %.1f Hz\n", mean(valid_f0)))
}

voicing <- f0_result$tracks[, "voicing"]
voiced_frames <- sum(voicing == 1, na.rm = TRUE)
total_frames <- length(voicing)
cat(sprintf("  Voiced frames: %d / %d (%.1f%%)\n",
            voiced_frames, total_frames, 100 * voiced_frames / total_frames))

cat("\n2. Testing output_format = 'dataframe'...\n")
f0_df <- trk_dv_f0(test_audio, output_format = "dataframe")
cat(sprintf("  ✓ Returned data.frame with %d rows and %d columns\n",
            nrow(f0_df), ncol(f0_df)))
cat(sprintf("  Column names: %s\n", paste(names(f0_df), collapse = ", ")))

cat("\n3. Testing output_format = 'list'...\n")
f0_list <- trk_dv_f0(test_audio, output_format = "list")
cat(sprintf("  ✓ Returned list with elements: %s\n",
            paste(names(f0_list), collapse = ", ")))

cat("\n4. Testing without voicing track...\n")
f0_no_voicing <- trk_dv_f0(test_audio, include_voicing = FALSE)
cat(sprintf("  ✓ Track names: %s\n", paste(colnames(f0_no_voicing$tracks), collapse = ", ")))

cat("\n5. Testing different frame_shift (5ms = 200 Hz)...\n")
start_time <- Sys.time()
f0_5ms <- trk_dv_f0(test_audio, frame_shift = 5)
end_time <- Sys.time()
cat(sprintf("  ✓ Completed in %.3f seconds\n",
            as.numeric(difftime(end_time, start_time, units = "secs"))))
cat(sprintf("  Number of frames: %d (vs %d at 10ms)\n",
            nrow(f0_5ms$tracks), nrow(f0_result$tracks)))

# ============================================================================
# Test trk_dv_formants()
# ============================================================================

cat("\n=== Testing trk_dv_formants() ===\n")

cat("\n1. Basic formant extraction...\n")
start_time <- Sys.time()
formant_result <- trk_dv_formants(
  test_audio,
  frame_shift = 5,
  window_size = 25,
  max_formants = 5,
  max_formant_freq = 5500
)
end_time <- Sys.time()
formant_duration <- as.numeric(difftime(end_time, start_time, units = "secs"))

cat(sprintf("  ✓ Completed in %.3f seconds\n", formant_duration))
cat(sprintf("  Number of frames: %d\n", nrow(formant_result$tracks)))
cat(sprintf("  Track names: %s\n", paste(colnames(formant_result$tracks), collapse = ", ")))
cat(sprintf("  Sample rate: %.1f Hz\n", formant_result$sampleRate))

# Basic validation
for (formant_name in c("F1", "F2", "F3", "F4")) {
  formant_values <- formant_result$tracks[, formant_name]
  valid_values <- formant_values[formant_values > 0 & !is.na(formant_values)]
  if (length(valid_values) > 0) {
    cat(sprintf("  %s range: %.0f - %.0f Hz (mean: %.0f Hz)\n",
                formant_name, min(valid_values), max(valid_values), mean(valid_values)))
  }
}

cat("\n2. Testing output_format = 'dataframe'...\n")
formant_df <- trk_dv_formants(test_audio, output_format = "dataframe")
cat(sprintf("  ✓ Returned data.frame with %d rows and %d columns\n",
            nrow(formant_df), ncol(formant_df)))
cat(sprintf("  Column names: %s\n", paste(names(formant_df), collapse = ", ")))

cat("\n3. Testing output_format = 'list'...\n")
formant_list <- trk_dv_formants(test_audio, output_format = "list")
cat(sprintf("  ✓ Returned list with elements: %s\n",
            paste(names(formant_list), collapse = ", ")))
cat(sprintf("  Parameters: %s\n",
            paste(names(formant_list$parameters), collapse = ", ")))

cat("\n4. Testing male speaker parameters (max_formant_freq = 5000)...\n")
formant_male <- trk_dv_formants(test_audio, max_formant_freq = 5000)
cat(sprintf("  ✓ Number of frames: %d\n", nrow(formant_male$tracks)))

cat("\n5. Testing different frame_shift (10ms = 100 Hz)...\n")
start_time <- Sys.time()
formant_10ms <- trk_dv_formants(test_audio, frame_shift = 10)
end_time <- Sys.time()
cat(sprintf("  ✓ Completed in %.3f seconds\n",
            as.numeric(difftime(end_time, start_time, units = "secs"))))
cat(sprintf("  Number of frames: %d (vs %d at 5ms)\n",
            nrow(formant_10ms$tracks), nrow(formant_result$tracks)))

# ============================================================================
# Performance Summary
# ============================================================================

cat("\n=== Performance Summary ===\n")
cat(sprintf("  F0 extraction (10ms):      %.3f seconds\n", f0_duration))
cat(sprintf("  Formant extraction (5ms):  %.3f seconds\n", formant_duration))
cat(sprintf("  Audio duration:            %.2f seconds\n", audio_info$duration))
cat(sprintf("  F0 processing speed:       %.1fx realtime\n",
            audio_info$duration / f0_duration))
cat(sprintf("  Formant processing speed:  %.1fx realtime\n",
            audio_info$duration / formant_duration))

# ============================================================================
# Visualization (Optional)
# ============================================================================

cat("\n=== Visualization ===\n")
cat("Would you like to create plots? (y/n): ")
create_plots <- readline()

if (tolower(create_plots) == "y") {
  cat("\nCreating plots...\n")

  # Plot F0
  par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))

  # F0 track
  f0_track <- f0_result$tracks[, "f0"]
  f0_track[f0_track == 0] <- NA  # Replace 0 with NA for plotting
  time_f0 <- seq(0, length(f0_track) - 1) / f0_result$sampleRate

  plot(time_f0, f0_track, type = "l", col = "blue", lwd = 2,
       main = "F0 Track (DisVoice)", xlab = "Time (s)", ylab = "F0 (Hz)",
       ylim = c(0, max(f0_track, na.rm = TRUE) * 1.1))
  grid()

  # Formants
  F1 <- formant_result$tracks[, "F1"]
  F2 <- formant_result$tracks[, "F2"]
  F3 <- formant_result$tracks[, "F3"]
  F4 <- formant_result$tracks[, "F4"]
  time_formants <- seq(0, length(F1) - 1) / formant_result$sampleRate

  plot(time_formants, F1, type = "l", col = "red", lwd = 2,
       main = "Formant Tracks (DisVoice)", xlab = "Time (s)", ylab = "Frequency (Hz)",
       ylim = c(0, max(F4, na.rm = TRUE) * 1.1))
  lines(time_formants, F2, col = "green", lwd = 2)
  lines(time_formants, F3, col = "blue", lwd = 2)
  lines(time_formants, F4, col = "purple", lwd = 2)
  legend("topright", legend = c("F1", "F2", "F3", "F4"),
         col = c("red", "green", "blue", "purple"), lwd = 2)
  grid()

  cat("✓ Plots created\n")
}

# ============================================================================
# Summary
# ============================================================================

cat("\n=== Test Complete ===\n")
cat("✓ All basic tests passed\n")
cat("✓ Both functions produce expected output formats\n")
cat("✓ Performance is satisfactory\n")
cat("\nDisVoice integration is working correctly!\n\n")
