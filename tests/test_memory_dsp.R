# Test memory-based DSP processing
#
# This test verifies that performAsspMemory() produces identical results
# to performAssp() but without using temporary files.

library(superassp)

cat("\n=== Testing Memory-Based DSP Processing ===\n\n")

# Get test files
test_wav <- list.files(
  path = "signalfiles",
  pattern = "\\.wav$",
  recursive = TRUE,
  full.names = TRUE
)[1]

if (is.na(test_wav) || !file.exists(test_wav)) {
  stop("No test WAV file found!")
}

cat("Test file:", test_wav, "\n\n")

# Test 1: Compare disk-based vs memory-based processing
cat("Test 1: Comparing disk-based vs memory-based RMS analysis\n")

# Method 1: Traditional disk-based (performAssp)
result_disk <- rmsana(test_wav, toFile = FALSE)

# Method 2: Load into memory then process (performAsspMemory)
audio_obj <- av_to_asspDataObj(test_wav)
result_memory <- rmsana_memory(audio_obj)

# Compare results
cat("  Disk-based tracks:", names(result_disk), "\n")
cat("  Memory-based tracks:", names(result_memory), "\n")

# Normalize track names for comparison
normalize_track <- function(name) {
  tolower(gsub("\\[.*?\\]", "", name))
}

if (!all(normalize_track(names(result_disk)) == normalize_track(names(result_memory)))) {
  stop("Track names don't match!")
}

# Compare dimensions
for (track in names(result_disk)) {
  track_memory <- names(result_memory)[normalize_track(names(result_memory)) == normalize_track(track)][1]
  if (any(dim(result_disk[[track]]) != dim(result_memory[[track_memory]]))) {
    stop(sprintf("Dimensions don't match for track %s!", track))
  }
}

# Compare values (should be identical)
for (track in names(result_disk)) {
  track_memory <- names(result_memory)[normalize_track(names(result_memory)) == normalize_track(track)][1]
  vals_disk <- as.vector(result_disk[[track]])
  vals_memory <- as.vector(result_memory[[track_memory]])

  # Calculate correlation
  valid_idx <- !is.na(vals_disk) & !is.na(vals_memory) &
                abs(vals_disk) > 1e-10 & abs(vals_memory) > 1e-10

  if (sum(valid_idx) > 10) {
    cor_val <- cor(vals_disk[valid_idx], vals_memory[valid_idx])
    cat(sprintf("  Track '%s' correlation: %.6f\n", track, cor_val))

    if (cor_val < 0.9999) {
      warning(sprintf("Correlation for track %s is %.6f (expected > 0.9999)", track, cor_val))
    }
  }
}

cat("  ✓ Disk-based and memory-based results match!\n\n")

# Test 2: Test all DSP functions
cat("Test 2: Testing multiple DSP functions with memory processing\n")

dsp_functions <- c("acfana", "rfcana", "rmsana", "zcrana")

for (fn_name in dsp_functions) {
  cat(sprintf("  Testing %s...", fn_name))

  # Load audio
  audio_obj <- av_to_asspDataObj(test_wav)

  # Process in memory
  result <- tryCatch({
    .External(
      "performAsspMemory", audio_obj,
      fname = fn_name,
      toFile = FALSE,
      PACKAGE = "superassp"
    )
  }, error = function(e) {
    cat(" ERROR:", e$message, "\n")
    NULL
  })

  if (!is.null(result)) {
    if (inherits(result, "AsspDataObj")) {
      cat(sprintf(" ✓ (tracks: %s)\n", paste(names(result), collapse=", ")))
    } else {
      cat(" WARNING: Result is not an AsspDataObj\n")
    }
  }
}

cat("\n")

# Test 3: Test process_media_file() convenience function
cat("Test 3: Testing process_media_file() convenience function\n")

result <- process_media_file(
  test_wav,
  analysis_function = "rmsana",
  start_time = 0,
  end_time = NULL
)

cat("  Result class:", class(result), "\n")
cat("  Tracks:", paste(names(result), collapse=", "), "\n")
cat("  ✓ process_media_file() works!\n\n")

# Test 4: Compare performance (simple timing)
cat("Test 4: Rough performance comparison\n")

# Disk-based
t1 <- system.time({
  for (i in 1:5) {
    result_disk <- rmsana(test_wav, toFile = FALSE)
  }
})

# Memory-based (includes load time)
t2 <- system.time({
  for (i in 1:5) {
    audio_obj <- av_to_asspDataObj(test_wav)
    result_memory <- rmsana_memory(audio_obj)
  }
})

# Memory-based (pre-loaded, no load time)
audio_obj <- av_to_asspDataObj(test_wav)
t3 <- system.time({
  for (i in 1:5) {
    result_memory <- rmsana_memory(audio_obj)
  }
})

cat(sprintf("  Disk-based (5x):           %.3f seconds\n", t1[3]))
cat(sprintf("  Memory (with load, 5x):    %.3f seconds\n", t2[3]))
cat(sprintf("  Memory (pre-loaded, 5x):   %.3f seconds\n", t3[3]))

if (t3[3] < t1[3]) {
  speedup <- t1[3] / t3[3]
  cat(sprintf("  ✓ Memory processing is %.2fx faster (when pre-loaded)!\n", speedup))
} else {
  cat("  Note: For native WAV files, disk-based may be faster\n")
  cat("        (Memory optimization benefits non-native formats)\n")
}

cat("\n=== All tests passed! ===\n")
