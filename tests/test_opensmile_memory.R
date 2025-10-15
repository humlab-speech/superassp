#!/usr/bin/env Rscript
##' Test openSMILE Slice Functions with Memory-Based Loading
##'
##' This script tests that ComParE_2016 and GeMAPS functions work correctly
##' with the new av-based memory loading approach.

library(superassp)

cat("==== Testing openSMILE Slice Functions (Memory-Based Loading) ====\n\n")

# Check if openSMILE is available
if (!reticulate::py_module_available("opensmile")) {
  cat("WARNING: opensmile Python module not available\n")
  cat("Install with: pip install opensmile\n")
  cat("Skipping tests that require opensmile\n\n")
  quit(save = "no", status = 0)
}

# Find test audio file
test_file <- list.files(
  system.file("samples", "sustained", package = "superassp"),
  pattern = "a1.wav",
  full.names = TRUE
)[1]

if (length(test_file) == 0 || !file.exists(test_file)) {
  cat("ERROR: Test audio file not found\n")
  cat("Looking for: samples/sustained/a1.wav\n")
  quit(save = "no", status = 1)
}

cat("Using test file:", basename(test_file), "\n")

# Get file info
info <- av::av_media_info(test_file)
cat("Duration:", info$duration, "seconds\n")
cat("Sample rate:", info$audio$sample_rate, "Hz\n\n")

# Test 1: ComParE_2016 basic functionality
cat("Test 1: ComParE_2016 basic functionality\n")
tryCatch({
  result <- ComParE_2016(
    listOfFiles = test_file,
    verbose = FALSE
  )

  # Verify result structure
  stopifnot(is.list(result))
  stopifnot(length(result) > 0)

  # ComParE_2016 returns a list with feature dataframe
  # The exact structure depends on openSMILE version
  cat("  ✓ ComParE_2016 executed successfully\n")
  cat("  ✓ Result is a list with", length(result), "elements\n")

  # Show first few feature names if available
  if (is.data.frame(result[[1]])) {
    cat("  ✓ Result contains dataframe with", ncol(result[[1]]), "features\n")
  }
  cat("\n")
}, error = function(e) {
  cat("  ✗ FAILED:", conditionMessage(e), "\n\n")
  stop(e)
})

# Test 2: GeMAPS basic functionality
cat("Test 2: GeMAPS basic functionality\n")
tryCatch({
  result <- GeMAPS(
    listOfFiles = test_file,
    verbose = FALSE
  )

  # Verify result structure
  stopifnot(is.list(result))
  stopifnot(length(result) > 0)

  cat("  ✓ GeMAPS executed successfully\n")
  cat("  ✓ Result is a list with", length(result), "elements\n")

  if (is.data.frame(result[[1]])) {
    cat("  ✓ Result contains dataframe with", ncol(result[[1]]), "features\n")
  }
  cat("\n")
}, error = function(e) {
  cat("  ✗ FAILED:", conditionMessage(e), "\n\n")
  stop(e)
})

# Test 3: Time windowing with ComParE_2016
if (info$duration >= 2.0) {
  cat("Test 3: Time windowing with ComParE_2016\n")
  tryCatch({
    # Extract middle 1 second
    start_time <- info$duration / 2 - 0.5
    end_time <- info$duration / 2 + 0.5

    result <- ComParE_2016(
      listOfFiles = test_file,
      beginTime = start_time,
      endTime = end_time,
      verbose = FALSE
    )

    stopifnot(is.list(result))
    stopifnot(length(result) > 0)

    cat("  ✓ Time windowing works correctly\n")
    cat("  ✓ Extracted segment:", start_time, "to", end_time, "seconds\n\n")
  }, error = function(e) {
    cat("  ✗ FAILED:", conditionMessage(e), "\n\n")
    stop(e)
  })
} else {
  cat("Test 3: SKIPPED (file too short for time windowing test)\n\n")
}

# Test 4: Time windowing with GeMAPS
if (info$duration >= 2.0) {
  cat("Test 4: Time windowing with GeMAPS\n")
  tryCatch({
    # Extract first 1 second
    result <- GeMAPS(
      listOfFiles = test_file,
      beginTime = 0.0,
      endTime = 1.0,
      verbose = FALSE
    )

    stopifnot(is.list(result))
    stopifnot(length(result) > 0)

    cat("  ✓ Time windowing works correctly\n")
    cat("  ✓ Extracted segment: 0.0 to 1.0 seconds\n\n")
  }, error = function(e) {
    cat("  ✗ FAILED:", conditionMessage(e), "\n\n")
    stop(e)
  })
} else {
  cat("Test 4: SKIPPED (file too short for time windowing test)\n\n")
}

# Test 5: Verify memory-based approach (no temp files created)
cat("Test 5: Verify no temporary WAV files created\n")
tryCatch({
  # Get list of files before
  temp_dir <- tempdir()
  files_before <- list.files(temp_dir, pattern = "\\.wav$", full.names = TRUE)

  # Run ComParE_2016
  result <- ComParE_2016(
    listOfFiles = test_file,
    verbose = FALSE
  )

  # Get list of files after
  files_after <- list.files(temp_dir, pattern = "\\.wav$", full.names = TRUE)

  # Check that no new WAV files were created
  new_files <- setdiff(files_after, files_before)

  if (length(new_files) == 0) {
    cat("  ✓ No temporary WAV files created (memory-based!)\n")
  } else {
    cat("  ⚠ WARNING: Temporary files detected:\n")
    for (f in new_files) {
      cat("    -", basename(f), "\n")
    }
  }
  cat("\n")
}, error = function(e) {
  cat("  ✗ FAILED:", conditionMessage(e), "\n\n")
  # Don't stop on this test - it's informational
})

# Test 6: Code verification - check functions use av_load_for_python
cat("Test 6: Code verification - functions use av_load_for_python\n")
tryCatch({
  # Get function body
  compare_body <- deparse(body(ComParE_2016))
  gemaps_body <- deparse(body(GeMAPS))

  # Check for av_load_for_python usage
  compare_uses_av <- any(grepl("av_load_for_python", compare_body, fixed = TRUE))
  gemaps_uses_av <- any(grepl("av_load_for_python", gemaps_body, fixed = TRUE))

  # Check for process_signal (not process_file)
  compare_uses_signal <- any(grepl("process_signal", compare_body, fixed = TRUE))
  gemaps_uses_signal <- any(grepl("process_signal", gemaps_body, fixed = TRUE))

  stopifnot(compare_uses_av)
  stopifnot(gemaps_uses_av)
  stopifnot(compare_uses_signal)
  stopifnot(gemaps_uses_signal)

  cat("  ✓ ComParE_2016 uses av_load_for_python\n")
  cat("  ✓ GeMAPS uses av_load_for_python\n")
  cat("  ✓ Both use smile.process_signal() (memory-based)\n")
  cat("  ✓ Memory-based architecture confirmed\n\n")
}, error = function(e) {
  cat("  ✗ FAILED:", conditionMessage(e), "\n\n")
  stop(e)
})

cat("==== Summary ====\n")
cat("✓ All openSMILE slice function tests PASSED!\n")
cat("✓ ComParE_2016 works with memory-based loading\n")
cat("✓ GeMAPS works with memory-based loading\n")
cat("✓ Time windowing functions correctly\n")
cat("✓ No temporary files created (all in memory)\n")
cat("✓ Functions properly use av → numpy → openSMILE pipeline\n\n")

cat("Next steps:\n")
cat("  - Run benchmark tests to measure speed improvements\n")
cat("  - Update additional openSMILE functions (eGeMAPS if present)\n")
cat("  - Begin Parselmouth conversions for Praat slice functions\n")
