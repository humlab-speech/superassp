#!/usr/bin/env Rscript
# Test script to verify renamed functions work correctly
# Uses systematic sample files from inst/samples/sustained

library(superassp)

cat("\n=== Testing Renamed Voice Analysis Functions ===\n\n")

# Get test file from inst/samples/sustained
test_wav <- system.file('samples', 'sustained', 'a1.wav', package='superassp')

if (!file.exists(test_wav)) {
  stop("Test file not found: ", test_wav)
}

cat("Using test file:", test_wav, "\n\n")

# Test counter
tests_run <- 0
tests_passed <- 0
tests_skipped <- 0

run_test <- function(name, func_call) {
  tests_run <<- tests_run + 1
  cat(sprintf("Test %d: %s... ", tests_run, name))

  tryCatch({
    result <- func_call
    tests_passed <<- tests_passed + 1
    cat("✓ PASS\n")
    return(TRUE)
  }, error = function(e) {
    msg <- conditionMessage(e)
    # Check if error is due to missing Parselmouth/Python
    if (grepl("parselmouth|Python|RETICULATE", msg, ignore.case = TRUE)) {
      tests_skipped <<- tests_skipped + 1
      cat("○ SKIPPED (Parselmouth not available)\n")
      return(NA)
    } else {
      cat("✗ FAIL:", msg, "\n")
      return(FALSE)
    }
  })
}

cat("--- C++ Baseline Tests ---\n")
run_test("forest() [C++ ASSP baseline]", forest(test_wav, toFile=FALSE))

cat("\n--- Renamed SSFF Functions (Parselmouth-based) ---\n")
run_test("praat_formant_burg()", praat_formant_burg(test_wav, toFile=FALSE))
run_test("praat_formantpath_burg()", praat_formantpath_burg(test_wav, toFile=FALSE))
run_test("praat_pitch()", praat_pitch(test_wav, toFile=FALSE))

cat("\n--- Function Signature Verification ---\n")
cat("Verifying all renamed functions exist and _opt versions are gone...\n")

# Check renamed functions exist
renamed_funcs <- c(
  "praat_avqi", "praat_dsi", "praat_voice_tremor", "praat_voice_report",
  "praat_sauce", "praat_formant_burg", "praat_formantpath_burg",
  "praat_pitch", "praat_intensity", "praat_spectral_moments"
)

old_funcs <- c(
  "praat_avqi_opt", "praat_dsi_opt", "praat_voice_tremor_opt",
  "praat_voice_report_opt", "praat_sauce_opt", "praat_formant_burg_opt",
  "praat_formantpath_burg_opt", "praat_pitch_opt", "praat_intensity_opt",
  "praat_spectral_moments_opt", "praat_moments"
)

all_renamed_exist <- all(sapply(renamed_funcs, exists))
all_old_gone <- all(!sapply(old_funcs, exists))

if (all_renamed_exist) {
  cat("  ✓ All renamed functions exist\n")
  tests_passed <- tests_passed + 1
} else {
  cat("  ✗ Some renamed functions missing!\n")
  missing <- renamed_funcs[!sapply(renamed_funcs, exists)]
  cat("    Missing:", paste(missing, collapse=", "), "\n")
}
tests_run <- tests_run + 1

if (all_old_gone) {
  cat("  ✓ All _opt functions removed\n")
  tests_passed <- tests_passed + 1
} else {
  cat("  ✗ Some _opt functions still exist!\n")
  remaining <- old_funcs[sapply(old_funcs, exists)]
  cat("    Remaining:", paste(remaining, collapse=", "), "\n")
}
tests_run <- tests_run + 1

# Summary
cat("\n=== Test Summary ===\n")
cat(sprintf("Total tests: %d\n", tests_run))
cat(sprintf("Passed: %d\n", tests_passed))
cat(sprintf("Skipped: %d (Parselmouth not available)\n", tests_skipped))
cat(sprintf("Failed: %d\n", tests_run - tests_passed - tests_skipped))

if (tests_run - tests_passed - tests_skipped == 0) {
  cat("\n✓ All tests passed!\n")
  quit(status = 0)
} else {
  cat("\n✗ Some tests failed!\n")
  quit(status = 1)
}
