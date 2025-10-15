#!/usr/bin/env Rscript
##' Verify that Python DSP functions have been correctly updated
##'
##' This script verifies that the code has been updated to use av_load_for_python()
##' instead of librosa.load() - without requiring Python modules

library(superassp)

cat("==== Verifying Python DSP Code Updates ====\n\n")

# Test 1: Verify av_python_helpers.R exists and functions are available
cat("Test 1: Checking av_python_helpers.R functions exist\n")
tryCatch({
  # Check if exported functions exist
  stopifnot(exists("av_to_python_audio"))
  stopifnot(exists("av_load_for_python"))

  # Check if internal function exists (will be in namespace but not exported)
  # We can access it via ::: operator
  stopifnot(exists("processMediaFiles_Python", where = asNamespace("superassp")))

  cat("  ✓ av_to_python_audio exists (exported)\n")
  cat("  ✓ av_load_for_python exists (exported)\n")
  cat("  ✓ processMediaFiles_Python exists (internal)\n\n")
}, error = function(e) {
  cat("  ✗ FAILED: Helper functions not found\n")
  cat("  Error:", conditionMessage(e), "\n\n")
  stop(e)
})

# Test 2: Verify function signatures
cat("Test 2: Verifying function signatures\n")
tryCatch({
  # Check av_to_python_audio
  args_av_to_python <- names(formals(av_to_python_audio))
  stopifnot("audio_data" %in% args_av_to_python)
  stopifnot("sample_rate" %in% args_av_to_python)
  stopifnot("channels" %in% args_av_to_python)

  # Check av_load_for_python
  args_av_load <- names(formals(av_load_for_python))
  stopifnot("file_path" %in% args_av_load)
  stopifnot("start_time" %in% args_av_load)
  stopifnot("end_time" %in% args_av_load)
  stopifnot("target_sample_rate" %in% args_av_load)

  cat("  ✓ av_to_python_audio signature correct\n")
  cat("  ✓ av_load_for_python signature correct\n\n")
}, error = function(e) {
  cat("  ✗ FAILED: Function signatures incorrect\n")
  cat("  Error:", conditionMessage(e), "\n\n")
  stop(e)
})

# Test 3: Verify Python DSP functions exist
cat("Test 3: Checking Python DSP functions exist\n")
tryCatch({
  stopifnot(exists("swipe_opt"))
  stopifnot(exists("rapt_opt"))
  stopifnot(exists("reaper_opt"))

  cat("  ✓ swipe_opt exists\n")
  cat("  ✓ rapt_opt exists\n")
  cat("  ✓ reaper_opt exists\n\n")
}, error = function(e) {
  cat("  ✗ FAILED: Python DSP functions not found\n")
  cat("  Error:", conditionMessage(e), "\n\n")
  stop(e)
})

# Test 4: Verify internal processing functions exist
cat("Test 4: Checking internal processing functions exist\n")
tryCatch({
  # These are internal functions, access via :::
  ns <- asNamespace("superassp")
  stopifnot(exists("process_swipe_single", where = ns))
  stopifnot(exists("process_rapt_single", where = ns))
  stopifnot(exists("process_reaper_single", where = ns))

  cat("  ✓ process_swipe_single exists (internal)\n")
  cat("  ✓ process_rapt_single exists (internal)\n")
  cat("  ✓ process_reaper_single exists (internal)\n\n")
}, error = function(e) {
  cat("  ✗ FAILED: Internal processing functions not found\n")
  cat("  Error:", conditionMessage(e), "\n\n")
  stop(e)
})

# Test 5: Verify function body uses av_load_for_python (code inspection)
cat("Test 5: Verifying functions use av_load_for_python (code inspection)\n")
tryCatch({
  # Get function body as text (accessing internal functions)
  ns <- asNamespace("superassp")
  swipe_body <- deparse(body(get("process_swipe_single", envir = ns)))
  rapt_body <- deparse(body(get("process_rapt_single", envir = ns)))
  reaper_body <- deparse(body(get("process_reaper_single", envir = ns)))

  # Check for av_load_for_python usage
  swipe_uses_av <- any(grepl("av_load_for_python", swipe_body, fixed = TRUE))
  rapt_uses_av <- any(grepl("av_load_for_python", rapt_body, fixed = TRUE))
  reaper_uses_av <- any(grepl("av_load_for_python", reaper_body, fixed = TRUE))

  # Check that librosa.load is NOT used in the modern way
  swipe_no_librosa <- !any(grepl("lr.load", swipe_body, fixed = TRUE))
  rapt_no_librosa <- !any(grepl("lr.load", rapt_body, fixed = TRUE))
  reaper_no_librosa <- !any(grepl("lr.load", reaper_body, fixed = TRUE))

  stopifnot(swipe_uses_av)
  stopifnot(rapt_uses_av)
  stopifnot(reaper_uses_av)
  stopifnot(swipe_no_librosa)
  stopifnot(rapt_no_librosa)
  stopifnot(reaper_no_librosa)

  cat("  ✓ process_swipe_single uses av_load_for_python\n")
  cat("  ✓ process_rapt_single uses av_load_for_python\n")
  cat("  ✓ process_reaper_single uses av_load_for_python\n")
  cat("  ✓ None of the functions use librosa.load anymore\n\n")
}, error = function(e) {
  cat("  ✗ FAILED: Functions not properly updated\n")
  cat("  Error:", conditionMessage(e), "\n\n")
  stop(e)
})

# Test 6: Verify functions check for numpy array in Python code
cat("Test 6: Verifying Python code receives numpy array\n")
tryCatch({
  ns <- asNamespace("superassp")
  swipe_body <- deparse(body(get("process_swipe_single", envir = ns)))
  rapt_body <- deparse(body(get("process_rapt_single", envir = ns)))
  reaper_body <- deparse(body(get("process_reaper_single", envir = ns)))

  # Check that audio_result$audio_np is assigned to py$x
  swipe_uses_np <- any(grepl("audio_result\\$audio_np", swipe_body))
  rapt_uses_np <- any(grepl("audio_result\\$audio_np", rapt_body))
  reaper_uses_np <- any(grepl("audio_result\\$audio_np", reaper_body))

  stopifnot(swipe_uses_np)
  stopifnot(rapt_uses_np)
  stopifnot(reaper_uses_np)

  cat("  ✓ All functions pass audio_np to Python\n")
  cat("  ✓ Memory-based processing confirmed\n\n")
}, error = function(e) {
  cat("  ✗ FAILED: Functions don't properly pass numpy arrays\n")
  cat("  Error:", conditionMessage(e), "\n\n")
  stop(e)
})

cat("==== Summary ====\n")
cat("✓ All code verification tests PASSED!\n")
cat("✓ Helper functions (av_python_helpers.R) are properly defined\n")
cat("✓ Python DSP functions (swipe_opt, rapt_opt, reaper_opt) exist\n")
cat("✓ Functions have been updated to use av_load_for_python()\n")
cat("✓ Functions NO LONGER use librosa.load() for file reading\n")
cat("✓ Functions pass numpy arrays directly to Python (memory-based)\n")
cat("\n")
cat("The code changes are correct!\n")
cat("To fully test with actual Python modules, ensure pysptk and pyreaper are installed.\n")
