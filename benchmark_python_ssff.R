#!/usr/bin/env Rscript
# Benchmark and test script for python_ssff optimized implementations
# This script compares the original and optimized implementations

library(superassp)
library(reticulate)

# Configure Python environment
Sys.setenv(RETICULATE_PYTHON = "/opt/miniconda3/bin/python3")

cat("=== Python SSFF Optimization Benchmark ===\n\n")

# Test file
test_file <- "tests/signalfiles/AVQI/input/sv1.wav"

if(!file.exists(test_file)) {
  stop("Test file not found: ", test_file)
}

cat("Test file:", test_file, "\n")
cat("File size:", file.size(test_file), "bytes\n\n")

# Check Python dependencies
cat("Checking Python dependencies...\n")
py_deps <- c("pysptk", "librosa", "numpy")
for(dep in py_deps) {
  available <- tryCatch({
    py_module_available(dep)
  }, error = function(e) FALSE)
  
  cat(sprintf("  %s: %s\n", dep, ifelse(available, "✓ Available", "✗ Missing")))
}
cat("\n")

# Function to run benchmarks
run_benchmark <- function(func_name, func, ...) {
  cat(sprintf("Benchmarking %s...\n", func_name))
  
  # Warm-up run
  tryCatch({
    suppressMessages(func(..., toFile = FALSE, verbose = FALSE))
  }, error = function(e) {
    cat("  Warm-up failed:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  # Benchmark runs
  times <- numeric(5)
  for(i in 1:5) {
    times[i] <- system.time({
      suppressMessages(func(..., toFile = FALSE, verbose = FALSE))
    })[3]
  }
  
  cat(sprintf("  Mean time: %.3f s (± %.3f)\n", mean(times), sd(times)))
  cat(sprintf("  Min time:  %.3f s\n", min(times)))
  cat(sprintf("  Max time:  %.3f s\n\n", max(times)))
  
  return(list(mean = mean(times), sd = sd(times), min = min(times), max = max(times)))
}

# Test consistency between implementations
test_consistency <- function(original_func, optimized_func, func_name, ...) {
  cat(sprintf("Testing consistency for %s...\n", func_name))
  
  orig_result <- tryCatch({
    suppressMessages(original_func(..., toFile = FALSE, verbose = FALSE))
  }, error = function(e) {
    cat("  Original function failed:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  opt_result <- tryCatch({
    suppressMessages(optimized_func(..., toFile = FALSE, verbose = FALSE))
  }, error = function(e) {
    cat("  Optimized function failed:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  if(is.null(orig_result) || is.null(opt_result)) {
    cat("  ✗ Cannot compare - one or both functions failed\n\n")
    return(FALSE)
  }
  
  # Compare track names
  orig_tracks <- names(orig_result)
  opt_tracks <- names(opt_result)
  
  if(!identical(orig_tracks, opt_tracks)) {
    cat("  ✗ Track names differ\n")
    cat("    Original:", paste(orig_tracks, collapse=", "), "\n")
    cat("    Optimized:", paste(opt_tracks, collapse=", "), "\n\n")
    return(FALSE)
  }
  
  # Compare dimensions
  all_match <- TRUE
  for(track in orig_tracks) {
    orig_dim <- dim(orig_result[[track]])
    opt_dim <- dim(opt_result[[track]])
    
    if(!identical(orig_dim, opt_dim)) {
      cat(sprintf("  ✗ Dimensions differ for track '%s'\n", track))
      cat("    Original:", paste(orig_dim, collapse=" x "), "\n")
      cat("    Optimized:", paste(opt_dim, collapse=" x "), "\n")
      all_match <- FALSE
    }
  }
  
  # Compare values (allowing for small numerical differences)
  for(track in orig_tracks) {
    max_diff <- max(abs(orig_result[[track]] - opt_result[[track]]), na.rm = TRUE)
    
    if(max_diff > 1) {  # Allow for rounding differences in INT16
      cat(sprintf("  ✗ Values differ significantly for track '%s' (max diff: %.2f)\n", 
                  track, max_diff))
      all_match <- FALSE
    }
  }
  
  if(all_match) {
    cat("  ✓ Results are consistent\n\n")
  } else {
    cat("\n")
  }
  
  return(all_match)
}

# Source the optimized implementations
source("R/python_ssff_optimized.R")

cat("=== SWIPE Tests ===\n\n")

if(exists("swipe") && exists("swipe_opt")) {
  # Test consistency
  swipe_consistent <- test_consistency(
    swipe, swipe_opt, "SWIPE",
    listOfFiles = test_file,
    windowShift = 5,
    minF = 70,
    maxF = 200
  )
  
  # Benchmark original
  cat("Original SWIPE:\n")
  swipe_orig_bench <- run_benchmark(
    "swipe (original)",
    swipe,
    listOfFiles = test_file,
    windowShift = 5,
    minF = 70,
    maxF = 200
  )
  
  # Benchmark optimized
  cat("Optimized SWIPE:\n")
  swipe_opt_bench <- run_benchmark(
    "swipe_opt",
    swipe_opt,
    listOfFiles = test_file,
    windowShift = 5,
    minF = 70,
    maxF = 200
  )
  
  if(!is.null(swipe_orig_bench) && !is.null(swipe_opt_bench)) {
    speedup <- swipe_orig_bench$mean / swipe_opt_bench$mean
    cat(sprintf("Speedup: %.2fx\n\n", speedup))
  }
} else {
  cat("SWIPE functions not available\n\n")
}

cat("=== RAPT Tests ===\n\n")

if(exists("rapt") && exists("rapt_opt")) {
  # Test consistency
  rapt_consistent <- test_consistency(
    rapt, rapt_opt, "RAPT",
    listOfFiles = test_file,
    windowShift = 5,
    minF = 70,
    maxF = 200
  )
  
  # Benchmark original
  cat("Original RAPT:\n")
  rapt_orig_bench <- run_benchmark(
    "rapt (original)",
    rapt,
    listOfFiles = test_file,
    windowShift = 5,
    minF = 70,
    maxF = 200
  )
  
  # Benchmark optimized
  cat("Optimized RAPT:\n")
  rapt_opt_bench <- run_benchmark(
    "rapt_opt",
    rapt_opt,
    listOfFiles = test_file,
    windowShift = 5,
    minF = 70,
    maxF = 200
  )
  
  if(!is.null(rapt_orig_bench) && !is.null(rapt_opt_bench)) {
    speedup <- rapt_orig_bench$mean / rapt_opt_bench$mean
    cat(sprintf("Speedup: %.2fx\n\n", speedup))
  }
} else {
  cat("RAPT functions not available\n\n")
}

cat("=== Batch Processing Test ===\n\n")

# Test with multiple files
test_files <- list.files("tests/signalfiles/AVQI/input/", 
                        pattern = "*.wav", full.names = TRUE)

if(length(test_files) > 1) {
  cat(sprintf("Testing batch processing with %d files\n\n", length(test_files)))
  
  if(exists("swipe_opt")) {
    cat("Batch SWIPE (optimized):\n")
    batch_time <- system.time({
      suppressMessages(swipe_opt(
        listOfFiles = test_files,
        toFile = FALSE,
        verbose = FALSE
      ))
    })
    
    cat(sprintf("  Total time: %.3f s\n", batch_time[3]))
    cat(sprintf("  Per file: %.3f s\n\n", batch_time[3] / length(test_files)))
  }
}

cat("=== Summary ===\n\n")
cat("Optimizations implemented:\n")
cat("  ✓ Batch file processing\n")
cat("  ✓ Efficient file format checking (Rcpp)\n")
cat("  ✓ Proper error handling with cli messages\n")
cat("  ✓ Media file conversion support\n")
cat("  ✓ Consistent interface with forest() template\n")
cat("  ✓ Proper logging and verbose output\n")
cat("  ✓ Time windowing support\n")
cat("\n")

cat("Recommendations:\n")
cat("  • Use reticulate for Python DSP (more stable than Rcpp Python)\n")
cat("  • Initialize Python environment once, reuse for all files\n")
cat("  • Use vectorized file processing where possible\n")
cat("  • Implement proper cleanup of temporary files\n")
cat("  • Add progress bars for large batch jobs\n")
