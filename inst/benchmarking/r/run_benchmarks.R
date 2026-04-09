#!/usr/bin/env Rscript
# Benchmark script for superassp README
# This script runs real benchmarks and generates plots

library(superassp)
library(microbenchmark)
library(ggplot2)

# Get test file from package installation (use larger file for better timing)
test_file <- system.file("samples", "sustained", "a32b.wav", package = "superassp")
if (!file.exists(test_file)) {
  # Fallback to development path
  test_file <- file.path(getwd(), "inst/samples/sustained/a32b.wav")
}

cat("Using test file:", test_file, "\n")
if (!file.exists(test_file)) {
  stop("Test file not found!")
}

# Get file duration for context
info <- av::av_media_info(test_file)
duration <- info$duration
cat(sprintf("File duration: %.2f seconds\n\n", duration))

cat("=== Running Formant Benchmarks ===\n")

# Benchmark 1: Formant analysis - compare multiple methods
formant_results <- data.frame(
  method = character(),
  median_ms = numeric(),
  stringsAsFactors = FALSE
)

# Benchmark all formant methods together for violin plot
cat("\nTesting formant analysis methods:\n")

formant_exprs <- list()

# Add wrassp::forest
tryCatch({
  test_result <- wrassp::forest(test_file, toFile = FALSE)
  formant_exprs[["wrassp::forest"]] <- quote(wrassp::forest(test_file, toFile = FALSE))
  cat("  ✓ wrassp::forest available\n")
}, error = function(e) {
  cat(sprintf("  ✗ wrassp::forest: %s\n", conditionMessage(e)))
})


# Add superassp::trk_forest
tryCatch({
  test_result <- trk_forest(test_file, toFile = FALSE, verbose = FALSE)
  formant_exprs[["superassp::trk_forest"]] <- quote(trk_forest(test_file, toFile = FALSE, verbose = FALSE))
  cat("  ✓ superassp::trk_forest available\n")
}, error = function(e) {
  cat(sprintf("  ✗ superassp::trk_forest: %s\n", conditionMessage(e)))
})

# Add trk_formant (Praat Burg method)
tryCatch({
  test_result <- trk_formant(test_file, toFile = FALSE)
  formant_exprs[["trk_formant"]] <- quote(trk_formant(test_file, toFile = FALSE))
  cat("  ✓ trk_formant available\n")
}, error = function(e) {
  cat(sprintf("  ✗ trk_formant: %s\n", conditionMessage(e)))
})

# Add trk_praat_sauce (VoiceSauce-style)
tryCatch({
  test_result <- trk_praat_sauce(test_file, toFile = FALSE)
  formant_exprs[["trk_praat_sauce"]] <- quote(trk_praat_sauce(test_file, toFile = FALSE))
  cat("  ✓ trk_praat_sauce available\n")
}, error = function(e) {
  cat(sprintf("  ✗ trk_praat_sauce: %s (skipping)\n", conditionMessage(e)))
})


# Add trk_deepformants (Deep learning) - skip if PyTorch not available
tryCatch({
  test_result <- trk_deepformants(test_file, toFile = FALSE, verbose = FALSE)
  formant_exprs[["trk_deepformants"]] <- quote(trk_deepformants(test_file, toFile = FALSE, verbose = FALSE))
  cat("  ✓ trk_deepformants available\n")
}, error = function(e) {
  cat(sprintf("  ✗ trk_deepformants: %s (skipping - PyTorch required)\n", conditionMessage(e)))


})

if (length(formant_exprs) > 0) {
  cat(sprintf("\nRunning formant benchmarks (100 iterations for %d methods)...\n",
              length(formant_exprs)))

  formant_bench <- microbenchmark(
    list = formant_exprs,
    times = 100,
    unit = "ms"
  )

  # Print summary
  cat("\nFormant Analysis Summary:\n")
  print(summary(formant_bench)[, c("expr", "min", "lq", "mean", "median", "uq", "max")])

  # Create violin plot
  p1 <- autoplot(formant_bench) +
    theme_minimal(base_size = 12) +
    labs(
      title = "Formant Analysis Performance Comparison",
      subtitle = sprintf("Processing %.2fs audio file (100 runs per method)", duration),
      x = "Time (milliseconds)",
      y = NULL
    ) +
    theme(
      plot.title = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      plot.margin = margin(10, 10, 10, 10)
    )

  ggsave("/tmp/benchmark_formant.png", p1, width = 10, height = 5, dpi = 150)
  cat("\nSaved plot: /tmp/benchmark_formant.png\n")
} else {
  cat("\nNo formant analysis functions available for benchmarking\n")
}

cat("\n=== Running Pitch Tracking Benchmarks ===\n")

# Benchmark 2: Pitch tracking - test all available algorithms
cat("\nTesting pitch tracking algorithms:\n")

# Load audio object once for C++ in-memory functions
audio_obj <- NULL
tryCatch({
  audio_obj <- av_to_asspDataObj(test_file)
  cat("  ✓ Audio loaded for C++ in-memory functions\n")
}, error = function(e) {
  cat(sprintf("  ✗ Failed to load audio: %s\n", conditionMessage(e)))
})

pitch_exprs <- list()


# Test KSV - ASSP fo function

tryCatch({
  test_result <- fo(test_file, toFile = FALSE, verbose = FALSE)
  pitch_exprs[["KSV (autocorrelation)"]] <- quote(fo(test_file, toFile = FALSE, verbose = FALSE))
  cat("  ✓ KSV (autocorrelation) available\n")
}, error = function(e) {
  cat(sprintf("  ✗ KSV: %s\n", conditionMessage(e)))
})



# Test MHS - ASSP pitch function

tryCatch({
  test_result <- pitch(test_file, toFile = FALSE, verbose = FALSE)
  pitch_exprs[["MHS (cepstrum)"]] <- quote(pitch(test_file, toFile = FALSE, verbose = FALSE))
  cat("  ✓ MHS (cepstrum) available\n")
}, error = function(e) {
  cat(sprintf("  ✗ MHS: %s\n", conditionMessage(e)))
})

# Test SPTK C++ functions (if audio loaded)


# Note: RAPT and DIO may fail with "Failed to initialize" errors on some systems
if (!is.null(audio_obj)) {
  # Test RAPT C++ - SKIP DUE TO KNOWN SEGFAULT ISSUE
  # RAPT works for single calls but crashes when called repeatedly in a loop
  # due to static variables in SPTK's Snack implementation (upstream bug)
  cat("  ⊘ RAPT C++ - skipped (segfault when called repeatedly - SPTK library bug)\n")


  # Test SWIPE C++
  tryCatch({
    test_result <- swipe_cpp(audio_obj, minF = 60, maxF = 400, windowShift = 10, verbose = FALSE)
    pitch_exprs[["SWIPE C++"]] <- quote(swipe_cpp(audio_obj, minF = 60, maxF = 400, windowShift = 10, verbose = FALSE))
    cat("  ✓ SWIPE C++ available\n")
  }, error = function(e) {
    cat(sprintf("  ✗ SWIPE C++: %s\n", conditionMessage(e)))
  })

  # Test REAPER C++
  tryCatch({
    test_result <- reaper_cpp(audio_obj, minF = 60, maxF = 400, windowShift = 10, verbose = FALSE)
    pitch_exprs[["REAPER C++"]] <- quote(reaper_cpp(audio_obj, minF = 60, maxF = 400, windowShift = 10, verbose = FALSE))
    cat("  ✓ REAPER C++ available\n")
  }, error = function(e) {
    cat(sprintf("  ✗ REAPER C++: %s\n", conditionMessage(e)))
  })

  # Test DIO C++
  tryCatch({
    test_result <- dio_cpp(audio_obj, minF = 60, maxF = 400, windowShift = 10, verbose = FALSE)
    pitch_exprs[["DIO C++ (WORLD)"]] <- quote(dio_cpp(audio_obj, minF = 60, maxF = 400, windowShift = 10, verbose = FALSE))
    cat("  ✓ DIO C++ available\n")
  }, error = function(e) {

    cat(sprintf("  ✗ DIO C++: %s (skipping - initialization issue)\n", conditionMessage(e)))
  })
}


# Test trk_straight_f0 (STRAIGHT baseline - important for accuracy comparison)
# Note: Known to have segfault issues, skip for now until fixed
tryCatch({
  # Skip STRAIGHT for now due to segfault issues
  cat("  ⊘ trk_straight_f0 (STRAIGHT baseline) - skipped (segfault issue)\n")
}, error = function(e) {
  cat(sprintf("  ✗ trk_straight_f0: %s\n", conditionMessage(e)))


})

if (length(pitch_exprs) > 0) {
  cat(sprintf("\nRunning pitch tracking benchmarks (100 iterations for %d algorithms)...\n",
              length(pitch_exprs)))

  # Suppress REAPER output messages during benchmarking
  pitch_bench <- suppressMessages({
    microbenchmark(
      list = pitch_exprs,
      times = 100,
      unit = "ms"
    )
  })

  # Print summary
  cat("\nPitch Tracking Summary:\n")
  print(summary(pitch_bench)[, c("expr", "min", "lq", "mean", "median", "uq", "max")])

  # Create violin plot
  p2 <- autoplot(pitch_bench) +
    theme_minimal(base_size = 12) +
    labs(
      title = "Pitch Tracking Algorithm Performance Comparison",
      subtitle = sprintf("Processing %.2fs audio file (100 runs per algorithm)", duration),
      x = "Time (milliseconds)",
      y = NULL
    ) +
    theme(
      plot.title = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      plot.margin = margin(10, 10, 10, 10)
    )

  ggsave("/tmp/benchmark_pitch.png", p2, width = 10, height = 5, dpi = 150)
  cat("\nSaved plot: /tmp/benchmark_pitch.png\n")
}

cat("\n=== Running Parallel Processing Benchmark ===\n")

# Benchmark 3: Parallel vs Sequential processing
test_files <- rep(test_file, 20)  # Use more files for dramatic difference

cat(sprintf("Testing with %d files (total %.1f seconds of audio)...\n",
           length(test_files), duration * length(test_files)))

cat("Running parallel processing benchmarks (100 iterations)...\n")

parallel_bench <- microbenchmark(
  "Sequential" = {

    lapply(test_files, function(f) trk_rmsana(f, toFile = FALSE, verbose = FALSE))
  },
  "Parallel" = {
    trk_rmsana(test_files, toFile = FALSE, verbose = FALSE)

  },
  times = 100,
  unit = "ms"
)

summary_stats <- summary(parallel_bench)
cat("\nParallel Processing Summary:\n")
print(summary_stats[, c("expr", "min", "lq", "mean", "median", "uq", "max")])

# Calculate speedup
seq_median <- summary_stats$median[summary_stats$expr == "Sequential"]
par_median <- summary_stats$median[summary_stats$expr == "Parallel"]
speedup <- seq_median / par_median

# Detect number of cores
n_cores <- parallel::detectCores()

cat(sprintf("\nSpeedup: %.2fx (%.1f ms → %.1f ms) on %d cores\n",
           speedup, seq_median, par_median, n_cores - 1))

# Create violin plot for parallel benchmark
p3 <- autoplot(parallel_bench) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Parallel Processing Performance",
    subtitle = sprintf("Processing %d files: %.2fx speedup on %d cores (100 runs)",
                      length(test_files), speedup, n_cores - 1),
    x = "Time (milliseconds)",
    y = NULL
  ) +
  theme(
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave("/tmp/benchmark_parallel.png", p3, width = 10, height = 4, dpi = 150)
cat("\nSaved plot: /tmp/benchmark_parallel.png\n")

cat("\n=== All benchmarks complete ===\n")
cat("\nGenerated plots:\n")
cat("  • benchmark_formant.png\n")
cat("  • benchmark_pitch.png\n")
cat("  • benchmark_parallel.png\n")
