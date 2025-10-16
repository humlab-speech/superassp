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

# Benchmark 1: Formant analysis
formant_bench <- microbenchmark(
  "forest" = forest(test_file, toFile = FALSE, verbose = FALSE),
  times = 30,
  unit = "ms"
)

print(summary(formant_bench)[, c("expr", "min", "lq", "mean", "median", "uq", "max")])

# Create plot for formant benchmark
p1 <- autoplot(formant_bench) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Formant Analysis Performance (superassp::forest)",
    subtitle = sprintf("Processing %.2fs audio file with av-based media loading", duration),
    x = "Time (milliseconds)",
    y = NULL
  ) +
  theme(
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave("/tmp/benchmark_formant.png", p1, width = 8, height = 3, dpi = 150)
cat("\nSaved plot: /tmp/benchmark_formant.png\n")

cat("\n=== Running Pitch Tracking Benchmarks ===\n")

# Benchmark 2: Pitch tracking comparison - only test available functions
pitch_results <- data.frame(
  algorithm = character(),
  median_ms = numeric(),
  stringsAsFactors = FALSE
)

test_pitch_functions <- list(
  "ksvfo" = list(fn = fo, label = "KSV (autocorrelation)"),
  "mhspitch" = list(fn = pitch, label = "MHS (cepstrum)"),
  "rapt" = list(fn = rapt, label = "RAPT"),
  "swipe" = list(fn = swipe, label = "SWIPE"),
  "reaper" = list(fn = reaper, label = "REAPER")
)

cat("\nTesting pitch tracking algorithms:\n")
for (algo_name in names(test_pitch_functions)) {
  algo_info <- test_pitch_functions[[algo_name]]
  tryCatch({
    # Run quick benchmark
    bench <- microbenchmark(
      algo_info$fn(test_file, toFile = FALSE, verbose = FALSE),
      times = 20,
      unit = "ms"
    )
    median_time <- median(bench$time) / 1e6  # Convert to ms
    pitch_results <- rbind(pitch_results, data.frame(
      algorithm = algo_info$label,
      median_ms = median_time
    ))
    cat(sprintf("  ✓ %-25s: %6.2f ms (median)\n", algo_info$label, median_time))
  }, error = function(e) {
    cat(sprintf("  ✗ %-25s: %s\n", algo_info$label, conditionMessage(e)))
  })
}

if (nrow(pitch_results) > 0) {
  # Sort by performance
  pitch_results <- pitch_results[order(pitch_results$median_ms), ]
  pitch_results$algorithm <- factor(pitch_results$algorithm,
                                   levels = rev(pitch_results$algorithm))

  # Create plot
  p2 <- ggplot(pitch_results, aes(x = median_ms, y = algorithm)) +
    geom_col(fill = "steelblue", alpha = 0.8) +
    geom_text(aes(label = sprintf("%.1f ms", median_ms)),
             hjust = -0.1, size = 3.5) +
    theme_minimal(base_size = 12) +
    labs(
      title = "Pitch Tracking Algorithm Performance Comparison",
      subtitle = sprintf("Processing %.2fs audio file (median time from 20 runs)", duration),
      x = "Time (milliseconds)",
      y = NULL
    ) +
    theme(
      plot.title = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      panel.grid.major.y = element_blank(),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    expand_limits(x = max(pitch_results$median_ms) * 1.15)

  ggsave("/tmp/benchmark_pitch.png", p2, width = 9, height = 4, dpi = 150)
  cat("\nSaved plot: /tmp/benchmark_pitch.png\n")
}

cat("\n=== Running Parallel Processing Benchmark ===\n")

# Benchmark 3: Parallel vs Sequential processing
test_files <- rep(test_file, 20)  # Use more files for dramatic difference

cat(sprintf("Testing with %d files (total %.1f seconds of audio)...\n",
           length(test_files), duration * length(test_files)))

parallel_bench <- microbenchmark(
  "Sequential" = {
    lapply(test_files, function(f) rmsana(f, toFile = FALSE, verbose = FALSE))
  },
  "Parallel" = {
    rmsana(test_files, toFile = FALSE, verbose = FALSE)
  },
  times = 10,
  unit = "ms"
)

summary_stats <- summary(parallel_bench)
print(summary_stats[, c("expr", "min", "lq", "mean", "median", "uq", "max")])

# Calculate speedup
seq_median <- summary_stats$median[summary_stats$expr == "Sequential"]
par_median <- summary_stats$median[summary_stats$expr == "Parallel"]
speedup <- seq_median / par_median

# Detect number of cores
n_cores <- parallel::detectCores()

cat(sprintf("\nSpeedup: %.2fx (%.1f ms → %.1f ms) on %d cores\n",
           speedup, seq_median, par_median, n_cores - 1))

# Create plot for parallel benchmark
p3 <- autoplot(parallel_bench) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Parallel Processing Performance",
    subtitle = sprintf("Processing %d files: %.2fx speedup on %d cores",
                      length(test_files), speedup, n_cores - 1),
    x = "Time (milliseconds)",
    y = NULL
  ) +
  theme(
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave("/tmp/benchmark_parallel.png", p3, width = 8, height = 3.5, dpi = 150)
cat("\nSaved plot: /tmp/benchmark_parallel.png\n")

cat("\n=== All benchmarks complete ===\n")
cat("\nGenerated plots:\n")
cat("  • benchmark_formant.png\n")
cat("  • benchmark_pitch.png\n")
cat("  • benchmark_parallel.png\n")
