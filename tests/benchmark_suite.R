#!/usr/bin/env Rscript
# Comprehensive Benchmark Suite for SuperASP DSP Functions
#
# This script benchmarks:
# 1. Parselmouth-optimized vs original Praat-calling functions
# 2. SuperASP vs wrassp implementations
#
# Results are saved as CSV and RDS for the Quarto report

library(superassp)
library(wrassp)
library(bench)
library(dplyr)
library(tidyr)

# Setup Python for Parselmouth functions
Sys.unsetenv("RETICULATE_PYTHON")
reticulate::use_python("/usr/bin/python3", required = FALSE)

# Check if parselmouth is available
has_parselmouth <- tryCatch({
  reticulate::py_module_available("parselmouth")
}, error = function(e) FALSE)

# Get test files
cat("Finding test audio files...\n")
test_files <- list.files(
  path = "signalfiles",
  pattern = "\\.wav$",
  recursive = TRUE,
  full.names = TRUE
)

# Select diverse test files
test_files_short <- test_files[grep("AVQI/input/(sv1|cs1)\\.wav$", test_files)]
test_files_medium <- test_files[grep("DSI/input/(fh1|im1)\\.wav$", test_files)]

if (length(test_files_short) == 0) {
  test_files_short <- test_files[1:min(2, length(test_files))]
}
if (length(test_files_medium) == 0) {
  test_files_medium <- test_files[1:min(2, length(test_files))]
}

all_test_files <- unique(c(test_files_short, test_files_medium))

cat(sprintf("Selected %d test files:\n", length(all_test_files)))
for (f in all_test_files) {
  cat(sprintf("  - %s\n", basename(f)))
}

# Initialize results list
benchmark_results <- list()

# =============================================================================
# Part 1: Parselmouth Optimized vs Original
# =============================================================================

if (has_parselmouth) {
  cat("\n=== Benchmarking Parselmouth-optimized functions ===\n")

  parselmouth_pairs <- list(
    list(
      name = "Formant Analysis (Burg)",
      orig = "praat_formant_burg",
      opt = "praat_formant_burg_opt",
      args = list(timeStep = 0.005, number_of_formants = 5)
    ),
    list(
      name = "Pitch Tracking",
      orig = "praat_pitch",
      opt = "praat_pitch_opt",
      args = list(time_step = 0.005)
    ),
    list(
      name = "Intensity Analysis",
      orig = "praat_intensity",
      opt = "praat_intensity_opt",
      args = list(time_step = 0.0)
    ),
    list(
      name = "Spectral Moments",
      orig = "praat_spectral_moments",
      opt = "praat_spectral_moments_opt",
      args = list(time_step = 0.005)
    ),
    list(
      name = "FormantPath (Burg)",
      orig = "praat_formantpath_burg",
      opt = "praat_formantpath_burg_opt",
      args = list(time_step = 0.005, track_formants = TRUE)
    )
  )

  for (pair in parselmouth_pairs) {
    # Check if functions exist
    if (!exists(pair$orig, mode = "function") ||
        !exists(pair$opt, mode = "function")) {
      cat(sprintf("  Skipping %s (function not found)\n", pair$name))
      next
    }

    cat(sprintf("\nBenchmarking: %s\n", pair$name))

    for (test_file in all_test_files) {
      cat(sprintf("  File: %s\n", basename(test_file)))

      # Prepare arguments
      args_orig <- c(list(listOfFiles = test_file, toFile = FALSE), pair$args)
      args_opt <- c(list(listOfFiles = test_file, toFile = FALSE), pair$args)

      # Run benchmark
      bm <- tryCatch({
        bench::mark(
          original = do.call(pair$orig, args_orig),
          optimized = do.call(pair$opt, args_opt),
          iterations = 10,
          check = FALSE,
          time_unit = "ms"
        )
      }, error = function(e) {
        cat(sprintf("    Error: %s\n", e$message))
        NULL
      })

      if (!is.null(bm)) {
        # Add metadata
        bm$function_name <- pair$name
        bm$file <- basename(test_file)
        bm$category <- "parselmouth"

        benchmark_results[[length(benchmark_results) + 1]] <- bm

        # Print summary
        speedup <- as.numeric(bm$median[1]) / as.numeric(bm$median[2])
        cat(sprintf("    Original: %.2f ms, Optimized: %.2f ms, Speedup: %.2fx\n",
                    as.numeric(bm$median[1]),
                    as.numeric(bm$median[2]),
                    speedup))
      }
    }
  }
} else {
  cat("\nParselmouth not available - skipping optimized function benchmarks\n")
}

# =============================================================================
# Part 2: SuperASP vs wrassp
# =============================================================================

cat("\n\n=== Benchmarking SuperASP vs wrassp ===\n")

common_functions <- list(
  list(name = "ACF Analysis", fn = "acfana", args = list()),
  list(name = "RFC Analysis", fn = "rfcana", args = list()),
  list(name = "RMS Analysis", fn = "rmsana", args = list()),
  list(name = "ZCR Analysis", fn = "zcrana", args = list())
)

for (fn_info in common_functions) {
  cat(sprintf("\nBenchmarking: %s\n", fn_info$name))

  for (test_file in all_test_files) {
    cat(sprintf("  File: %s\n", basename(test_file)))

    args_list <- c(list(listOfFiles = test_file, toFile = FALSE), fn_info$args)

    # Run benchmark
    bm <- tryCatch({
      bench::mark(
        superassp = do.call(get(fn_info$fn, envir = asNamespace("superassp")),
                           args_list),
        wrassp = do.call(get(fn_info$fn, envir = asNamespace("wrassp")),
                        args_list),
        iterations = 10,
        check = FALSE,
        time_unit = "ms"
      )
    }, error = function(e) {
      cat(sprintf("    Error: %s\n", e$message))
      NULL
    })

    if (!is.null(bm)) {
      # Add metadata
      bm$function_name <- fn_info$name
      bm$file <- basename(test_file)
      bm$category <- "superassp_vs_wrassp"

      benchmark_results[[length(benchmark_results) + 1]] <- bm

      # Print summary
      ratio <- as.numeric(bm$median[1]) / as.numeric(bm$median[2])
      comparison <- if (ratio < 1) {
        sprintf("SuperASP %.2fx faster", 1/ratio)
      } else {
        sprintf("wrassp %.2fx faster", ratio)
      }
      cat(sprintf("    SuperASP: %.2f ms, wrassp: %.2f ms (%s)\n",
                  as.numeric(bm$median[1]),
                  as.numeric(bm$median[2]),
                  comparison))
    }
  }
}

# =============================================================================
# Save Results
# =============================================================================

if (length(benchmark_results) > 0) {
  cat("\n\nSaving results...\n")

  # Combine all benchmarks
  all_benchmarks <- bind_rows(benchmark_results)

  # Save as RDS for detailed analysis
  saveRDS(all_benchmarks, "benchmark_results.rds")
  cat("  Saved: benchmark_results.rds\n")

  # Create summary CSV
  summary_df <- all_benchmarks %>%
    select(category, function_name, file, expression, median, mem_alloc) %>%
    mutate(
      median_ms = as.numeric(median) * 1000,
      mem_mb = as.numeric(mem_alloc) / 1024^2
    ) %>%
    select(-median, -mem_alloc)

  write.csv(summary_df, "benchmark_summary.csv", row.names = FALSE)
  cat("  Saved: benchmark_summary.csv\n")

  # Print overall summary
  cat("\n=== Overall Summary ===\n")

  # Parselmouth speedups
  if (has_parselmouth && any(all_benchmarks$category == "parselmouth")) {
    parselmouth_summary <- all_benchmarks %>%
      filter(category == "parselmouth") %>%
      group_by(function_name, file) %>%
      summarise(
        speedup = as.numeric(median[expression == "original"]) /
                  as.numeric(median[expression == "optimized"]),
        .groups = "drop"
      )

    cat("\nParselmouth Optimization Speedups:\n")
    print(parselmouth_summary, n = Inf)

    avg_speedup <- mean(parselmouth_summary$speedup, na.rm = TRUE)
    cat(sprintf("\nAverage speedup: %.2fx\n", avg_speedup))
  }

  # SuperASP vs wrassp comparison
  if (any(all_benchmarks$category == "superassp_vs_wrassp")) {
    comparison_summary <- all_benchmarks %>%
      filter(category == "superassp_vs_wrassp") %>%
      group_by(function_name, file) %>%
      summarise(
        superassp_ms = as.numeric(median[expression == "superassp"]),
        wrassp_ms = as.numeric(median[expression == "wrassp"]),
        ratio = superassp_ms / wrassp_ms,
        .groups = "drop"
      )

    cat("\nSuperASP vs wrassp Performance:\n")
    print(comparison_summary, n = Inf)
  }

  cat("\n✓ Benchmark complete!\n")
} else {
  cat("\nNo benchmark results generated.\n")
}
