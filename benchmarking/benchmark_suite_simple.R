#!/usr/bin/env Rscript
# Simplified Benchmark Suite for SuperASP DSP Functions
#
# Benchmarks SuperASP vs wrassp implementations

library(superassp)
library(wrassp)
library(bench)
library(dplyr)

cat("SuperASP vs wrassp Benchmark Suite\n")
cat("===================================\n\n")

# Get test files
test_files <- list.files(
  path = "signalfiles",
  pattern = "\\.wav$",
  recursive = TRUE,
  full.names = TRUE
)

# Select a representative subset
test_files <- test_files[grep("(sv1|cs1|fh1)\\.wav$", test_files)]

if (length(test_files) == 0) {
  test_files <- test_files[1:min(3, length(test_files))]
}

cat(sprintf("Testing with %d audio files\n\n", length(test_files)))

# Common functions to benchmark
functions_to_test <- list(
  list(name = "ACF Analysis", fn = "trk_acfana"),
  list(name = "RFC Analysis", fn = "rfcana"),
  list(name = "RMS Analysis", fn = "trk_rmsana"),
  list(name = "ZCR Analysis", fn = "trk_zcrana")
)

results <- list()

for (fn_info in functions_to_test) {
  cat(sprintf("Benchmarking: %s\n", fn_info$name))
  cat(strrep("-", 50), "\n")

  for (test_file in test_files) {
    cat(sprintf("  %s... ", basename(test_file)))

    bm <- tryCatch({
      bench::mark(
        superassp = do.call(
          get(fn_info$fn, envir = asNamespace("superassp")),
          list(listOfFiles = test_file, toFile = FALSE)
        ),
        wrassp = do.call(
          get(fn_info$fn, envir = asNamespace("wrassp")),
          list(listOfFiles = test_file, toFile = FALSE)
        ),
        iterations = 5,
        check = FALSE,
        time_unit = "ms"
      )
    }, error = function(e) {
      cat(sprintf("ERROR: %s\n", e$message))
      NULL
    })

    if (!is.null(bm)) {
      superassp_time <- as.numeric(bm$median[1])
      wrassp_time <- as.numeric(bm$median[2])
      ratio <- superassp_time / wrassp_time

      if (ratio < 1) {
        msg <- sprintf("SuperASP %.2fx faster", 1/ratio)
      } else {
        msg <- sprintf("wrassp %.2fx faster", ratio)
      }

      cat(sprintf("%.2f ms vs %.2f ms (%s)\n",
                  superassp_time, wrassp_time, msg))

      # Store results
      bm$function_name <- fn_info$name
      bm$file <- basename(test_file)
      results[[length(results) + 1]] <- bm
    }
  }
  cat("\n")
}

# Combine and save results
if (length(results) > 0) {
  all_results <- bind_rows(results)

  # Save results
  saveRDS(all_results, "benchmark_results.rds")
  cat("\n✓ Results saved to benchmark_results.rds\n")

  # Print summary
  summary <- all_results %>%
    group_by(function_name) %>%
    summarise(
      superassp_median = median(as.numeric(median[expression == "superassp"])),
      wrassp_median = median(as.numeric(median[expression == "wrassp"])),
      ratio = superassp_median / wrassp_median,
      .groups = "drop"
    )

  cat("\n")
  cat("Summary Statistics\n")
  cat(strrep("=", 50), "\n")
  print(summary, n = Inf)

  cat("\n✓ Benchmark complete!\n")
} else {
  cat("\n✗ No benchmarks completed\n")
}
