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
      orig = "trk_formant",
      opt = "praat_formant_burg_opt",
      args = list(timeStep = 0.005, number_of_formants = 5)
    ),
    list(
      name = "Pitch Tracking",
      orig = "trk_pitch_cc",
      opt = "praat_pitch_opt",
      args = list(time_step = 0.005)
    ),
    list(
      name = "Intensity Analysis",
      orig = "trk_intensity",
      opt = "praat_intensity_opt",
      args = list(time_step = 0.0)
    ),
    list(
      name = "Spectral Moments",
      orig = "trk_spectral_moments",
      opt = "praat_spectral_moments_opt",
      args = list(time_step = 0.005)
    ),
    list(
      name = "FormantPath (Burg)",
      orig = "trk_formantathp",
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
  list(name = "ACF Analysis", fn = "trk_acfana", args = list()),
  list(name = "RFC Analysis", fn = "rfcana", args = list()),
  list(name = "RMS Analysis", fn = "trk_rmsana", args = list()),
  list(name = "ZCR Analysis", fn = "trk_zcrana", args = list())
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
# Part 3: Comprehensive Pitch Tracking Comparison
# =============================================================================

cat("\n\n=== Comprehensive Pitch Tracking Benchmarks ===\n")

# Check availability of pitch tracking methods
has_rapt <- exists("trk_rapt", mode = "function")
has_swipe <- exists("trk_swipe", mode = "function")
has_reaper <- exists("trk_reaper", mode = "function")
has_dio <- exists("trk_dio", mode = "function")
has_rapt_cpp <- exists("rapt_cpp", mode = "function")
has_swipe_cpp <- exists("swipe_cpp", mode = "function")
has_reaper_cpp <- exists("reaper_cpp", mode = "function")
has_dio_cpp <- exists("dio_cpp", mode = "function")
has_estk_pda <- exists("estk_pda_cpp", mode = "function")
has_praat_pitch <- exists("praat_pitch_opt", mode = "function") || exists("trk_pitch_cc", mode = "function")

# Check Python environment for SPTK methods
has_pysptk <- tryCatch({
  reticulate::py_module_available("pysptk")
}, error = function(e) FALSE)

cat("\nAvailable pitch tracking methods:\n")
cat(sprintf("  - KSV F0 (ASSP): %s\n", "✓"))
cat(sprintf("  - MHS Pitch (ASSP): %s\n", "✓"))
cat(sprintf("  - ESTK PDA (C++): %s\n", ifelse(has_estk_pda, "✓", "✗")))
cat(sprintf("  - RAPT C++ (SPTK): %s\n", ifelse(has_rapt, "✓", "✗")))
cat(sprintf("  - SWIPE C++ (SPTK): %s\n", ifelse(has_swipe, "✓", "✗")))
cat(sprintf("  - REAPER C++ (SPTK): %s\n", ifelse(has_reaper, "✓", "✗")))
cat(sprintf("  - DIO C++ (SPTK): %s\n", ifelse(has_dio, "✓", "✗")))
cat(sprintf("  - Praat Pitch: %s\n", ifelse(has_praat_pitch && has_parselmouth, "✓", "✗")))

# Run comprehensive pitch tracking benchmark
if (length(all_test_files) > 0) {
  cat("\nBenchmarking pitch tracking algorithms...\n")

  for (test_file in all_test_files) {
    cat(sprintf("\n  File: %s\n", basename(test_file)))

    # Build list of expressions to benchmark
    pitch_expressions <- list(
      ksv_f0 = quote(trk_ksvfo(test_file, toFile = FALSE, minF = 60, maxF = 400,
                           windowShift = 10.0, verbose = FALSE)),
      mhs_pitch = quote(trk_mhspitch(test_file, toFile = FALSE, minF = 60, maxF = 400,
                                  windowShift = 10.0, verbose = FALSE))
    )

    # Add C++ wrapper methods (preferred - fast and full-featured)
    if (has_rapt) {
      pitch_expressions$rapt_cpp_wrapper <- quote(trk_rapt(test_file, toFile = FALSE, minF = 60,
                                                       maxF = 400, windowShift = 10.0, verbose = FALSE))
    }

    if (has_swipe) {
      pitch_expressions$swipe_cpp_wrapper <- quote(trk_swipe(test_file, toFile = FALSE, minF = 60,
                                                         maxF = 400, windowShift = 10.0, verbose = FALSE))
    }

    if (has_reaper) {
      pitch_expressions$reaper_cpp_wrapper <- quote(trk_reaper(test_file, toFile = FALSE, minF = 60,
                                                           maxF = 400, windowShift = 10.0, verbose = FALSE))
    }

    if (has_dio) {
      pitch_expressions$dio_cpp_wrapper <- quote(trk_dio(test_file, toFile = FALSE, minF = 60,
                                                     maxF = 400, windowShift = 10.0, verbose = FALSE))
    }

    # Add ESTK PDA if available
    if (has_estk_pda) {
      # ESTK PDA needs AsspDataObj
      pitch_expressions$estk_pda <- quote({
        audio_obj <- av_to_asspDataObj(test_file)
        estk_pda_cpp(audio_obj, minF = 60, maxF = 400, windowShift = 10.0)
      })
    }

    # Add Praat pitch if available
    if (has_praat_pitch && has_parselmouth) {
      if (exists("praat_pitch_opt", mode = "function")) {
        pitch_expressions$praat_pitch <- quote(praat_pitch_opt(test_file, toFile = FALSE,
                                                                 time_step = 0.01,
                                                                 pitch_floor = 60,
                                                                 pitch_ceiling = 400))
      } else {
        pitch_expressions$praat_pitch <- quote(trk_pitch_cc(test_file, toFile = FALSE,
                                                            time_step = 0.01,
                                                            pitch_floor = 60,
                                                            pitch_ceiling = 400))
      }
    }

    # Run benchmark
    bm <- tryCatch({
      do.call(bench::mark, c(pitch_expressions,
                             list(iterations = 5,  # Fewer iterations for slower Python methods
                                  check = FALSE,
                                  time_unit = "ms")))
    }, error = function(e) {
      cat(sprintf("    Error: %s\n", e$message))
      NULL
    })

    if (!is.null(bm)) {
      # Add metadata
      bm$function_name <- "Pitch Tracking (All Methods)"
      bm$file <- basename(test_file)
      bm$category <- "pitch_tracking"

      benchmark_results[[length(benchmark_results) + 1]] <- bm

      # Print summary
      for (i in 1:nrow(bm)) {
        method_name <- as.character(bm$expression[i])
        cat(sprintf("    %s: %.2f ms\n", method_name, as.numeric(bm$median[i])))
      }

      # Calculate speedups relative to fastest
      fastest_time <- min(as.numeric(bm$median))
      fastest_idx <- which.min(as.numeric(bm$median))
      fastest_method <- as.character(bm$expression[fastest_idx])

      cat(sprintf("\n    Fastest: %s (%.2f ms)\n", fastest_method, fastest_time))

      # Show speedups for other methods
      for (i in 1:nrow(bm)) {
        if (i != fastest_idx) {
          method_name <- as.character(bm$expression[i])
          speedup <- as.numeric(bm$median[i]) / fastest_time
          cat(sprintf("    %s is %.2fx slower\n", method_name, speedup))
        }
      }
    }
  }
} else {
  cat("\nNo test files available for pitch tracking benchmarks\n")
}

# =============================================================================
# Part 4: ESTK Algorithms
# =============================================================================

cat("\n\n=== Benchmarking ESTK Algorithms ===\n")

# Check if ESTK functions are available
has_estk_pda <- exists("estk_pda_cpp", mode = "function")
has_estk_pitchmark <- exists("estk_pitchmark_cpp", mode = "function")

if (has_estk_pda || has_estk_pitchmark) {
  # ESTK needs AsspDataObj, so we'll convert files first
  cat("\nLoading test files as AsspDataObj...\n")

  estk_test_objs <- lapply(all_test_files, function(f) {
    tryCatch({
      av_to_asspDataObj(f)
    }, error = function(e) {
      cat(sprintf("  Failed to load %s: %s\n", basename(f), e$message))
      NULL
    })
  })
  names(estk_test_objs) <- basename(all_test_files)
  estk_test_objs <- estk_test_objs[!sapply(estk_test_objs, is.null)]

  # ESTK PDA benchmarks
  if (has_estk_pda) {
    cat("\nBenchmarking: ESTK PDA (Super Resolution Pitch Detection)\n")

    # Compare with other pitch tracking methods
    pitch_methods <- list(
      list(name = "ESTK PDA", fn = function(obj) {
        estk_pda_cpp(obj, minF = 60, maxF = 400, windowShift = 10.0)
      }),
      list(name = "MHS Pitch", fn = function(file) {
        trk_mhspitch(file, toFile = FALSE, minF = 60, maxF = 400, windowShift = 10.0)
      }),
      list(name = "KSV F0", fn = function(file) {
        trk_ksvfo(file, toFile = FALSE, minF = 60, maxF = 400, windowShift = 10.0)
      })
    )

    for (i in seq_along(estk_test_objs)) {
      test_file <- all_test_files[i]
      test_obj <- estk_test_objs[[i]]

      cat(sprintf("  File: %s\n", basename(test_file)))

      # Run benchmark
      bm <- tryCatch({
        bench::mark(
          estk_pda = estk_pda_cpp(test_obj, minF = 60, maxF = 400, windowShift = 10.0),
          mhs_pitch = trk_mhspitch(test_file, toFile = FALSE, minF = 60, maxF = 400,
                               windowShift = 10.0),
          ksv_f0 = trk_ksvfo(test_file, toFile = FALSE, minF = 60, maxF = 400,
                         windowShift = 10.0),
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
        bm$function_name <- "Pitch Detection (ESTK PDA vs others)"
        bm$file <- basename(test_file)
        bm$category <- "estk_pitch"

        benchmark_results[[length(benchmark_results) + 1]] <- bm

        # Print summary
        cat(sprintf("    ESTK PDA: %.2f ms\n", as.numeric(bm$median[1])))
        cat(sprintf("    MHS Pitch: %.2f ms\n", as.numeric(bm$median[2])))
        cat(sprintf("    KSV F0: %.2f ms\n", as.numeric(bm$median[3])))

        # Calculate speedups
        estk_vs_mhs <- as.numeric(bm$median[2]) / as.numeric(bm$median[1])
        estk_vs_ksv <- as.numeric(bm$median[3]) / as.numeric(bm$median[1])

        if (estk_vs_mhs > 1) {
          cat(sprintf("    ESTK PDA %.2fx faster than MHS\n", estk_vs_mhs))
        } else {
          cat(sprintf("    MHS %.2fx faster than ESTK PDA\n", 1/estk_vs_mhs))
        }

        if (estk_vs_ksv > 1) {
          cat(sprintf("    ESTK PDA %.2fx faster than KSV\n", estk_vs_ksv))
        } else {
          cat(sprintf("    KSV %.2fx faster than ESTK PDA\n", 1/estk_vs_ksv))
        }
      }
    }
  }

  # ESTK Pitchmark benchmarks
  if (has_estk_pitchmark) {
    cat("\n\nBenchmarking: ESTK Pitchmark (Glottal Closure Detection)\n")

    for (i in seq_along(estk_test_objs)) {
      test_file <- all_test_files[i]
      test_obj <- estk_test_objs[[i]]

      cat(sprintf("  File: %s\n", basename(test_file)))

      # Benchmark different parameter settings
      bm <- tryCatch({
        bench::mark(
          default = estk_pitchmark_cpp(test_obj),
          with_fill = estk_pitchmark_cpp(test_obj, fill = TRUE,
                                         min_period = 0.003, max_period = 0.02),
          with_f0 = estk_pitchmark_cpp(test_obj, to_f0 = TRUE),
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
        bm$function_name <- "Pitchmark Detection (ESTK)"
        bm$file <- basename(test_file)
        bm$category <- "trk_pitchmark"

        benchmark_results[[length(benchmark_results) + 1]] <- bm

        # Print summary
        cat(sprintf("    Default: %.2f ms\n", as.numeric(bm$median[1])))
        cat(sprintf("    With fill: %.2f ms\n", as.numeric(bm$median[2])))
        cat(sprintf("    With F0: %.2f ms\n", as.numeric(bm$median[3])))
      }
    }
  }

} else {
  cat("\nESTK algorithms not available - skipping ESTK benchmarks\n")
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

  # Comprehensive pitch tracking comparison
  if (any(all_benchmarks$category == "pitch_tracking")) {
    pitch_summary <- all_benchmarks %>%
      filter(category == "pitch_tracking") %>%
      group_by(file) %>%
      summarise(
        n_methods = n(),
        fastest_method = as.character(expression[which.min(median)]),
        fastest_time_ms = min(as.numeric(median)),
        slowest_method = as.character(expression[which.max(median)]),
        slowest_time_ms = max(as.numeric(median)),
        speed_range = max(as.numeric(median)) / min(as.numeric(median)),
        .groups = "drop"
      )

    cat("\nPitch Tracking Methods Comparison:\n")
    print(pitch_summary, n = Inf)

    # Overall average by method across all files
    method_averages <- all_benchmarks %>%
      filter(category == "pitch_tracking") %>%
      group_by(expression) %>%
      summarise(
        avg_time_ms = mean(as.numeric(median)),
        .groups = "drop"
      ) %>%
      arrange(avg_time_ms)

    cat("\nAverage Performance by Method (across all files):\n")
    print(method_averages, n = Inf)
  }

  # ESTK pitch detection comparison
  if (any(all_benchmarks$category == "estk_pitch")) {
    estk_pitch_summary <- all_benchmarks %>%
      filter(category == "estk_pitch") %>%
      group_by(function_name, file) %>%
      summarise(
        estk_pda_ms = as.numeric(median[expression == "estk_pda"]),
        mhs_ms = as.numeric(median[expression == "mhs_pitch"]),
        ksv_ms = as.numeric(median[expression == "ksv_f0"]),
        vs_mhs = mhs_ms / estk_pda_ms,
        vs_ksv = ksv_ms / estk_pda_ms,
        .groups = "drop"
      )

    cat("\nESTK PDA vs Other Pitch Detection Methods:\n")
    print(estk_pitch_summary, n = Inf)

    avg_vs_mhs <- mean(estk_pitch_summary$vs_mhs, na.rm = TRUE)
    avg_vs_ksv <- mean(estk_pitch_summary$vs_ksv, na.rm = TRUE)

    cat(sprintf("\nAverage performance vs MHS: %.2fx\n", avg_vs_mhs))
    cat(sprintf("Average performance vs KSV: %.2fx\n", avg_vs_ksv))
  }

  # ESTK pitchmark performance
  if (any(all_benchmarks$category == "trk_pitchmark")) {
    estk_pm_summary <- all_benchmarks %>%
      filter(category == "trk_pitchmark") %>%
      group_by(function_name, file) %>%
      summarise(
        default_ms = as.numeric(median[expression == "default"]),
        with_fill_ms = as.numeric(median[expression == "with_fill"]),
        with_f0_ms = as.numeric(median[expression == "with_f0"]),
        .groups = "drop"
      )

    cat("\nESTK Pitchmark Performance:\n")
    print(estk_pm_summary, n = Inf)

    cat(sprintf("\nAverage pitchmark detection time: %.2f ms\n",
                mean(estk_pm_summary$default_ms, na.rm = TRUE)))
  }

  cat("\n✓ Benchmark complete!\n")
} else {
  cat("\nNo benchmark results generated.\n")
}

# =============================================================================
# Part 4: Voice Quality Analysis Benchmarks (lst_ functions)
# =============================================================================

cat("\n\n=== Voice Quality Analysis Benchmarks ===\n")

# Check availability of voice quality analysis functions
has_covarep_vq <- exists("lst_covarep_vq", mode = "function") &&
  tryCatch(covarep_available(), error = function(e) FALSE)

has_voice_sauce <- exists("lst_voice_sauce", mode = "function") &&
  tryCatch(voice_sauce_available(), error = function(e) FALSE)

cat("\nAvailable voice quality analysis methods:\n")
cat(sprintf("  - COVAREP VQ (8 params): %s\n", ifelse(has_covarep_vq, "✓", "✗")))
cat(sprintf("  - VoiceSauce (40+ params): %s\n", ifelse(has_voice_sauce, "✓", "✗")))

if (has_covarep_vq || has_voice_sauce) {
  cat("\nBenchmarking voice quality analysis...\n")

  # Use sustained vowel files (ideal for voice quality)
  vq_test_files <- test_files[grep("sustained|vowel|sv1", test_files, ignore.case = TRUE)]

  if (length(vq_test_files) == 0) {
    vq_test_files <- all_test_files[1:min(2, length(all_test_files))]
  }

  for (test_file in vq_test_files) {
    cat(sprintf("\n  File: %s\n", basename(test_file)))

    vq_expressions <- list()

    # Add COVAREP VQ if available
    if (has_covarep_vq) {
      vq_expressions$covarep_vq_basic <- quote(
        lst_covarep_vq(test_file, verbose = FALSE)
      )

      vq_expressions$covarep_vq_with_f0 <- quote(
        lst_covarep_vq(test_file, f0 = 150, verbose = FALSE)
      )
    }

    # Add VoiceSauce if available
    if (has_voice_sauce) {
      vq_expressions$voice_sauce_reaper <- quote(
        lst_voice_sauce(test_file, f0_method = "reaper",
                       frame_shift = 1.0, verbose = FALSE)
      )

      vq_expressions$voice_sauce_praat <- quote(
        lst_voice_sauce(test_file, f0_method = "praat",
                       frame_shift = 1.0, verbose = FALSE)
      )

      vq_expressions$voice_sauce_coarse <- quote(
        lst_voice_sauce(test_file, f0_method = "reaper",
                       frame_shift = 10.0, verbose = FALSE)
      )
    }

    if (length(vq_expressions) > 0) {
      # Run benchmark
      bm <- bench::mark(
        !!!vq_expressions,
        iterations = 5,
        check = FALSE
      )

      # Store results
      for (i in seq_len(nrow(bm))) {
        benchmark_results[[length(benchmark_results) + 1]] <- list(
          category = "voice_quality",
          function_name = names(vq_expressions)[i],
          file = basename(test_file),
          median_ms = as.numeric(bm$median[i]) * 1000,
          expression = as.character(bm$expression[i])
        )
      }

      # Print summary
      cat("\n    Results:\n")
      for (i in seq_len(nrow(bm))) {
        cat(sprintf("      %s: %.2f ms\n",
                    names(vq_expressions)[i],
                    as.numeric(bm$median[i]) * 1000))
      }
    }
  }

  # Summary of voice quality benchmarks
  if (length(benchmark_results) > 0) {
    vq_benchmarks <- do.call(rbind, lapply(benchmark_results, function(x) {
      if (x$category == "voice_quality") {
        data.frame(
          function_name = x$function_name,
          file = x$file,
          median_ms = x$median_ms,
          stringsAsFactors = FALSE
        )
      } else {
        NULL
      }
    }))

    if (!is.null(vq_benchmarks) && nrow(vq_benchmarks) > 0) {
      cat("\n\nVoice Quality Analysis Performance Summary:\n")

      vq_summary <- vq_benchmarks %>%
        group_by(function_name) %>%
        summarise(
          avg_ms = mean(median_ms, na.rm = TRUE),
          min_ms = min(median_ms, na.rm = TRUE),
          max_ms = max(median_ms, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        arrange(avg_ms)

      print(vq_summary, n = Inf)

      if (has_covarep_vq && has_voice_sauce) {
        cat("\n\nComparison (COVAREP vs VoiceSauce):\n")
        covarep_avg <- mean(vq_benchmarks$median_ms[grep("covarep", vq_benchmarks$function_name)])
        vs_avg <- mean(vq_benchmarks$median_ms[grep("voice_sauce", vq_benchmarks$function_name)])
        cat(sprintf("  COVAREP average: %.2f ms (8 parameters)\n", covarep_avg))
        cat(sprintf("  VoiceSauce average: %.2f ms (40+ parameters)\n", vs_avg))
        cat(sprintf("  VoiceSauce overhead: %.2fx for %.1fx more parameters\n",
                    vs_avg / covarep_avg,
                    40 / 8))
      }
    }
  }
} else {
  cat("\nNo voice quality analysis functions available.\n")
  cat("Install with:\n")
  cat("  install_covarep()\n")
  cat("  install_voice_sauce()\n")
}
