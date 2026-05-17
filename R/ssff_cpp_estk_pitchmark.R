##' Detect glottal closure instants in laryngograph signals using ESTk pitchmark
##'
##' Finds pitchmark times in laryngograph (EGG) or speech waveforms using the
##' Edinburgh Speech Tools pitchmark algorithm (zero-crossing detection on the
##' filtered, differentiated signal). Prefer \code{protoscribe::draft_pitchmark()}
##' for new code; this function is retained for backwards compatibility only.
##'
##' @note **DEPRECATED**. Use \code{protoscribe::draft_pitchmark()} instead, which
##'   follows the \code{draft_} naming convention and integrates with reindeer workflows.
##'
##' @note Pitchmarking accuracy is highest with lossless formats (WAV, FLAC). Lossy
##'   formats (MP3, AAC) may degrade detection.
##'
##' @param listOfFiles Character vector of audio or EGG file paths. Any format supported
##'   by \pkg{av} is accepted; non-native inputs are transcoded automatically.
##' @param beginTime Numeric. Start of analysis window in seconds. Default 0 (file start).
##' @param endTime Numeric. End of analysis window in seconds. Default 0 (file end).
##' @param lx_low_frequency Numeric. Low-pass cutoff in Hz for initial denoising filter.
##'   Default 400 Hz.
##' @param lx_low_order Integer. Order of the initial low-pass FIR filter. Default 19.
##' @param lx_high_frequency Numeric. High-pass cutoff in Hz to remove low-frequency
##'   swell. Default 40 Hz.
##' @param lx_high_order Integer. Order of the high-pass FIR filter. Default 19.
##' @param df_low_frequency Numeric. Low-pass cutoff in Hz applied to the differentiated
##'   signal. Default 1000 Hz.
##' @param df_low_order Integer. Order of the differentiated-signal low-pass filter.
##'   Set to 0 to disable. Default 19.
##' @param median_order Integer. Order of the median smoother on the differentiated signal.
##'   Set to 0 to disable. Default 19.
##' @param fill Logical. If \code{TRUE}, post-process pitchmarks: remove marks closer than
##'   \code{min_period}, interpolate gaps larger than \code{max_period}. Default \code{FALSE}.
##' @param min_period Numeric. Minimum pitch period in seconds (used when \code{fill = TRUE}).
##'   Default 0.003 s (≈333 Hz max F0).
##' @param max_period Numeric. Maximum pitch period in seconds (used when \code{fill = TRUE}).
##'   Default 0.02 s (≈50 Hz min F0).
##' @param def_period Numeric. Default pitch period for interpolated marks (used when
##'   \code{fill = TRUE}). Default 0.01 s (100 Hz).
##' @param invert Logical. Invert signal polarity before processing (use for upside-down
##'   EGG recordings). Default \code{FALSE}.
##' @param to_f0 Logical. If \code{TRUE}, return F0 values derived from pitchmark intervals
##'   instead of raw pitchmark times. Default \code{FALSE}.
##' @param toFile Logical. If \code{TRUE}, write output files and return the count written
##'   invisibly. If \code{FALSE}, return results as R objects. Default \code{TRUE}.
##' @param explicitExt Character. Output file extension. Default \code{"pm"} (or
##'   \code{"f0"} when \code{to_f0 = TRUE}).
##' @param outputDirectory Character. Directory for output files. \code{NULL} (default)
##'   writes alongside the input file.
##' @param verbose Logical. Print per-file progress. Default \code{TRUE}.
##' @param parallel Logical. Use parallel processing for multiple files. \code{NULL}
##'   (default) enables automatically for 2+ files.
##' @param n_cores Integer. Number of cores for parallel processing.
##'   \code{NULL} (default) uses \code{detectCores() - 1}.
##' @param use_cpp Logical. Use C++ implementation (default \code{TRUE}). Setting
##'   \code{FALSE} falls back to the ESTK binary (slower, requires temp files).
##'
##' @return If \code{toFile = TRUE}: integer count of files written, returned invisibly.
##'   If \code{toFile = FALSE} and \code{to_f0 = FALSE}: a data frame (single file) or
##'   list of data frames with pitchmark times in seconds.
##'   If \code{toFile = FALSE} and \code{to_f0 = TRUE}: a data frame (or list) with F0
##'   values derived from inter-pitchmark intervals.
##'
##' @references
##' \insertCite{EdinburghSpeechTools2020}{superassp}
##'
##' \insertCite{Macon1997Pitchmark}{superassp}
##'
##' @export
##' @examples
##' \dontrun{
##' # Basic pitchmarking of laryngograph file
##' trk_pitchmark_estk("recording.egg")
##'
##' # Pitchmark with filling (interpolate unvoiced regions)
##' trk_pitchmark_estk("recording.wav", fill = TRUE,
##'                min = 0.003, max = 0.02, def = 0.01)
##'
##' # Extract F0 from pitchmarks
##' f0_data <- trk_pitchmark_estk("recording.egg", to_f0 = TRUE, toFile = FALSE)
##'
##' # Process video file (extracts audio automatically)
##' trk_pitchmark_estk("interview.mp4")
##'
##' # Batch processing with parallel execution
##' files <- c("rec1.wav", "rec2.wav", "rec3.wav")
##' trk_pitchmark_estk(files, parallel = TRUE, n_cores = 4)
##'
##' # Inverted laryngograph signal
##' trk_pitchmark_estk("recording.egg", invert = TRUE)
##' }
##'
trk_pitchmark_estk <- function(listOfFiles,
                           beginTime = 0.0,
                           endTime = 0.0,
                           lx_low_frequency = 400,
                           lx_low_order = 19,
                           lx_high_frequency = 40,
                           lx_high_order = 19,
                           df_low_frequency = 1000,
                           df_low_order = 19,
                           median_order = 19,
                           fill = FALSE,
                           min_period = 0.003,
                           max_period = 0.02,
                           def_period = 0.01,
                           invert = FALSE,
                           to_f0 = FALSE,
                           toFile = TRUE,
                           explicitExt = NULL,
                           outputDirectory = NULL,
                           verbose = TRUE,
                           parallel = NULL,
                           n_cores = NULL,
                           use_cpp = TRUE) {

  # Deprecation warning
  if (verbose) {
    cli::cli_alert_warning(c(
      "{.fn trk_pitchmark_estk} is deprecated",
      "i" = "Use {.fn protoscribe::draft_pitchmark} instead",
      "i" = "New function follows draft_ pattern and integrates with reindeer"
    ))
  }

  # Initial validation and setup
  if (is.null(listOfFiles) || length(listOfFiles) == 0) {
    cli::cli_abort("No input files specified in {.arg listOfFiles}")
  }

  # Determine output extension
  if (is.null(explicitExt)) {
    explicitExt <- if (to_f0) "f0" else "pm"
  }

  # Get function metadata
  currCall <- rlang::current_call()
  funName <- "trk_pitchmark_estk"

  # Supported formats - ESTK pitchmark works with WAV files
  # We'll use av to convert other formats to WAV
  nativeFiletypes <- c("wav", "egg")  # EGG files are typically WAV format

  # Normalize paths
  listOfFiles <- fast_strip_file_protocol(listOfFiles)
  listOfFiles <- normalizePath(path.expand(listOfFiles), mustWork = FALSE)

  # Check file existence
  files_exist <- file.exists(listOfFiles)
  if (!all(files_exist)) {
    missing_files <- listOfFiles[!files_exist]
    cli::cli_abort(c(
      "!" = "Some files do not exist:",
      "x" = "{.file {fast_basename(missing_files)}}"
    ))
  }

  n_files <- length(listOfFiles)

  # Normalize time parameters
  beginTime <- if (is.null(beginTime)) 0.0 else beginTime
  endTime <- if (is.null(endTime)) 0.0 else endTime

  # Recycle time parameters to match file count
  if (length(beginTime) == 1) beginTime <- rep(beginTime, n_files)
  if (length(endTime) == 1) endTime <- rep(endTime, n_files)

  # Validate time parameters
  if (length(beginTime) != n_files || length(endTime) != n_files) {
    cli::cli_abort("Time parameters must have length 1 or match number of files")
  }

  # Warn about lossless formats
  known_lossless <- knownLossless()
  file_exts <- fast_file_ext(listOfFiles)
  is_lossless <- tolower(file_exts) %in% tolower(known_lossless)
  not_lossless <- listOfFiles[!is_lossless]

  if (length(not_lossless) > 0 && verbose) {
    cli::cli_warn(c(
      "!" = "Found {length(not_lossless)} recording{?s} in lossy format{?s}",
      "i" = "Lossy compression may affect pitchmark detection accuracy",
      "x" = "For accurate results, use lossless formats: {.val {intersect(nativeFiletypes, known_lossless)}}"
    ))
  }

  # Setup output directory
  makeOutputDirectory(outputDirectory, FALSE, funName)

  # Auto-enable parallel for batches
  if (is.null(parallel)) {
    parallel <- n_files > 1
  }

  # Determine number of cores
  if (is.null(n_cores)) {
    n_cores <- parallel::detectCores() - 1
    if (is.na(n_cores) || n_cores < 1) n_cores <- 1
  }

  # Disable parallel for single file
  use_parallel <- parallel && n_files > 1 && n_cores > 1

  if (verbose) {
    format_apply_msg(funName, n_files, beginTime, endTime)
    if (use_cpp) {
      cli::cli_inform("Using fast C++ implementation (in-memory processing)")
    } else {
      cli::cli_inform("Using ESTK binary (requires temporary files)")
    }
    if (use_parallel) {
      cli::cli_inform("Using parallel processing on {n_cores} core{?s}")
    }
  }

  # Check for ESTK pitchmark binary (only if not using C++)
  estk_binary <- NULL
  if (!use_cpp) {
    estk_binary <- system.file("ESTK", "bin", "pitchmark", package = "superassp")
    if (estk_binary == "" || !file.exists(estk_binary)) {
      # Try alternative location (development)
      estk_binary <- file.path(getwd(), "src", "ESTK", "bin", "pitchmark")
      if (!file.exists(estk_binary)) {
        estk_binary <- "/Users/frkkan96/Documents/src/superassp/src/ESTK/bin/pitchmark"
      }
    }

    if (!file.exists(estk_binary)) {
      cli::cli_abort(c(
        "x" = "ESTK pitchmark binary not found",
        "i" = "Expected location: {.path {system.file('ESTK', 'bin', package = 'superassp')}}"
      ))
    }

    # Make binary executable
    Sys.chmod(estk_binary, mode = "0755", use_umask = TRUE)
  }

  # Define processing function for a single file
  process_single_file <- function(i) {
    file_path <- listOfFiles[i]
    bt <- beginTime[i]
    et <- endTime[i]

    tryCatch({
      # Load audio with av (handles all formats and time windowing)
      audio_obj <- read_audio(
        file_path,
        begin = bt,
        end   = et
      )

      if (use_cpp) {
        # ===== C++ Implementation (in-memory) =====
        # Call C++ pitchmark function directly
        pm_result <- estk_pitchmark_cpp(
          audio_obj = audio_obj,
          lx_low_frequency = as.integer(lx_low_frequency),
          lx_low_order = as.integer(lx_low_order),
          lx_high_frequency = as.integer(lx_high_frequency),
          lx_high_order = as.integer(lx_high_order),
          df_low_frequency = as.integer(df_low_frequency),
          df_low_order = as.integer(df_low_order),
          median_order = as.integer(median_order),
          fill = fill,
          min_period = min_period,
          max_period = max_period,
          def_period = def_period,
          invert = invert,
          to_f0 = to_f0,
          verbose = FALSE
        )

        # Handle output
        if (toFile) {
          # Prepare output file path
          out_dir <- if (is.null(outputDirectory)) {
            dirname(file_path)
          } else {
            outputDirectory
          }
          base_name <- fast_file_path_sans_ext(c(basename(file_path)))[1]
          output_file <- file.path(out_dir, paste0(base_name, ".", explicitExt))

          # Create AsspDataObj for writing
          # Structure depends on whether we're writing pitchmarks or F0
          if (to_f0 && !is.null(pm_result$f0)) {
            # Write F0 track
            out_obj <- list(
              f0 = pm_result$f0
            )
            attr(out_obj, "sampleRate") <- pm_result$sample_rate
            attr(out_obj, "startTime") <- 0.0
            attr(out_obj, "startRecord") <- 1L
            attr(out_obj, "endRecord") <- nrow(pm_result$f0)
            class(out_obj) <- c("AsspDataObj", "list")

            # Write to file
            write.AsspDataObj(out_obj, output_file)
          } else {
            # Write pitchmark track
            # Pitchmarks are event times, create single-column matrix
            pm_matrix <- matrix(pm_result$pitchmarks, ncol = 1)
            out_obj <- list(
              pm = pm_matrix
            )
            attr(out_obj, "sampleRate") <- pm_result$sample_rate
            attr(out_obj, "startTime") <- 0.0
            attr(out_obj, "startRecord") <- 1L
            attr(out_obj, "endRecord") <- nrow(pm_matrix)
            class(out_obj) <- c("AsspDataObj", "list")

            # Write to file
            write.AsspDataObj(out_obj, output_file)
          }

          return(TRUE)
        } else {
          # Return result as R object
          return(pm_result)
        }

      } else {
        # ===== Binary Implementation (requires temp files) =====
        # Write temporary WAV file for ESTK pitchmark
        temp_dir <- tempdir()
        temp_wav <- file.path(temp_dir, paste0("estk_pm_", i, "_", basename(file_path), ".wav"))

        # Convert AsspDataObj to WAV file
        write.AsspDataObj(audio_obj, temp_wav)

        # Prepare output file
        if (toFile) {
          out_dir <- if (is.null(outputDirectory)) {
            dirname(file_path)
          } else {
            outputDirectory
          }
          base_name <- fast_file_path_sans_ext(c(basename(file_path)))[1]
          output_file <- file.path(out_dir, paste0(base_name, ".", explicitExt))
        } else {
          output_file <- file.path(temp_dir, paste0("estk_pm_", i, "_output.", explicitExt))
        }

        # Build pitchmark command
        cmd_args <- c(
          temp_wav,
          "-o", output_file,
          "-otype", "est"
        )

        # Add filter parameters
        cmd_args <- c(cmd_args,
                      "-lx_lf", as.character(lx_low_frequency),
                      "-lx_lo", as.character(lx_low_order),
                      "-lx_hf", as.character(lx_high_frequency),
                      "-lx_ho", as.character(lx_high_order))

        if (df_low_order > 0) {
          cmd_args <- c(cmd_args,
                        "-df_lf", as.character(df_low_frequency),
                        "-df_lo", as.character(df_low_order))
        }

        if (median_order > 0) {
          cmd_args <- c(cmd_args, "-mean_o", as.character(median_order))
        }

        # Add filling parameters
        if (fill) {
          cmd_args <- c(cmd_args,
                        "-fill",
                        "-min", as.character(min_period),
                        "-max", as.character(max_period),
                        "-def", as.character(def_period),
                        "-wave_end")
        }

        # Add inversion if requested
        if (invert) {
          cmd_args <- c(cmd_args, "-inv")
        }

        # Add F0 conversion if requested
        if (to_f0) {
          f0_file <- sub("\\.pm$", ".f0", output_file)
          f0_file <- sub("\\.f0$", ".f0", f0_file)  # ensure .f0 extension
          cmd_args <- c(cmd_args, "-f0", f0_file)
        }

        # Execute ESTK pitchmark
        result <- system2(estk_binary, args = cmd_args, stdout = TRUE, stderr = TRUE)

        # Clean up temporary WAV
        unlink(temp_wav)

        # Read results if not writing to file
        if (!toFile) {
          if (to_f0 && file.exists(f0_file)) {
            # Read F0 track
            pm_data <- wrassp::read.AsspDataObj(f0_file)
            unlink(f0_file)
          } else if (file.exists(output_file)) {
            # Read pitchmark track
            pm_data <- wrassp::read.AsspDataObj(output_file)
            unlink(output_file)
          } else {
            cli::cli_warn("Failed to generate pitchmark output for {.file {basename(file_path)}}")
            return(NULL)
          }

          return(pm_data)
        } else {
          # Verify output file was created
          if (file.exists(output_file)) {
            return(TRUE)
          } else {
            cli::cli_warn("Failed to write pitchmark output for {.file {basename(file_path)}}")
            return(FALSE)
          }
        }
      }

    }, error = function(e) {
      cli::cli_warn("Error processing {.file {basename(file_path)}}: {conditionMessage(e)}")
      return(if (toFile) FALSE else NULL)
    })
  }

  # Process files (parallel or sequential)
  if (use_parallel) {
    if (.Platform$OS.type == "windows") {
      # Windows: socket cluster
      cl <- parallel::makeCluster(n_cores)
      on.exit(parallel::stopCluster(cl), add = TRUE)

      parallel::clusterExport(cl, c(
        "listOfFiles", "beginTime", "endTime", "lx_low_frequency", "lx_low_order",
        "lx_high_frequency", "lx_high_order", "df_low_frequency", "df_low_order",
        "median_order", "fill", "min_period", "max_period", "def_period",
        "invert", "to_f0", "toFile", "explicitExt", "outputDirectory",
        "use_cpp", "estk_binary", "process_single_file", "av_to_asspDataObj"
      ), envir = environment())

      parallel::clusterEvalQ(cl, {
        library(superassp)
      })

      if (verbose) {
        results <- pbapply::pblapply(seq_along(listOfFiles), process_single_file, cl = cl)
      } else {
        results <- parallel::parLapply(cl, seq_along(listOfFiles), process_single_file)
      }
    } else {
      # Unix/Mac: fork-based
      if (verbose && requireNamespace("pbmcapply", quietly = TRUE)) {
        results <- pbmcapply::pbmclapply(
          seq_along(listOfFiles),
          process_single_file,
          mc.cores = n_cores,
          mc.preschedule = TRUE
        )
      } else {
        results <- parallel::mclapply(
          seq_along(listOfFiles),
          process_single_file,
          mc.cores = n_cores,
          mc.preschedule = TRUE
        )
      }
    }
  } else {
    # Sequential processing
    results <- vector("list", n_files)
    if (verbose && n_files > 1) {
      cli::cli_progress_bar(
        "Processing files",
        total = n_files,
        format = "{cli::pb_spin} {cli::pb_current}/{cli::pb_total} | ETA: {cli::pb_eta}"
      )
      for (i in seq_along(listOfFiles)) {
        results[[i]] <- process_single_file(i)
        cli::cli_progress_update()
      }
      cli::cli_progress_done()
    } else {
      for (i in seq_along(listOfFiles)) {
        results[[i]] <- process_single_file(i)
      }
    }
  }

  # Process results
  if (toFile) {
    n_success <- sum(unlist(results), na.rm = TRUE)
    if (verbose) {
      cli::cli_inform("Successfully processed {n_success} of {n_files} file{?s}")
    }
    return(invisible(n_success))
  } else {
    # Remove NULL results
    results <- results[!sapply(results, is.null)]

    # Simplify for single file
    if (length(results) == 1) {
      return(results[[1]])
    } else {
      return(results)
    }
  }
}

# Function attributes
attr(trk_pitchmark_estk, "ext") <- "pm"
attr(trk_pitchmark_estk, "tracks") <- "pitchmarks"
attr(trk_pitchmark_estk, "outputType") <- "EST_Track"
attr(trk_pitchmark_estk, "nativeFiletypes") <- c("wav", "egg")
attr(trk_pitchmark_estk, "suggestCaching") <- FALSE
