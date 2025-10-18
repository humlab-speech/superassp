##' ESTK Pitchmark - Find glottal closure instants in laryngograph signals
##'
##' @description This function finds instants of glottal closure in laryngograph
##'   (EGG/electroglottograph) waveforms using the Edinburgh Speech Tools pitchmark
##'   algorithm. It can also process regular speech WAV files, though it is optimized
##'   for laryngograph signals.
##'
##'   The algorithm performs the following operations:
##'   \enumerate{
##'     \item Double low-pass filter the signal to remove noise
##'     \item Double high-pass filter to remove low-frequency swell
##'     \item Calculate delta (differentiated) signal
##'     \item Low-pass filter the delta signal
##'     \item Pick negative zero crossings as pitchmarks
##'   }
##'
##'   Input files in formats not natively supported by ESTK will be loaded via the
##'   av package, allowing processing of any media format including video files.
##'
##'   The results will be written to an SSFF formatted file with the base name of
##'   the input file and extension *.pm* containing pitchmark times.
##'
##' @details This function wraps the Edinburgh Speech Tools pitchmark algorithm,
##'   which was designed for finding glottal closure instants in laryngograph
##'   waveforms. The algorithm uses carefully tuned filtering operations to isolate
##'   the pitch pulses. All input media formats are supported via av package integration.
##'
##' @note Pitchmarking is most accurate with lossless audio formats (WAV, FLAC, etc.).
##'   Lossy compression (MP3, OGG, AAC) may affect the accuracy of pitchmark detection.
##'
##' @param listOfFiles Vector of file paths to process. Can be audio files (WAV, EGG, FLAC, etc.)
##'   or video files (MP4, MKV, AVI, etc.). The av package will extract audio automatically.
##' @param beginTime Start time in seconds for analysis window (default: 0.0)
##' @param endTime End time in seconds for analysis window (default: 0.0 = end of file)
##' @param lx_low_frequency Low-pass cutoff frequency for initial filtering (default: 400 Hz).
##'   This removes high-frequency noise from the laryngograph signal.
##' @param lx_low_order Order of the low-pass FIR filter (default: 19).
##'   Higher values give sharper cutoff but more computation.
##' @param lx_high_frequency High-pass cutoff frequency for initial filtering (default: 40 Hz).
##'   This removes the low-frequency swell often seen in laryngograph signals.
##' @param lx_high_order Order of the high-pass FIR filter (default: 19).
##' @param df_low_frequency Low-pass cutoff for the differentiated signal (default: 1000 Hz).
##'   Applied after differentiation to smooth the signal.
##' @param df_low_order Order of the differentiated signal low-pass filter (default: 19).
##'   Set to 0 to disable this filtering stage.
##' @param median_order Order of median smoother for the differentiated signal (default: 19).
##'   Set to 0 to disable median smoothing.
##' @param fill Logical, whether to post-process pitchmarks (default: FALSE).
##'   If TRUE, ensures minimum/maximum pitch periods and fills unvoiced regions.
##' @param min_period Minimum allowed pitch period in seconds (default: 0.003 = ~333 Hz max F0).
##'   Used when fill=TRUE to remove spurious close pitchmarks.
##' @param max_period Maximum allowed pitch period in seconds (default: 0.02 = ~50 Hz min F0).
##'   Used when fill=TRUE to insert interpolated pitchmarks in unvoiced regions.
##' @param def_period Default pitch period for interpolation in seconds (default: 0.01 = 100 Hz).
##'   Used when fill=TRUE to determine spacing of interpolated pitchmarks.
##' @param invert Logical, invert the signal polarity (default: FALSE).
##'   Sometimes laryngograph signals are recorded upside-down; use this to correct.
##' @param to_f0 Logical, convert pitchmarks to F0 contour (default: FALSE).
##'   If TRUE, returns F0 values derived from pitch period instead of pitchmark times.
##' @param toFile Logical, write results to file (default: TRUE).
##'   If FALSE, returns results as R objects.
##' @param explicitExt Output file extension (default: "pm" for pitchmarks, "f0" if to_f0=TRUE)
##' @param outputDirectory Optional directory for output files (default: NULL = same as input)
##' @param verbose Logical, show progress messages (default: TRUE)
##' @param parallel Logical, use parallel processing for multiple files (default: NULL = auto).
##'   Automatically enabled for batches of 2+ files.
##' @param n_cores Number of CPU cores to use for parallel processing (default: NULL = auto).
##'   Defaults to detectCores() - 1.
##' @param use_cpp Logical, use C++ implementation (default: TRUE).
##'   If FALSE, falls back to calling ESTK binary (slower, requires temp files).
##'
##' @return If toFile=TRUE, returns the number of successfully processed files.
##'   If toFile=FALSE, returns a list of data frames with pitchmark times (and F0 values if to_f0=TRUE).
##'   For single file input, returns a single data frame instead of a list.
##'
##' @references
##' Edinburgh Speech Tools: \url{https://www.cstr.ed.ac.uk/projects/speech_tools/}
##'
##' Pitchmarking algorithm developed by Mike Macon and Paul Taylor,
##' Centre for Speech Technology Research, University of Edinburgh.
##'
##' @export
##' @examples
##' \dontrun{
##' # Basic pitchmarking of laryngograph file
##' estk_pitchmark("recording.egg")
##'
##' # Pitchmark with filling (interpolate unvoiced regions)
##' estk_pitchmark("recording.wav", fill = TRUE,
##'                min = 0.003, max = 0.02, def = 0.01)
##'
##' # Extract F0 from pitchmarks
##' f0_data <- estk_pitchmark("recording.egg", to_f0 = TRUE, toFile = FALSE)
##'
##' # Process video file (extracts audio automatically)
##' estk_pitchmark("interview.mp4")
##'
##' # Batch processing with parallel execution
##' files <- c("rec1.wav", "rec2.wav", "rec3.wav")
##' estk_pitchmark(files, parallel = TRUE, n_cores = 4)
##'
##' # Inverted laryngograph signal
##' estk_pitchmark("recording.egg", invert = TRUE)
##' }
##'
estk_pitchmark <- function(listOfFiles = NULL,
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
  funName <- "estk_pitchmark"

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
    cli::cli_inform("Applying {.fun {funName}} to {cli::no(n_files)} recording{?s}")
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
      audio_obj <- av_to_asspDataObj(
        file_path,
        start_time = bt,
        end_time = if (et == 0.0) NULL else et,
        target_sample_rate = NULL
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
attr(estk_pitchmark, "ext") <- "pm"
attr(estk_pitchmark, "tracks") <- "pitchmarks"
attr(estk_pitchmark, "outputType") <- "EST_Track"
attr(estk_pitchmark, "nativeFiletypes") <- c("wav", "egg")
attr(estk_pitchmark, "suggestCaching") <- FALSE
