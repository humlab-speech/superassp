##' RAPT Pitch Tracking (C++ implementation)
##'
##' @description Extract F0 (fundamental frequency) using the Robust Algorithm
##'   for Pitch Tracking (RAPT) from SPTK. This is a high-performance C++
##'   implementation that is 2-3x faster than the Python version and requires
##'   no Python dependencies.
##'
##'   The algorithm uses dynamic programming to find the optimal F0 contour,
##'   making it robust to noise and reliable for speech analysis.
##'
##'   All input media formats are supported via the av package, including video
##'   files from which audio will be automatically extracted.
##'
##' @param listOfFiles Vector of file paths to process
##' @param beginTime Start time in seconds (default: 0.0)
##' @param endTime End time in seconds (default: 0.0 = end of file)
##' @param windowShift Frame shift in milliseconds (default: 10.0)
##' @param minF Minimum F0 in Hz (default: 60.0)
##' @param maxF Maximum F0 in Hz (default: 400.0)
##' @param voicing_threshold Voicing threshold (default: 0.9)
##' @param toFile Write results to file (default: TRUE)
##' @param explicitExt Output file extension (default: "f0")
##' @param outputDirectory Output directory (default: NULL = same as input)
##' @param verbose Show progress messages (default: TRUE)
##'
##' @return If toFile=TRUE, returns the number of successfully processed files.
##'   If toFile=FALSE, returns AsspDataObj or list of AsspDataObj objects.
##'
##' @export
##' @examples
##' \dontrun{
##' # Extract F0 from audio file
##' rapt("recording.wav")
##'
##' # Process with custom F0 range
##' rapt("speech.mp3", minF = 75, maxF = 300)
##'
##' # Return data without writing file
##' f0_data <- rapt("audio.wav", toFile = FALSE)
##'
##' # Process video file (extracts audio)
##' rapt("interview.mp4")
##' }
rapt <- function(listOfFiles = NULL,
                 beginTime = 0.0,
                 endTime = 0.0,
                 windowShift = 10.0,
                 minF = 60.0,
                 maxF = 400.0,
                 voicing_threshold = 0.9,
                 toFile = TRUE,
                 explicitExt = "f0",
                 outputDirectory = NULL,
                 verbose = TRUE) {

  # Validate inputs
  if (is.null(listOfFiles) || length(listOfFiles) == 0) {
    cli::cli_abort("No input files specified in {.arg listOfFiles}")
  }

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

  # Recycle time parameters
  if (length(beginTime) == 1) beginTime <- rep(beginTime, n_files)
  if (length(endTime) == 1) endTime <- rep(endTime, n_files)

  # Setup output directory
  makeOutputDirectory(outputDirectory, FALSE, "rapt")

  if (verbose) {
    cli::cli_inform("Applying {.fun rapt} (C++) to {cli::no(n_files)} recording{?s}")
  }

  # Process each file
  results <- vector("list", n_files)

  if (verbose && n_files > 1) {
    cli::cli_progress_bar(
      "Processing files",
      total = n_files,
      format = "{cli::pb_spin} {cli::pb_current}/{cli::pb_total} | ETA: {cli::pb_eta}"
    )
  }

  for (i in seq_len(n_files)) {
    file_path <- listOfFiles[i]
    bt <- beginTime[i]
    et <- endTime[i]

    tryCatch({
      # Load audio with av
      audio_obj <- av_to_asspDataObj(
        file_path,
        start_time = bt,
        end_time = if (et == 0.0) NULL else et
      )

      # Call C++ RAPT
      rapt_result <- rapt_cpp(
        audio_obj = audio_obj,
        minF = minF,
        maxF = maxF,
        windowShift = windowShift,
        voicing_threshold = voicing_threshold,
        verbose = FALSE
      )

      # Convert to AsspDataObj
      out_obj <- create_f0_asspobj(rapt_result, windowShift)

      # Handle output
      if (toFile) {
        out_file <- generate_output_path(file_path, explicitExt, outputDirectory)
        write.AsspDataObj(out_obj, out_file)
        results[[i]] <- TRUE
      } else {
        results[[i]] <- out_obj
      }

    }, error = function(e) {
      cli::cli_warn("Error processing {.file {basename(file_path)}}: {conditionMessage(e)}")
      results[[i]] <- if (toFile) FALSE else NULL
    })

    if (verbose && n_files > 1) {
      cli::cli_progress_update()
    }
  }

  if (verbose && n_files > 1) {
    cli::cli_progress_done()
  }

  # Return results
  if (toFile) {
    n_success <- sum(unlist(results), na.rm = TRUE)
    if (verbose) {
      cli::cli_inform("Successfully processed {n_success} of {n_files} file{?s}")
    }
    return(invisible(n_success))
  } else {
    results <- results[!sapply(results, is.null)]
    if (length(results) == 1) {
      return(results[[1]])
    } else {
      return(results)
    }
  }
}

attr(rapt, "ext") <- "f0"
attr(rapt, "tracks") <- c("f0")
attr(rapt, "outputType") <- "SSFF"
attr(rapt, "nativeFiletypes") <- c("wav", "flac", "mp3", "mp4", "mkv", "avi")  # Via av
attr(rapt, "suggestCaching") <- FALSE


##' SWIPE Pitch Tracking (C++ implementation)
##'
##' @description Extract F0 using the Sawtooth Waveform Inspired Pitch Estimator
##'   (SWIPE) from SPTK. This is a high-performance C++ implementation that is
##'   2-3x faster than the Python version.
##'
##'   SWIPE uses spectral pattern matching and is particularly effective for noisy
##'   speech or challenging recording conditions.
##'
##' @inheritParams rapt
##' @param voicing_threshold Voicing threshold (default: 0.3, lower than RAPT)
##'
##' @return If toFile=TRUE, returns the number of successfully processed files.
##'   If toFile=FALSE, returns AsspDataObj or list of AsspDataObj objects.
##'
##' @export
##' @examples
##' \dontrun{
##' # Extract F0 from audio file
##' swipe("recording.wav")
##'
##' # Process with custom parameters
##' swipe("speech.wav", minF = 100, maxF = 500, voicing_threshold = 0.4)
##' }
swipe <- function(listOfFiles = NULL,
                  beginTime = 0.0,
                  endTime = 0.0,
                  windowShift = 10.0,
                  minF = 60.0,
                  maxF = 400.0,
                  voicing_threshold = 0.3,
                  toFile = TRUE,
                  explicitExt = "f0",
                  outputDirectory = NULL,
                  verbose = TRUE) {

  if (is.null(listOfFiles) || length(listOfFiles) == 0) {
    cli::cli_abort("No input files specified in {.arg listOfFiles}")
  }

  listOfFiles <- fast_strip_file_protocol(listOfFiles)
  listOfFiles <- normalizePath(path.expand(listOfFiles), mustWork = FALSE)

  files_exist <- file.exists(listOfFiles)
  if (!all(files_exist)) {
    missing_files <- listOfFiles[!files_exist]
    cli::cli_abort(c(
      "!" = "Some files do not exist:",
      "x" = "{.file {fast_basename(missing_files)}}"
    ))
  }

  n_files <- length(listOfFiles)

  beginTime <- if (is.null(beginTime)) 0.0 else beginTime
  endTime <- if (is.null(endTime)) 0.0 else endTime

  if (length(beginTime) == 1) beginTime <- rep(beginTime, n_files)
  if (length(endTime) == 1) endTime <- rep(endTime, n_files)

  makeOutputDirectory(outputDirectory, FALSE, "swipe")

  if (verbose) {
    cli::cli_inform("Applying {.fun swipe} (C++) to {cli::no(n_files)} recording{?s}")
  }

  results <- vector("list", n_files)

  if (verbose && n_files > 1) {
    cli::cli_progress_bar(
      "Processing files",
      total = n_files,
      format = "{cli::pb_spin} {cli::pb_current}/{cli::pb_total} | ETA: {cli::pb_eta}"
    )
  }

  for (i in seq_len(n_files)) {
    file_path <- listOfFiles[i]
    bt <- beginTime[i]
    et <- endTime[i]

    tryCatch({
      audio_obj <- av_to_asspDataObj(
        file_path,
        start_time = bt,
        end_time = if (et == 0.0) NULL else et
      )

      swipe_result <- swipe_cpp(
        audio_obj = audio_obj,
        minF = minF,
        maxF = maxF,
        windowShift = windowShift,
        voicing_threshold = voicing_threshold,
        verbose = FALSE
      )

      out_obj <- create_f0_asspobj(swipe_result, windowShift)

      if (toFile) {
        out_file <- generate_output_path(file_path, explicitExt, outputDirectory)
        write.AsspDataObj(out_obj, out_file)
        results[[i]] <- TRUE
      } else {
        results[[i]] <- out_obj
      }

    }, error = function(e) {
      cli::cli_warn("Error processing {.file {basename(file_path)}}: {conditionMessage(e)}")
      results[[i]] <- if (toFile) FALSE else NULL
    })

    if (verbose && n_files > 1) {
      cli::cli_progress_update()
    }
  }

  if (verbose && n_files > 1) {
    cli::cli_progress_done()
  }

  if (toFile) {
    n_success <- sum(unlist(results), na.rm = TRUE)
    if (verbose) {
      cli::cli_inform("Successfully processed {n_success} of {n_files} file{?s}")
    }
    return(invisible(n_success))
  } else {
    results <- results[!sapply(results, is.null)]
    if (length(results) == 1) {
      return(results[[1]])
    } else {
      return(results)
    }
  }
}

attr(swipe, "ext") <- "f0"
attr(swipe, "tracks") <- c("f0")
attr(swipe, "outputType") <- "SSFF"
attr(swipe, "nativeFiletypes") <- c("wav", "flac", "mp3", "mp4", "mkv", "avi")
attr(swipe, "suggestCaching") <- FALSE


##' REAPER Pitch and Epoch Tracking (C++ implementation)
##'
##' @description Extract F0 and glottal closure instants using the Robust Epoch
##'   And Pitch EstimatoR (REAPER) from SPTK. This C++ implementation is 2-3x
##'   faster than the Python version.
##'
##'   REAPER provides both F0 estimates and epoch marks (glottal closure instants),
##'   making it useful for voice source analysis.
##'
##' @inheritParams rapt
##' @param voicing_threshold Voicing threshold (default: 0.9)
##'
##' @return If toFile=TRUE, returns the number of successfully processed files.
##'   If toFile=FALSE, returns AsspDataObj or list of AsspDataObj objects with
##'   f0 track and optionally epochs.
##'
##' @export
##' @examples
##' \dontrun{
##' # Extract F0 and epochs
##' reaper("recording.wav")
##'
##' # Get epochs for voice source analysis
##' result <- reaper("speech.wav", toFile = FALSE)
##' epochs <- attr(result, "epochs")
##' }
reaper <- function(listOfFiles = NULL,
                   beginTime = 0.0,
                   endTime = 0.0,
                   windowShift = 10.0,
                   minF = 60.0,
                   maxF = 400.0,
                   voicing_threshold = 0.9,
                   toFile = TRUE,
                   explicitExt = "f0",
                   outputDirectory = NULL,
                   verbose = TRUE) {

  if (is.null(listOfFiles) || length(listOfFiles) == 0) {
    cli::cli_abort("No input files specified in {.arg listOfFiles}")
  }

  listOfFiles <- fast_strip_file_protocol(listOfFiles)
  listOfFiles <- normalizePath(path.expand(listOfFiles), mustWork = FALSE)

  files_exist <- file.exists(listOfFiles)
  if (!all(files_exist)) {
    missing_files <- listOfFiles[!files_exist]
    cli::cli_abort(c(
      "!" = "Some files do not exist:",
      "x" = "{.file {fast_basename(missing_files)}}"
    ))
  }

  n_files <- length(listOfFiles)

  beginTime <- if (is.null(beginTime)) 0.0 else beginTime
  endTime <- if (is.null(endTime)) 0.0 else endTime

  if (length(beginTime) == 1) beginTime <- rep(beginTime, n_files)
  if (length(endTime) == 1) endTime <- rep(endTime, n_files)

  makeOutputDirectory(outputDirectory, FALSE, "reaper")

  if (verbose) {
    cli::cli_inform("Applying {.fun reaper} (C++) to {cli::no(n_files)} recording{?s}")
  }

  results <- vector("list", n_files)

  if (verbose && n_files > 1) {
    cli::cli_progress_bar(
      "Processing files",
      total = n_files,
      format = "{cli::pb_spin} {cli::pb_current}/{cli::pb_total} | ETA: {cli::pb_eta}"
    )
  }

  for (i in seq_len(n_files)) {
    file_path <- listOfFiles[i]
    bt <- beginTime[i]
    et <- endTime[i]

    tryCatch({
      audio_obj <- av_to_asspDataObj(
        file_path,
        start_time = bt,
        end_time = if (et == 0.0) NULL else et
      )

      reaper_result <- reaper_cpp(
        audio_obj = audio_obj,
        minF = minF,
        maxF = maxF,
        windowShift = windowShift,
        voicing_threshold = voicing_threshold,
        verbose = FALSE
      )

      out_obj <- create_f0_asspobj(reaper_result, windowShift)

      # Store epochs and polarity as attributes
      attr(out_obj, "epochs") <- reaper_result$epochs
      attr(out_obj, "n_epochs") <- reaper_result$n_epochs
      attr(out_obj, "polarity") <- reaper_result$polarity

      if (toFile) {
        out_file <- generate_output_path(file_path, explicitExt, outputDirectory)
        write.AsspDataObj(out_obj, out_file)
        results[[i]] <- TRUE
      } else {
        results[[i]] <- out_obj
      }

    }, error = function(e) {
      cli::cli_warn("Error processing {.file {basename(file_path)}}: {conditionMessage(e)}")
      results[[i]] <- if (toFile) FALSE else NULL
    })

    if (verbose && n_files > 1) {
      cli::cli_progress_update()
    }
  }

  if (verbose && n_files > 1) {
    cli::cli_progress_done()
  }

  if (toFile) {
    n_success <- sum(unlist(results), na.rm = TRUE)
    if (verbose) {
      cli::cli_inform("Successfully processed {n_success} of {n_files} file{?s}")
    }
    return(invisible(n_success))
  } else {
    results <- results[!sapply(results, is.null)]
    if (length(results) == 1) {
      return(results[[1]])
    } else {
      return(results)
    }
  }
}

attr(reaper, "ext") <- "f0"
attr(reaper, "tracks") <- c("f0")
attr(reaper, "outputType") <- "SSFF"
attr(reaper, "nativeFiletypes") <- c("wav", "flac", "mp3", "mp4", "mkv", "avi")
attr(reaper, "suggestCaching") <- FALSE


##' DIO Pitch Tracking (C++ implementation)
##'
##' @description Extract F0 using the DIO algorithm from the WORLD vocoder (via SPTK).
##'   DIO is designed for high-quality pitch extraction for speech synthesis applications.
##'
##' @inheritParams rapt
##' @param voicing_threshold Voicing threshold (default: 0.85)
##'
##' @return If toFile=TRUE, returns the number of successfully processed files.
##'   If toFile=FALSE, returns AsspDataObj or list of AsspDataObj objects.
##'
##' @export
##' @examples
##' \dontrun{
##' # Extract F0 using DIO
##' dio("recording.wav")
##'
##' # Process with custom F0 range
##' dio("speech.wav", minF = 80, maxF = 350)
##' }
dio <- function(listOfFiles = NULL,
                beginTime = 0.0,
                endTime = 0.0,
                windowShift = 10.0,
                minF = 60.0,
                maxF = 400.0,
                voicing_threshold = 0.85,
                toFile = TRUE,
                explicitExt = "f0",
                outputDirectory = NULL,
                verbose = TRUE) {

  if (is.null(listOfFiles) || length(listOfFiles) == 0) {
    cli::cli_abort("No input files specified in {.arg listOfFiles}")
  }

  listOfFiles <- fast_strip_file_protocol(listOfFiles)
  listOfFiles <- normalizePath(path.expand(listOfFiles), mustWork = FALSE)

  files_exist <- file.exists(listOfFiles)
  if (!all(files_exist)) {
    missing_files <- listOfFiles[!files_exist]
    cli::cli_abort(c(
      "!" = "Some files do not exist:",
      "x" = "{.file {fast_basename(missing_files)}}"
    ))
  }

  n_files <- length(listOfFiles)

  beginTime <- if (is.null(beginTime)) 0.0 else beginTime
  endTime <- if (is.null(endTime)) 0.0 else endTime

  if (length(beginTime) == 1) beginTime <- rep(beginTime, n_files)
  if (length(endTime) == 1) endTime <- rep(endTime, n_files)

  makeOutputDirectory(outputDirectory, FALSE, "dio")

  if (verbose) {
    cli::cli_inform("Applying {.fun dio} (C++) to {cli::no(n_files)} recording{?s}")
  }

  results <- vector("list", n_files)

  if (verbose && n_files > 1) {
    cli::cli_progress_bar(
      "Processing files",
      total = n_files,
      format = "{cli::pb_spin} {cli::pb_current}/{cli::pb_total} | ETA: {cli::pb_eta}"
    )
  }

  for (i in seq_len(n_files)) {
    file_path <- listOfFiles[i]
    bt <- beginTime[i]
    et <- endTime[i]

    tryCatch({
      audio_obj <- av_to_asspDataObj(
        file_path,
        start_time = bt,
        end_time = if (et == 0.0) NULL else et
      )

      dio_result <- dio_cpp(
        audio_obj = audio_obj,
        minF = minF,
        maxF = maxF,
        windowShift = windowShift,
        voicing_threshold = voicing_threshold,
        verbose = FALSE
      )

      out_obj <- create_f0_asspobj(dio_result, windowShift)

      if (toFile) {
        out_file <- generate_output_path(file_path, explicitExt, outputDirectory)
        write.AsspDataObj(out_obj, out_file)
        results[[i]] <- TRUE
      } else {
        results[[i]] <- out_obj
      }

    }, error = function(e) {
      cli::cli_warn("Error processing {.file {basename(file_path)}}: {conditionMessage(e)}")
      results[[i]] <- if (toFile) FALSE else NULL
    })

    if (verbose && n_files > 1) {
      cli::cli_progress_update()
    }
  }

  if (verbose && n_files > 1) {
    cli::cli_progress_done()
  }

  if (toFile) {
    n_success <- sum(unlist(results), na.rm = TRUE)
    if (verbose) {
      cli::cli_inform("Successfully processed {n_success} of {n_files} file{?s}")
    }
    return(invisible(n_success))
  } else {
    results <- results[!sapply(results, is.null)]
    if (length(results) == 1) {
      return(results[[1]])
    } else {
      return(results)
    }
  }
}

attr(dio, "ext") <- "f0"
attr(dio, "tracks") <- c("f0")
attr(dio, "outputType") <- "SSFF"
attr(dio, "nativeFiletypes") <- c("wav", "flac", "mp3", "mp4", "mkv", "avi")
attr(dio, "suggestCaching") <- FALSE


##' Convert SPTK C++ pitch result to AsspDataObj
##'
##' @param pitch_result List returned from rapt_cpp, swipe_cpp, reaper_cpp, or dio_cpp
##' @param windowShift Frame shift in milliseconds
##' @return AsspDataObj with f0 track
##' @keywords internal
create_f0_asspobj <- function(pitch_result, windowShift) {
  # Extract F0 matrix and metadata
  f0_matrix <- pitch_result$f0
  sample_rate <- pitch_result$sample_rate
  n_frames <- pitch_result$n_frames

  # Create AsspDataObj
  out_obj <- list(
    f0 = f0_matrix
  )

  # Set attributes
  attr(out_obj, "sampleRate") <- sample_rate
  attr(out_obj, "startTime") <- 0.0
  attr(out_obj, "startRecord") <- 1L
  attr(out_obj, "endRecord") <- n_frames
  attr(out_obj, "windowShift") <- windowShift / 1000.0  # Convert ms to seconds

  class(out_obj) <- c("AsspDataObj", "list")

  return(out_obj)
}


##' Generate output file path for SSFF files
##'
##' @param input_file Input file path
##' @param ext Output extension
##' @param output_dir Output directory (optional)
##' @return Output file path
##' @keywords internal
generate_output_path <- function(input_file, ext, output_dir = NULL) {
  if (is.null(output_dir)) {
    out_dir <- dirname(input_file)
  } else {
    out_dir <- output_dir
  }

  base_name <- fast_file_path_sans_ext(c(basename(input_file)))[1]
  output_file <- file.path(out_dir, paste0(base_name, ".", ext))

  return(output_file)
}
