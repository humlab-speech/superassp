#' Convert audio data to AsspDataObj with automatic format fallback
#'
#' This function reads audio data using the av package and converts it to
#' an AsspDataObj that can be processed by superassp analysis functions without
#' writing intermediate files to disk.
#'
#' For niche audio formats not supported by av (au, kay, nist, nsp, csre, ssff),
#' the function automatically falls back to wrassp's read.AsspDataObj().
#'
#' @param file_path Path to the audio/video file
#' @param start_time Start time in seconds (default 0)
#' @param end_time End time in seconds (default NULL for end of file)
#' @param target_sample_rate Target sample rate (default NULL to keep original).
#'   Note: Resampling is not supported for wrassp fallback formats.
#'
#' @return An AsspDataObj containing the audio data
#'
#' @details
#' \strong{Format Support:}
#' \itemize{
#'   \item \strong{via av}: wav, mp3, mp4, m4a, flac, ogg, aac, opus, wma, and most video formats
#'   \item \strong{via wrassp fallback}: au, kay, nist, nsp, csre, ssff
#' }
#'
#' The function first attempts to read using av for maximum format support and features.
#' If av fails (typically for niche formats), it automatically falls back to wrassp.
#'
#' \strong{Feature Support by Method:}
#' \itemize{
#'   \item \strong{av}: Full support for time windowing and resampling
#'   \item \strong{wrassp}: Time windowing supported, resampling not available
#' }
#'
#' @examples
#' \dontrun{
#' # Read common format using av
#' audio_obj <- av_to_asspDataObj("myfile.mp3")
#'
#' # Read niche format (automatic wrassp fallback)
#' audio_obj <- av_to_asspDataObj("myfile.nist")
#'
#' # With time windowing
#' audio_obj <- av_to_asspDataObj("myfile.wav", start_time = 1.0, end_time = 3.0)
#'
#' # Perform RMS analysis without intermediate files
#' rms_result <- rmsana_memory(audio_obj)
#' }
av_to_asspDataObj <- function(file_path, start_time = 0, end_time = NULL,
                               target_sample_rate = NULL) {

  if (!requireNamespace("av", quietly = TRUE)) {
    stop("Package 'av' is required but not installed. Please install it with: install.packages('av')")
  }

  # Try av first, fall back to wrassp for niche formats
  av_success <- FALSE

  tryCatch({
    # Read audio info first
    info <- av::av_media_info(file_path)

    if (length(info$audio) == 0) {
      stop("No audio stream found in file: ", file_path)
    }

    audio_info <- info$audio
    original_sample_rate <- audio_info$sample_rate
    channels <- audio_info$channels

    # Use original sample rate if not specified
    if (is.null(target_sample_rate)) {
      target_sample_rate <- original_sample_rate
    }

    # Calculate time window
    duration <- info$duration
    if (is.null(end_time)) {
      end_time <- duration
    }

    # Validate time window
    if (start_time < 0) start_time <- 0
    if (end_time > duration) end_time <- duration
    if (start_time >= end_time) {
      stop("Invalid time window: start_time (", start_time,
           ") >= end_time (", end_time, ")")
    }

    av_success <- TRUE

  }, error = function(e) {
    # av failed - will try wrassp fallback below
    NULL
  })

  # If av succeeded, read with av
  if (av_success) {
    # Read audio data using av
    # av::read_audio_bin returns 32-bit signed integers (s32le format)
    audio_data <- av::read_audio_bin(file_path,
                                      channels = channels,
                                      start_time = start_time,
                                      end_time = end_time,
                                      sample_rate = target_sample_rate)

    # audio_data is an integer vector with interleaved samples
    # av returns s32le (32-bit signed integers)
    # We need to convert to 16-bit for AsspDataObj
    # Scale from 32-bit range to 16-bit range
    samples_int16 <- as.integer(audio_data / 65536)

    # De-interleave channels if multi-channel
    if (channels > 1) {
      n_frames <- length(samples_int16) / channels
      sample_matrix <- matrix(samples_int16, nrow = n_frames, ncol = channels, byrow = TRUE)
    } else {
      sample_matrix <- matrix(samples_int16, ncol = 1)
    }

    # Create AsspDataObj matching the structure from read.AsspDataObj
    n_frames <- nrow(sample_matrix)

    result <- list()
    result$audio <- sample_matrix

    # Set attributes with correct types (numeric for rates/times, integer for records)
    attr(result, "sampleRate") <- as.numeric(target_sample_rate)
    attr(result, "trackFormats") <- c("INT16")  # Must be a vector, not a single string
    attr(result, "startTime") <- as.numeric(start_time)
    attr(result, "startRecord") <- 1L
    attr(result, "endRecord") <- as.integer(n_frames)
    attr(result, "origFreq") <- 0.0  # For WAVE files, origFreq should be 0
    attr(result, "filePath") <- file_path
    attr(result, "fileInfo") <- c(21L, 2L)  # FF_WAVE=21, FDF_BIN=2

    class(result) <- "AsspDataObj"

    return(result)

  } else {
    # av failed - try wrassp fallback for niche formats (au, kay, nist, nsp)
    # Check if wrassp can handle this format
    file_ext <- tolower(tools::file_ext(file_path))

    # Formats supported by wrassp but not by av
    wrassp_formats <- c("au", "kay", "nist", "nsp", "csre", "ssff")

    if (!(file_ext %in% wrassp_formats)) {
      stop("Audio format '", file_ext, "' is not supported by av package. ",
           "Supported niche formats via wrassp: ", paste(wrassp_formats, collapse = ", "),
           call. = FALSE)
    }

    # Use wrassp to read the file
    tryCatch({
      result <- read.AsspDataObj(file_path)

      # Handle time windowing if requested
      # wrassp doesn't support time windowing directly, so we need to extract frames
      if (start_time > 0 || !is.null(end_time)) {
        sample_rate <- attr(result, "sampleRate")

        if (is.null(end_time)) {
          # Get file duration from result
          n_samples <- nrow(result$audio)
          end_time <- n_samples / sample_rate
        }

        # Convert times to sample indices
        start_sample <- max(1, floor(start_time * sample_rate) + 1)
        end_sample <- min(nrow(result$audio), floor(end_time * sample_rate))

        # Extract subset
        result$audio <- result$audio[start_sample:end_sample, , drop = FALSE]

        # Update attributes
        attr(result, "startTime") <- as.numeric(start_time)
        attr(result, "startRecord") <- 1L
        attr(result, "endRecord") <- as.integer(nrow(result$audio))
      }

      # Handle resampling if requested
      if (!is.null(target_sample_rate)) {
        current_rate <- attr(result, "sampleRate")

        if (current_rate != target_sample_rate) {
          warning("Resampling from ", current_rate, " Hz to ", target_sample_rate,
                  " Hz is not supported for wrassp formats. ",
                  "Using original sample rate.", call. = FALSE)
        }
      }

      return(result)

    }, error = function(e) {
      stop("Failed to read audio file with both av and wrassp: ", e$message, call. = FALSE)
    })
  }
}


#' Perform RMS analysis on AsspDataObj in memory
#'
#' This function performs RMS analysis on an AsspDataObj without writing
#' intermediate files. It's useful when processing audio data read from
#' video files or other non-WAV formats using the av package.
#'
#' @param audio_obj AsspDataObj containing audio data
#' @param ... Additional parameters passed to rmsana
#'
#' @return AsspDataObj with RMS analysis results
#' @examples
#' \dontrun{
#' audio_obj <- av_to_asspDataObj("myfile.mp4")
#' rms_result <- rmsana_memory(audio_obj, windowShift = 5)
#' }
rmsana_memory <- function(audio_obj, ...) {

  if (!inherits(audio_obj, "AsspDataObj")) {
    stop("audio_obj must be an AsspDataObj. Use av_to_asspDataObj() to convert.")
  }

  # Call performAsspMemory for true in-memory processing (no temp files!)
  result <- .External(
    "performAsspMemory", audio_obj,
    fname = "trk_rmsana",
    toFile = FALSE,
    ...,
    PACKAGE = "superassp"
  )

  return(result)
}


#' Process audio from any media file format
#'
#' This function provides a convenient interface to process audio from any
#' media file format supported by av (including video files), perform
#' superassp analysis, and return results without creating intermediate files.
#'
#' @param file_path Path to the media file
#' @param analysis_function Name of analysis function (e.g., "trk_rmsana", "trk_forest")
#' @param start_time Start time in seconds (default 0)
#' @param end_time End time in seconds (default NULL for end)
#' @param target_sample_rate Target sample rate (default NULL for original)
#' @param ... Additional parameters for the analysis function
#'
#' @return Result from the analysis function
#' @examples
#' \dontrun{
#' # Analyze RMS from a video file
#' rms <- process_media_file("video.mp4", "trk_rmsana", windowShift = 5)
#'
#' # Analyze formants from audio segment
#' formants <- process_media_file("audio.m4a", "trk_forest",
#'                                 start_time = 10, end_time = 20)
#' }
process_media_file <- function(file_path, analysis_function = "trk_rmsana",
                                start_time = 0, end_time = NULL,
                                target_sample_rate = NULL, ...) {

  # Convert media to AsspDataObj
  audio_obj <- av_to_asspDataObj(file_path, start_time, end_time, target_sample_rate)

  # Call performAsspMemory for true in-memory processing (no temp files!)
  result <- .External(
    "performAsspMemory", audio_obj,
    fname = analysis_function,
    toFile = FALSE,
    ...,
    PACKAGE = "superassp"
  )

  return(result)
}


#' Process media files with performAssp using av package (load-and-process pattern)
#'
#' This internal function replaces the convert-then-process pattern with a
#' load-and-process pattern. Instead of converting media files to WAV on disk,
#' it loads them directly into memory using av, then processes them with performAssp.
#'
#' Supports parallel processing for batch operations on multi-core systems.
#'
#' @param listOfFiles Character vector of input file paths
#' @param beginTime Numeric vector of begin times (seconds)
#' @param endTime Numeric vector of end times (seconds, 0 = end of file)
#' @param nativeFiletypes Character vector of natively supported formats
#' @param fname Character name of performAssp function to call
#' @param toFile Logical whether to write output files
#' @param verbose Logical whether to show progress messages
#' @param parallel Logical whether to use parallel processing (default TRUE for >1 files)
#' @param n_cores Integer number of cores to use (default: detectCores() - 1)
#' @param ... Additional parameters to pass to performAssp
#'
#' @return List with:
#'   - externalRes: Results from performAssp
#'   - listOfFilesDF: Data frame with file processing information
#'   - processed_native: Logical vector indicating which files were native
#'
#' @keywords internal
processMediaFiles_LoadAndProcess <- function(listOfFiles, beginTime, endTime,
                                             nativeFiletypes, fname,
                                             toFile = TRUE, verbose = TRUE, ...) {

  # Extract parallel processing parameters from ... and any other DSP parameters
  dots <- list(...)
  parallel <- dots$parallel
  n_cores <- dots$n_cores

  # Remove parallel-specific params so they don't get passed to .External
  dsp_params <- dots
  dsp_params$parallel <- NULL
  dsp_params$n_cores <- NULL

  n_files <- length(listOfFiles)

  # Normalize paths
  listOfFiles <- fast_strip_file_protocol(listOfFiles)
  listOfFiles <- normalizePath(path.expand(listOfFiles))

  # Ensure time vectors match file count
  if(length(beginTime) == 1) beginTime <- rep(beginTime, n_files)
  if(length(endTime) == 1) endTime <- rep(endTime, n_files)

  # Determine which files are native format
  file_exts <- fast_file_ext(listOfFiles)
  is_native <- fast_is_native(file_exts, nativeFiletypes)

  # Check for lossless formats - ESSENTIAL for accurate DSP
  known_lossless <- knownLossless()
  is_lossless <- tolower(file_exts) %in% tolower(known_lossless)
  not_lossless <- listOfFiles[!is_lossless]

  if(length(not_lossless) > 0 && verbose) {
    lossless_msg <- if(all(beginTime == 0) && all(endTime == 0)) {
      c(
        "!" = "Found {length(not_lossless)} recording{?s} in lossy format{?s}",
        "i" = "Lossy compression may affect {.fun {fname}} accuracy",
        "x" = "For accurate DSP, use lossless formats: {.val {intersect(nativeFiletypes, known_lossless)}}"
      )
    } else {
      c(
        "!" = "Found {length(not_lossless)} recording{?s} with lossy compression",
        "i" = "May affect time window extraction and {.fun {fname}} accuracy",
        "x" = "For accurate DSP, use lossless formats: {.val {intersect(nativeFiletypes, known_lossless)}}"
      )
    }

    cli::cli_warn(lossless_msg)
  }

  # Auto-enable parallel processing for batches (unless explicitly disabled)
  if(is.null(parallel)) {
    parallel <- n_files > 1
  }

  # Determine number of cores
  if(is.null(n_cores)) {
    n_cores <- parallel::detectCores() - 1
    if(is.na(n_cores) || n_cores < 1) n_cores <- 1
  }

  # Disable parallel for single file or when explicitly set to FALSE
  use_parallel <- parallel && n_files > 1 && n_cores > 1

  # Define the processing function for a single file
  # This function is thread-safe as it creates independent objects
  process_single_file <- function(i) {
    file_path <- listOfFiles[i]
    bt <- beginTime[i]
    et <- endTime[i]

    # Load audio with av (handles all formats and time windowing)
    # av::read_audio_bin is thread-safe - each call reads independently
    audio_obj <- av_to_asspDataObj(
      file_path,
      start_time = bt,
      end_time = if(et == 0.0) NULL else et,
      target_sample_rate = NULL
    )

    # Call performAsspMemory for true in-memory processing
    # performAsspMemory is thread-safe - uses stack-allocated AOPTS
    # and creates independent DOBJ structures per call
    result <- do.call(
      .External,
      c(list("performAsspMemory", audio_obj,
             fname = fname,
             toFile = toFile,
             progressBar = NULL,
             PACKAGE = "superassp"),
        dsp_params)
    )

    return(result)
  }

  # Process files (parallel or sequential)
  if(use_parallel) {
    if(verbose) {
      cli::cli_inform(c(
        "i" = "Processing {n_files} file{?s} in parallel on {n_cores} core{?s}"
      ))
    }

    # Use platform-appropriate parallel backend
    if(.Platform$OS.type == "windows") {
      # Windows: use parallel socket cluster (parLapply)
      cl <- parallel::makeCluster(n_cores)
      on.exit(parallel::stopCluster(cl), add = TRUE)

      # Export necessary functions and variables to cluster
      parallel::clusterExport(cl, c(
        "listOfFiles", "beginTime", "endTime", "fname", "toFile",
        "av_to_asspDataObj", "process_single_file", "dsp_params"
      ), envir = environment())

      # Load required packages on each worker
      parallel::clusterEvalQ(cl, {
        library(superassp)
      })

      # Run in parallel with progress bar if verbose
      if(verbose) {
        externalRes <- pbapply::pblapply(
          seq_along(listOfFiles),
          process_single_file,
          cl = cl
        )
      } else {
        externalRes <- parallel::parLapply(
          cl,
          seq_along(listOfFiles),
          process_single_file
        )
      }
    } else {
      # Unix/Mac: use fork-based parallelism (mclapply)
      # Fork is more efficient as it shares memory (copy-on-write)
      if(verbose) {
        # pbmclapply provides progress bar for mclapply
        if(requireNamespace("pbmcapply", quietly = TRUE)) {
          externalRes <- pbmcapply::pbmclapply(
            seq_along(listOfFiles),
            process_single_file,
            mc.cores = n_cores,
            mc.preschedule = TRUE  # Better load balancing
          )
        } else {
          # Fallback to mclapply without progress bar
          externalRes <- parallel::mclapply(
            seq_along(listOfFiles),
            process_single_file,
            mc.cores = n_cores,
            mc.preschedule = TRUE
          )
        }
      } else {
        externalRes <- parallel::mclapply(
          seq_along(listOfFiles),
          process_single_file,
          mc.cores = n_cores,
          mc.preschedule = TRUE
        )
      }
    }

  } else {
    # Sequential processing with optional progress bar
    if(verbose && n_files > 1) {
      cli::cli_inform("Processing {n_files} file{?s} sequentially")

      # Use cli progress bar for sequential processing
      externalRes <- vector("list", n_files)
      cli::cli_progress_bar(
        "Processing files",
        total = n_files,
        format = "{cli::pb_spin} {cli::pb_current}/{cli::pb_total} | ETA: {cli::pb_eta}"
      )

      for (i in seq_along(listOfFiles)) {
        externalRes[[i]] <- process_single_file(i)
        cli::cli_progress_update()
      }

      cli::cli_progress_done()
    } else {
      # Simple sequential loop (single file or verbose=FALSE)
      externalRes <- vector("list", n_files)
      for (i in seq_along(listOfFiles)) {
        externalRes[[i]] <- process_single_file(i)
      }
    }
  }

  # Build info dataframe
  listOfFilesDF <- data.frame(
    audio = listOfFiles,
    dsp_input = listOfFiles,  # For compatibility
    beginTime = beginTime,
    endTime = endTime,
    lossless = is_lossless,
    stringsAsFactors = FALSE
  )

  return(list(
    externalRes = externalRes,
    listOfFilesDF = listOfFilesDF,
    processed_native = is_native
  ))
}
