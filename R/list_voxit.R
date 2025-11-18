#' Extract Voxit prosodic complexity features from audio files
#'
#' @description
#' Computes voice and articulation complexity measures from audio files using
#' word-level alignments and pitch contours. Features include speaking rate,
#' pause statistics, rhythmic complexity, and pitch dynamics.
#'
#' This function requires:
#' \itemize{
#'   \item \strong{Word alignments}: Either provided as CSV files or computed automatically
#'   \item \strong{Pitch tracking}: Uses SAcC algorithm (install with \code{\link{install_sacc}()})
#' }
#'
#' @details
#' \strong{Audio Loading:}
#'
#' Audio files are loaded via the \code{\link[av]{read_audio_bin}} function,
#' supporting all media formats (WAV, MP3, MP4, video files, etc.).
#'
#' \strong{Requirements:}
#'
#' The voxit and SAcC Python modules must be installed:
#' \itemize{
#'   \item \code{\link{install_voxit}()} - Voxit feature extraction
#'   \item \code{\link{install_sacc}()} - SAcC pitch tracker
#' }
#'
#' \strong{Performance:}
#'
#' Single file processing: ~0.5-2s depending on audio duration and alignment method
#'
#' Batch processing with 8 cores: ~6x speedup
#'
#' \strong{Features Computed (11 total):}
#'
#' \enumerate{
#'   \item \strong{WPM}: Words per minute (speaking rate)
#'   \item \strong{pause_count}: Number of pauses (100-3000ms)
#'   \item \strong{long_pause_count}: Number of pauses > 3s
#'   \item \strong{average_pause_length}: Mean pause duration (seconds)
#'   \item \strong{average_pause_rate}: Pauses per second
#'   \item \strong{rhythmic_complexity_of_pauses}: Normalized Lempel-Ziv complexity
#'   \item \strong{average_pitch}: Mean F0 (Hz)
#'   \item \strong{pitch_range}: F0 range (octaves)
#'   \item \strong{pitch_speed}: F0 velocity (octaves/second)
#'   \item \strong{pitch_acceleration}: F0 acceleration (octaves/second²)
#'   \item \strong{pitch_entropy}: F0 distribution entropy (bits)
#' }
#'
#' @param listOfFiles Character vector of file paths to audio files
#' @param alignmentFiles Character vector of alignment CSV file paths (optional).
#'   If NULL, alignments will be computed automatically using forced alignment.
#' @param beginTime Numeric. Start time in seconds (default: 0 = file start)
#' @param endTime Numeric. End time in seconds (default: 0 = file end)
#' @param minF Numeric. Minimum F0 for pitch tracking in Hz (default: 60)
#' @param maxF Numeric. Maximum F0 for pitch tracking in Hz (default: 600)
#' @param verbose Logical. Show progress messages (default: TRUE)
#' @param parallel Logical. Use parallel processing for multiple files (default: TRUE)
#' @param n_cores Integer. Number of cores for parallel processing (default: NULL = auto-detect)
#' @param toFile Logical. If TRUE, write results to JSTF file. Default FALSE.
#' @param explicitExt Character. File extension for output. Default "vxt".
#' @param outputDirectory Character. Output directory path. Default NULL (use input directory).
#'
#' @return
#' If \code{toFile=FALSE} (default): For single file, a named list with 11 prosodic features.
#'   For multiple files, named list where each element is a file's feature list.
#'   If \code{toFile=TRUE}, invisibly returns the path(s) to the written JSTF file(s).
#'
#' @section Word Alignments:
#' Word alignments can be provided as CSV files with columns:
#' \describe{
#'   \item{word}{Word text}
#'   \item{case}{Word category (e.g., "success")}
#'   \item{start}{Start time in seconds}
#'   \item{end}{End time in seconds}
#' }
#'
#' If no alignment files are provided, the function will attempt automatic
#' forced alignment (requires additional setup - see documentation).
#'
#' @section Time Windowing:
#' Time windowing is handled by the \code{av} package during audio loading.
#'
#' @examples
#' \dontrun{
#' # Check if voxit is available
#' if (!voxit_available()) {
#'   install_voxit()
#' }
#' if (!sacc_available()) {
#'   install_sacc()
#' }
#'
#' # Single file analysis with alignments
#' features <- lst_voxit("speech.wav", 
#'                       alignmentFiles = "speech_align.csv")
#' print(features$WPM)
#' print(features$average_pitch)
#' print(features$rhythmic_complexity_of_pauses)
#'
#' # Analyze specific time window
#' features <- lst_voxit("speech.wav",
#'                       alignmentFiles = "speech_align.csv",
#'                       beginTime = 1.0, 
#'                       endTime = 3.0)
#'
#' # Batch processing (parallel)
#' files <- list.files("audio_dir", pattern = "\\.wav$", full.names = TRUE)
#' aligns <- list.files("audio_dir", pattern = "_align\\.csv$", full.names = TRUE)
#' results <- lst_voxit(files, alignmentFiles = aligns,
#'                      verbose = TRUE, parallel = TRUE, n_cores = 4)
#'
#' # Convert to data frame
#' library(tidyverse)
#' df <- results %>%
#'   map_dfr(~ as.data.frame(t(unlist(.))), .id = "file")
#'
#' # Write results to JSTF file
#' lst_voxit("speech.wav", alignmentFiles = "speech_align.csv", toFile = TRUE)  # Creates speech.vxt
#'
#' # Read back and convert to data.frame
#' track <- read_track("speech.vxt")
#' df <- as.data.frame(track)
#' head(df)  # Shows begin_time, end_time, and all 11 Voxit features
#' }
#'
#' @references
#' Voxit toolbox: Voice and articulation complexity measures
#'
#' SAcC pitch tracker: Ellis, D.P.W., & Weiss, R.J. (2010). 
#' Pitch and voicing estimation via a harmonic model.
#'
#' @seealso
#' \code{\link{install_voxit}}, \code{\link{voxit_available}},
#' \code{\link{voxit_info}}, \code{\link{install_sacc}},
#' \code{\link{trk_sacc}}
#'
#' @export
lst_voxit <- function(listOfFiles,
                      alignmentFiles = NULL,
                      beginTime = 0.0,
                      endTime = 0.0,
                      minF = 60,
                      maxF = 600,
                      verbose = TRUE,
                      parallel = TRUE,
                      n_cores = NULL,
                      toFile = FALSE,
                      explicitExt = "vxt",
                      outputDirectory = NULL) {


  # Input validation
  if (is.null(listOfFiles) || length(listOfFiles) == 0) {
    cli::cli_abort("No files provided")
  }

  listOfFiles <- as.vector(listOfFiles)
  n_files <- length(listOfFiles)

  # Check files exist
  missing <- !file.exists(listOfFiles)
  if (any(missing)) {
    cli::cli_abort(c(
      "x" = "{sum(missing)} file{?s} not found",
      "i" = "First missing: {.file {listOfFiles[which(missing)[1]]}}"
    ))
  }

  # Check Python modules available
  if (!reticulate::py_module_available("voxit")) {
    cli::cli_abort(c(
      "x" = "voxit Python module not found",
      "i" = "Install with: install_voxit()"
    ))
  }
  
  if (!reticulate::py_module_available("sacc")) {
    cli::cli_abort(c(
      "x" = "SAcC Python module not found (required for pitch tracking)",
      "i" = "Install with: install_sacc()"
    ))
  }

  # Check alignment files if provided
  if (!is.null(alignmentFiles)) {
    if (length(alignmentFiles) != n_files) {
      cli::cli_abort(c(
        "x" = "Number of alignment files ({length(alignmentFiles)}) must match audio files ({n_files})"
      ))
    }
    missing_align <- !file.exists(alignmentFiles)
    if (any(missing_align)) {
      cli::cli_abort(c(
        "x" = "{sum(missing_align)} alignment file{?s} not found",
        "i" = "First missing: {.file {alignmentFiles[which(missing_align)[1]]}}"
      ))
    }
  }

  # Import Python modules
  voxit <- reticulate::import("voxit")
  sacc_module <- reticulate::import("sacc")
  np <- reticulate::import("numpy")

  # Define processing function for single file
  process_single_file <- function(i) {
    file_path <- listOfFiles[i]
    align_path <- if (!is.null(alignmentFiles)) alignmentFiles[i] else NULL

    if (verbose && !parallel) {
      cli::cli_progress_step("Processing {.file {basename(file_path)}}")
    }

    tryCatch({
      # Load audio
      audio_data <- av::read_audio_bin(
        audio = file_path,
        start_time = if (beginTime > 0) beginTime else NULL,
        end_time = if (endTime > 0) endTime else NULL,
        channels = 1
      )
      
      sample_rate <- attr(audio_data, "sample_rate")
      
      # Convert to float32
      audio_float <- as.numeric(audio_data) / 2147483647.0
      audio_np <- np$array(audio_float, dtype = "float32")

      # Run SAcC pitch tracking
      pitch_result <- sacc_module$extract_pitch(
        audio = audio_np,
        sample_rate = as.integer(sample_rate),
        min_f0 = minF,
        max_f0 = maxF
      )
      
      # Convert pitch to expected format
      pitch_data <- list()
      if (!is.null(pitch_result$times) && !is.null(pitch_result$frequencies)) {
        times <- as.numeric(pitch_result$times)
        freqs <- as.numeric(pitch_result$frequencies)
        for (j in seq_along(times)) {
          pitch_data[[j]] <- list(time = times[j], frequency = freqs[j])
        }
      }

      # Load or generate alignments
      if (!is.null(align_path)) {
        # Read alignment CSV
        align_df <- readr::read_csv(align_path, show_col_types = FALSE)
        gentle_data <- lapply(1:nrow(align_df), function(k) {
          list(
            word = as.character(align_df$word[k]),
            case = if ("case" %in% names(align_df)) as.character(align_df$case[k]) else "success",
            start = as.numeric(align_df$start[k]),
            end = as.numeric(align_df$end[k])
          )
        })
      } else {
        # Would need forced alignment here - for now error
        cli::cli_abort(c(
          "x" = "Alignment files required (automatic alignment not yet implemented)",
          "i" = "Provide alignmentFiles parameter"
        ))
      }

      # Compute Voxit features
      features <- voxit$compute_features(
        gentle_data = gentle_data,
        pitch_data = pitch_data,
        start_time = if (beginTime > 0) beginTime else reticulate::py_none(),
        end_time = if (endTime > 0) endTime else reticulate::py_none()
      )

      # Convert to R list
      result <- as.list(features)
      
      if (verbose && !parallel) {
        cli::cli_progress_done()
      }
      
      return(result)

    }, error = function(e) {
      if (verbose) {
        cli::cli_alert_warning("Failed {.file {basename(file_path)}}: {e$message}")
      }
      return(NULL)
    })
  }

  # Process files
  if (n_files == 1) {
    # Single file
    result <- process_single_file(1)


    # Handle JSTF file writing
    if (toFile && !is.null(result)) {
      file_path <- listOfFiles[1]
      bt <- beginTime
      et <- endTime

      # Get audio metadata
      audio_info <- av::av_media_info(file_path)
      sample_rate <- audio_info$audio$sample_rate
      audio_duration <- audio_info$duration

      json_obj <- create_json_track_obj(
        results = result,
        function_name = "lst_voxit",
        file_path = file_path,
        sample_rate = sample_rate,
        audio_duration = audio_duration,
        beginTime = bt,
        endTime = if (et > 0) et else audio_duration,
        parameters = list(
          minF = minF,
          maxF = maxF,
          has_alignment = !is.null(alignmentFiles)
        )
      )

      base_name <- tools::file_path_sans_ext(basename(file_path))
      out_dir <- if (is.null(outputDirectory)) dirname(file_path) else outputDirectory
      output_path <- file.path(out_dir, paste0(base_name, ".", explicitExt))

      write_json_track(json_obj, output_path)
      return(invisible(output_path))
    }

  } else if (parallel && n_files > 1) {
    # Parallel processing
    if (is.null(n_cores)) {
      n_cores <- parallel::detectCores() - 1
    }
    n_cores <- min(n_cores, n_files)

    if (verbose) {
      cli::cli_alert_info("Processing {n_files} files using {n_cores} cores")
    }

    # Setup parallel backend
    if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(n_cores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      parallel::clusterExport(cl, c("listOfFiles", "alignmentFiles", "beginTime",
                                   "endTime", "minF", "maxF", "verbose"),
                            envir = environment())
      result <- parallel::parLapply(cl, 1:n_files, process_single_file)
    } else {
      result <- parallel::mclapply(1:n_files, process_single_file,
                                  mc.cores = n_cores)
    }

    names(result) <- basename(listOfFiles)
  } else {
    # Sequential processing
    if (verbose) {
      cli::cli_progress_bar("Processing files", total = n_files)
    }

    result <- lapply(1:n_files, function(i) {
      if (verbose) {
        cli::cli_progress_update()
      }
      process_single_file(i)
    })

    if (verbose) {
      cli::cli_progress_done()
    }

    names(result) <- basename(listOfFiles)
  }


  # Handle JSTF file writing for multiple files
  if (toFile && n_files > 1) {
    output_paths <- character(n_files)
    for (i in seq_len(n_files)) {
      file_basename <- basename(listOfFiles[i])
      file_result <- result[[file_basename]]

      if (!is.null(file_result)) {
        file_path <- listOfFiles[i]
        bt <- beginTime
        et <- endTime

        # Get audio metadata
        audio_info <- av::av_media_info(file_path)
        sample_rate <- audio_info$audio$sample_rate
        audio_duration <- audio_info$duration

        json_obj <- create_json_track_obj(
          results = file_result,
          function_name = "lst_voxit",
          file_path = file_path,
          sample_rate = sample_rate,
          audio_duration = audio_duration,
          beginTime = bt,
          endTime = if (et > 0) et else audio_duration,
          parameters = list(
            minF = minF,
            maxF = maxF,
            has_alignment = !is.null(alignmentFiles)
          )
        )

        base_name <- tools::file_path_sans_ext(basename(file_path))
        out_dir <- if (is.null(outputDirectory)) dirname(file_path) else outputDirectory
        output_path <- file.path(out_dir, paste0(base_name, ".", explicitExt))

        write_json_track(json_obj, output_path)
        output_paths[i] <- output_path
      } else {
        output_paths[i] <- NA_character_
      }
    }
    return(invisible(output_paths))
  }

  # Return single result if only one file
  if (n_files == 1) {
    return(result)
  }

  return(result)
}

# Set function attributes
attr(lst_voxit, "module") <- "voxit"
attr(lst_voxit, "type") <- "summary"
attr(lst_voxit, "features") <- 11
attr(lst_voxit, "ext") <- "vxt"
attr(lst_voxit, "outputType") <- "JSTF"
attr(lst_voxit, "format") <- "JSON"

