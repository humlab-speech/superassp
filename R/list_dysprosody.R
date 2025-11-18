#' Extract dysprosody features from audio files
#'
#' @description
#' Computes 193 prosodic features using the dysprosody model described in
#' Nylén (2025). Features include MOMEL-INTSINT pitch targets, tone labels,
#' spectral tilt measures, and statistical summaries.
#'
#' This function implements the prosodic assessment model described in
#' \insertCite{Nylen.2025.10.3389/fnhum.2025.1566274}{superassp}.
#'
#' The dysprosody model extracts:
#' \itemize{
#'   \item \strong{MOMEL targets} - Pitch targets via quadratic spline modeling
#'   \item \strong{INTSINT coding} - Optimal tone labels (M, T, B, H, L, U, D, S)
#'   \item \strong{Spectral tilt} - Harmonic amplitude measures with formant correction
#'   \item \strong{Statistical summaries} - Mean, std, variation, IQR, min, max for all features
#'   \item \strong{Differential features} - Inter-label changes (suffix: _diff)
#' }
#'
#' @details
#' \strong{Audio Loading:}
#'
#' Audio files are loaded via the \code{\link[av]{read_audio_bin}} function,
#' supporting all media formats (WAV, MP3, MP4, video files, etc.). The audio
#' is converted to a parselmouth Sound object in memory using
#' \code{\link{av_load_for_parselmouth}}, eliminating the need for temporary
#' WAV files. This provides true in-memory processing with improved performance.
#'
#' \strong{Requirements:}
#'
#' The dysprosody Python module must be installed. Install with:
#' \code{\link{install_dysprosody}()}
#'
#' \strong{Performance:}
#'
#' Single file processing: ~0.16-0.44s (14x realtime for 2-6s audio)
#'
#' Batch processing with 8 cores: ~5x speedup
#'
#' Files < 1 second are skipped automatically.
#'
#' \strong{Algorithms:}
#'
#' \enumerate{
#'   \item Automatic F0 range estimation (two-pass method)
#'   \item MOMEL pitch target extraction (quadratic spline modeling)
#'   \item INTSINT tone coding (optimal label assignment)
#'   \item Spectral tilt computation (Iseli-Alwan harmonic correction)
#'   \item Statistical aggregation across INTSINT labels
#' }
#'
#' @param listOfFiles Character vector of file paths to audio files
#' @param beginTime Numeric. Start time in seconds (default: 0 = file start)
#' @param endTime Numeric. End time in seconds (default: 0 = file end)
#' @param minF Numeric. Minimum F0 for pitch extraction in Hz (default: 60)
#' @param maxF Numeric. Maximum F0 for pitch extraction in Hz (default: 750)
#' @param windowShift Numeric. Window shift for intensity in milliseconds (default: 1.0)
#' @param verbose Logical. Show progress messages (default: TRUE)
#' @param parallel Logical. Use parallel processing for multiple files (default: TRUE)
#' @param n_cores Integer. Number of cores for parallel processing (default: NULL = auto-detect)
#' @param toFile Logical. If TRUE, write results to JSTF file. Default FALSE.
#' @param explicitExt Character. File extension for output. Default "dyp".
#' @param outputDirectory Character. Output directory path. Default NULL (use input directory).
#'
#' @return
#' If \code{toFile=FALSE} (default): For single file, a named list with 193 prosodic features.
#'   For multiple files, named list where each element is a file's feature list.
#'   If \code{toFile=TRUE}, invisibly returns the path(s) to the written JSTF file(s).
#' Feature categories:
#' \describe{
#'   \item{Prosodic metadata}{Duration, PitchKey, PitchRange, PitchMean, IntsIntLabels, UniqueIntsInt}
#'   \item{Spectral features}{L2L1, L2cL1c, L1cLF3c, L1LF3, SLF, C1, SpectralBalance, SLF6D coefficients}
#'   \item{Statistical summaries}{_mean, _std, _var, _iqr, _max, _min for all time-varying features}
#'   \item{Differential features}{_diff versions showing inter-INTSINT-label changes}
#' }
#'
#' @section Time Windowing:
#' Time windowing is handled by the \code{av} package during audio loading,
#' ensuring efficient in-memory processing. Audio is converted directly to
#' parselmouth Sound objects via \code{\link{av_load_for_parselmouth}},
#' eliminating the need for temporary files.
#'
#' @examples
#' \dontrun{
#' # Check if dysprosody is available
#' if (!dysprosody_available()) {
#'   install_dysprosody()
#' }
#'
#' # Single file analysis
#' features <- lst_dysprosody("speech.wav")
#' print(names(features))  # Show all 193 features
#' print(features$Duration)
#' print(features$PitchMean)
#' print(features$IntsIntLabels)
#'
#' # Analyze specific time window
#' features <- lst_dysprosody("speech.wav", beginTime = 1.0, endTime = 3.0)
#'
#' # Batch processing (parallel)
#' files <- list.files("audio_dir", pattern = "\\.wav$", full.names = TRUE)
#' results <- lst_dysprosody(files, verbose = TRUE, parallel = TRUE, n_cores = 4)
#'
#' # Convert to data frame (rows = files, columns = features)
#' library(tidyverse)
#' df <- results %>%
#'   map_dfr(~ as.data.frame(t(unlist(.))), .id = "file")
#'
#' # Write results to JSTF file
#' lst_dysprosody("speech.wav", toFile = TRUE)  # Creates speech.dyp
#'
#' # Read back and convert to data.frame
#' track <- read_track("speech.dyp")
#' df <- as.data.frame(track)
#' head(df)  # Shows begin_time, end_time, and all 193 prosodic features
#' }
#'
#' @references
#' \insertRef{Nylen.2025.10.3389/fnhum.2025.1566274}{superassp}
#'
#' Hirst, D., & Espesser, R. (1993). Automatic Modelling Of Fundamental
#' Frequency Using A Quadratic Spline Function. \emph{Travaux de l'Institut
#' de Phonétique d'Aix}, 15, 71-85.
#'
#' Hirst, D. (2019). INTSINT: a new algorithm using the OMe scale.
#' \emph{ExLing 2018: Proceedings of 9th Tutorial and Research Workshop on
#' Experimental Linguistics}, 53-56. \doi{10.36505/exling-2018/09/0012/000345}
#'
#' @seealso
#' \code{\link{install_dysprosody}}, \code{\link{dysprosody_available}},
#' \code{\link{dysprosody_info}}
#'
#' @export
lst_dysprosody <- function(listOfFiles,
                           beginTime = 0.0,
                           endTime = 0.0,
                           minF = 60,
                           maxF = 750,
                           windowShift = 1.0,
                           verbose = TRUE,
                           parallel = TRUE,
                           n_cores = NULL,
                           toFile = FALSE,
                           explicitExt = "dyp",
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

  # Check dysprosody availability
  if (!dysprosody_available()) {
    cli::cli_abort(c(
      "x" = "Dysprosody module not available",
      "i" = "Install with: install_dysprosody()"
    ))
  }

  # Import dysprosody module
  dysprosody <- reticulate::import("dysprosody", delay_load = FALSE)

  # Normalize time parameters
  beginTime <- if (is.null(beginTime)) rep(0.0, n_files) else beginTime
  endTime <- if (is.null(endTime)) rep(0.0, n_files) else endTime

  if (length(beginTime) == 1 && n_files > 1) {
    beginTime <- rep(beginTime, n_files)
  }
  if (length(endTime) == 1 && n_files > 1) {
    endTime <- rep(endTime, n_files)
  }

  if (length(beginTime) != n_files || length(endTime) != n_files) {
    cli::cli_abort("beginTime and endTime must be same length as listOfFiles")
  }

  # Processing function for a single file
  process_file <- function(i) {
    file_path <- listOfFiles[i]
    bt <- beginTime[i]
    et <- endTime[i]

    tryCatch({
      # Load audio using av package and convert to parselmouth Sound object
      # This uses pure in-memory processing with NO temp files
      sound <- av_load_for_parselmouth(
        file_path = file_path,
        start_time = if (bt > 0) bt else NULL,
        end_time = if (et > 0 && et > bt) et else NULL,
        channels = 1
      )

      # Call Python function with the Sound object directly
      result <- dysprosody$prosody_measures_from_sound(
        sound = sound,
        minF = as.double(minF),
        maxF = as.double(maxF),
        windowShift = as.double(windowShift)
      )

      # Convert pandas Series to R list
      if (!is.null(result)) {
        # Convert to list and ensure numeric types
        features <- as.list(result)
        names(features) <- names(result)
        return(features)
      } else {
        if (verbose) {
          cli::cli_alert_warning("File skipped (< 1 second): {.file {basename(file_path)}}")
        }
        return(NULL)
      }

    }, error = function(e) {
      if (verbose) {
        cli::cli_alert_danger("Error processing {.file {basename(file_path)}}: {e$message}")
      }
      return(NULL)
    })
  }

  # Process files
  if (n_files == 1) {
    # Single file - no parallel
    if (verbose) {
      cli::cli_alert_info("Processing 1 file...")
    }
    result <- process_file(1)


    # Handle JSTF file writing
    if (toFile && !is.null(result)) {
      file_path <- listOfFiles[1]
      bt <- beginTime[1]
      et <- endTime[1]

      # Get audio metadata
      audio_info <- av::av_media_info(file_path)
      sample_rate <- audio_info$audio$sample_rate
      audio_duration <- audio_info$duration

      json_obj <- create_json_track_obj(
        results = result,
        function_name = "lst_dysprosody",
        file_path = file_path,
        sample_rate = sample_rate,
        audio_duration = audio_duration,
        beginTime = bt,
        endTime = if (et > 0) et else audio_duration,
        parameters = list(
          minF = minF,
          maxF = maxF,
          windowShift = windowShift
        )
      )

      base_name <- tools::file_path_sans_ext(basename(file_path))
      out_dir <- if (is.null(outputDirectory)) dirname(file_path) else outputDirectory
      output_path <- file.path(out_dir, paste0(base_name, ".", explicitExt))

      write_json_track(json_obj, output_path)
      return(invisible(output_path))
    }

    return(result)

  } else {
    # Multiple files - optionally parallel
    if (verbose) {
      cli::cli_alert_info("Processing {n_files} files...")
    }

    if (parallel && n_files > 1) {
      # Parallel processing - NOTE: Python dysprosody module batch_process
      # expects file paths, so we cannot use pure in-memory approach for
      # batch operations. For now, fall back to sequential processing.
      # TODO: Modify Python dysprosody to accept Sound objects in batch_process
      if (verbose) {
        cli::cli_alert_warning("Parallel processing not yet supported with in-memory approach")
        cli::cli_alert_info("Using sequential processing...")
      }

      # Fall through to sequential processing below
      parallel <- FALSE
    }

    if (!parallel) {
      # Sequential processing
      if (verbose) {
        cli::cli_progress_bar("Processing files", total = n_files)
      }

      results <- list()
      for (i in seq_len(n_files)) {
        result <- process_file(i)
        file_basename <- tools::file_path_sans_ext(basename(listOfFiles[i]))
        results[[file_basename]] <- result

        if (verbose) {
          cli::cli_progress_update()
        }
      }

      if (verbose) {
        cli::cli_progress_done()
      }
    }

    # Report success
    n_success <- sum(sapply(results, function(x) !is.null(x)))
    if (verbose) {
      cli::cli_alert_success("Processed {n_success}/{n_files} files successfully")
    }


    # Handle JSTF file writing for multiple files
    if (toFile) {
      output_paths <- character(n_files)
      for (i in seq_len(n_files)) {
        file_basename <- tools::file_path_sans_ext(basename(listOfFiles[i]))
        result <- results[[file_basename]]

        if (!is.null(result)) {
          file_path <- listOfFiles[i]
          bt <- beginTime[i]
          et <- endTime[i]

          # Get audio metadata
          audio_info <- av::av_media_info(file_path)
          sample_rate <- audio_info$audio$sample_rate
          audio_duration <- audio_info$duration

          json_obj <- create_json_track_obj(
            results = result,
            function_name = "lst_dysprosody",
            file_path = file_path,
            sample_rate = sample_rate,
            audio_duration = audio_duration,
            beginTime = bt,
            endTime = if (et > 0) et else audio_duration,
            parameters = list(
              minF = minF,
              maxF = maxF,
              windowShift = windowShift
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


    return(results)
  }
}


# Set function attributes for consistency with other DSP functions
attr(lst_dysprosody, "type") <- "list"
attr(lst_dysprosody, "module") <- "dysprosody"
attr(lst_dysprosody, "ext") <- "dyp"
attr(lst_dysprosody, "outputType") <- "JSTF"
attr(lst_dysprosody, "format") <- "JSON"

