#' JSTF Helper Functions
#'
#' Internal helper functions for writing LST function results to JSTF files.
#' These helpers reduce code duplication across lst_* functions.
#'
#' @name jstf_helpers
#' @keywords internal
NULL

#' Write LST Results to JSTF Files
#'
#' Helper function to write list-based DSP results to JSTF format files.
#' This function standardizes the file writing logic used across all lst_* functions.
#'
#' @param results List or list of lists containing DSP results
#' @param file_paths Character vector of input file paths
#' @param beginTime Numeric vector of start times (seconds)
#' @param endTime Numeric vector of end times (seconds, 0 = full file)
#' @param function_name Character string identifying the calling function
#' @param parameters Named list of DSP parameters to store in metadata
#' @param explicitExt Character string for output file extension
#' @param outputDirectory Character string for output directory (NULL = same as input)
#' @param speaker_id Optional character vector of speaker IDs for output filenames
#' @param verbose Logical, show progress messages (default FALSE)
#'
#' @return Character vector of output file paths (single path for single file input)
#'
#' @details
#' This function handles:
#' - Audio metadata extraction via av package
#' - Time range calculation
#' - JSTF object creation
#' - Output filename determination (speaker ID or file basename)
#' - File writing with error handling
#' - Single vs multiple file return value formatting
#'
#' Results with errors (result$error not NULL) are written as NA_character_.
#'
#' @keywords internal
#' @noRd
write_lst_results_to_jstf <- function(results,
                                      file_paths,
                                      beginTime,
                                      endTime,
                                      function_name,
                                      parameters,
                                      explicitExt,
                                      outputDirectory = NULL,
                                      speaker_id = NULL,
                                      verbose = FALSE) {

  # Ensure results is a list
  if (!is.list(results)) {
    stop("results must be a list or list of lists", call. = FALSE)
  }

  # Handle single result wrapped in list vs list of results
  n_files <- length(file_paths)

  # Detect structure: is results a single result (named list) or list of results?
  if (n_files == 1) {
    # Single file case: check if results is already wrapped
    if (is.null(names(results)) && length(results) == 1 && is.list(results[[1]])) {
      # Already wrapped: list(list(...))
      results_list <- results
    } else if (!is.null(names(results))) {
      # Named list: list(field1=..., field2=...)
      results_list <- list(results)
    } else {
      # Assume already wrapped
      results_list <- results
    }
  } else {
    # Multiple files case
    results_list <- results
  }

  # Ensure time vectors match file count
  if (length(beginTime) == 1 && n_files > 1) {
    beginTime <- rep(beginTime, n_files)
  }
  if (length(endTime) == 1 && n_files > 1) {
    endTime <- rep(endTime, n_files)
  }

  # Initialize output paths
  output_paths <- character(n_files)

  # Process each file
  for (i in seq_len(n_files)) {
    result <- results_list[[i]]
    file_path <- file_paths[i]

    # Skip NULL or error results
    if (is.null(result) || (!is.null(result$error))) {
      output_paths[i] <- NA_character_
      if (verbose) {
        warning("Skipping JSTF write for ", basename(file_path),
                " (NULL or error result)", call. = FALSE)
      }
      next
    }

    tryCatch({
      # Get audio metadata
      audio_info <- av::av_media_info(file_path)
      sample_rate <- audio_info$audio$sample_rate
      audio_duration <- audio_info$duration

      # Calculate analysis time range
      analysis_begin <- beginTime[i]
      analysis_end <- if (endTime[i] > 0) endTime[i] else audio_duration

      # Create JSTF object
      json_obj <- create_json_track_obj(
        results = result,
        function_name = function_name,
        file_path = file_path,
        sample_rate = sample_rate,
        audio_duration = audio_duration,
        beginTime = analysis_begin,
        endTime = analysis_end,
        parameters = parameters
      )

      # Determine output filename
      if (!is.null(speaker_id) && length(speaker_id) >= i && !is.null(speaker_id[i]) && nchar(speaker_id[i]) > 0) {
        base_name <- as.character(speaker_id[i])
      } else {
        base_name <- tools::file_path_sans_ext(basename(file_path))
      }

      # Determine output directory
      out_dir <- if (is.null(outputDirectory)) dirname(file_path) else outputDirectory
      output_path <- file.path(out_dir, paste0(base_name, ".", explicitExt))

      # Write JSTF file
      write_json_track(json_obj, output_path)
      output_paths[i] <- output_path

      if (verbose) {
        message("Wrote JSTF: ", output_path)
      }

    }, error = function(e) {
      warning("Failed to write JSTF for ", basename(file_path), ": ",
              e$message, call. = FALSE)
      output_paths[i] <- NA_character_
    })
  }

  # Return single path for single file, vector for multiple files
  if (n_files == 1) {
    return(output_paths[1])
  } else {
    return(output_paths)
  }
}
