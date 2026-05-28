#' JSTF Helper Functions
#'
#' Internal helper functions for writing LST function results to JSTF files.
#' These helpers reduce code duplication across lst_* functions.
#'
#' @name jstf_helpers
#' @keywords internal
NULL

#' Build JSTF Objects From LST Results
#'
#' Creates a list of JsonTrackObj objects from DSP results without writing to
#' disk. Use this when you need the objects in memory (e.g., `return_jstf=TRUE`).
#' `write_lst_results_to_jstf()` calls this internally.
#'
#' @param results List or list of lists containing DSP results
#' @param file_paths Character vector of input file paths
#' @param beginTime Numeric vector of start times (seconds)
#' @param endTime Numeric vector of end times (seconds, 0 = full file)
#' @param function_name Character string identifying the calling function
#' @param parameters Named list of DSP parameters to store in metadata
#' @param verbose Logical, emit warnings on per-file errors (default FALSE)
#'
#' @return List of JsonTrackObj (one per file; NULL entries on error)
#'
#' @keywords internal
#' @noRd
build_lst_jstf_objects <- function(results,
                                   file_paths,
                                   beginTime,
                                   endTime,
                                   function_name,
                                   parameters,
                                   verbose = FALSE) {

  if (!is.list(results)) {
    stop("results must be a list or list of lists", call. = FALSE)
  }

  n_files <- length(file_paths)

  # Normalise to list-of-results
  if (n_files == 1) {
    if (is.null(names(results)) && length(results) == 1 && is.list(results[[1]])) {
      results_list <- results
    } else if (!is.null(names(results))) {
      results_list <- list(results)
    } else {
      results_list <- results
    }
  } else {
    results_list <- results
  }

  if (length(beginTime) == 1 && n_files > 1) beginTime <- rep(beginTime, n_files)
  if (length(endTime)   == 1 && n_files > 1) endTime   <- rep(endTime,   n_files)

  jstf_objs <- vector("list", n_files)

  for (i in seq_len(n_files)) {
    result    <- results_list[[i]]
    file_path <- file_paths[i]

    if (is.null(result) || (!is.null(result$error))) {
      jstf_objs[[i]] <- NULL
      next
    }

    tryCatch({
      audio_info     <- media_info(file_path)
      sr             <- audio_info$audio$sample_rate
      audio_duration <- audio_info$duration
      analysis_end   <- if (endTime[i] > 0) endTime[i] else audio_duration

      jstf_objs[[i]] <- create_json_track_obj(
        results        = result,
        function_name  = function_name,
        file_path      = file_path,
        sample_rate    = sr,
        audio_duration = audio_duration,
        beginTime      = beginTime[i],
        endTime        = analysis_end,
        parameters     = parameters
      )
    }, error = function(e) {
      if (verbose) {
        warning("Failed to build JSTF for ", basename(file_path),
                ": ", e$message, call. = FALSE)
      }
      jstf_objs[[i]] <<- NULL
    })
  }

  jstf_objs
}


#' Write LST Results to JSTF Files
#'
#' Helper function to write list-based DSP results to JSTF format files.
#' This function standardizes the file writing logic used across all lst_*
#' functions. Internally calls `build_lst_jstf_objects()` then writes each
#' object to disk.
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

  n_files   <- length(file_paths)
  jstf_objs <- build_lst_jstf_objects(
    results = results, file_paths = file_paths,
    beginTime = beginTime, endTime = endTime,
    function_name = function_name, parameters = parameters,
    verbose = verbose
  )

  output_paths <- character(n_files)

  for (i in seq_len(n_files)) {
    obj       <- jstf_objs[[i]]
    file_path <- file_paths[i]

    if (is.null(obj)) {
      output_paths[i] <- NA_character_
      next
    }

    tryCatch({
      base_name <- if (!is.null(speaker_id) && length(speaker_id) >= i &&
                       !is.null(speaker_id[i]) && nchar(speaker_id[i]) > 0) {
        as.character(speaker_id[i])
      } else {
        tools::file_path_sans_ext(basename(file_path))
      }
      out_dir     <- if (is.null(outputDirectory)) dirname(file_path) else outputDirectory
      output_path <- file.path(out_dir, paste0(base_name, ".", explicitExt))
      write_jstf(obj, output_path)
      output_paths[i] <- output_path
      if (verbose) message("Wrote JSTF: ", output_path)
    }, error = function(e) {
      warning("Failed to write JSTF for ", basename(file_path), ": ",
              e$message, call. = FALSE)
      output_paths[i] <<- NA_character_
    })
  }

  if (n_files == 1) output_paths[1] else output_paths
}
