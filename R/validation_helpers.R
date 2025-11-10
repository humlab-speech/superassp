#' Parameter Validation Helpers
#'
#' @description
#' This file provides validation helper functions used throughout superassp
#' for consistent parameter checking and error messaging.
#'
#' @name validation_helpers
#' @keywords internal
NULL

#' Validate JSTF output parameters
#'
#' Validates the standard JSTF file output parameters (`toFile`, `explicitExt`,
#' `outputDirectory`) used across all DSP functions that support JSTF output.
#'
#' @param toFile Logical value indicating whether to write to file
#' @param explicitExt Character string specifying the file extension
#' @param outputDirectory Character string or NULL specifying output directory
#' @param function_name Character string with the calling function name (for error messages)
#'
#' @return Invisible TRUE if all validations pass
#' @keywords internal
#' @noRd
#'
#' @examples
#' \dontrun{
#' # In a DSP function
#' validate_jstf_parameters(toFile, explicitExt, outputDirectory, "lst_covarep_vq")
#' }
validate_jstf_parameters <- function(toFile, explicitExt, outputDirectory,
                                     function_name = "unknown function") {

  # Validate toFile
  if (!is.logical(toFile)) {
    stop(
      function_name, "(): toFile must be logical (TRUE/FALSE), not ",
      class(toFile)[1],
      call. = FALSE
    )
  }

  if (length(toFile) != 1) {
    stop(
      function_name, "(): toFile must be a single logical value, not a vector of length ",
      length(toFile),
      call. = FALSE
    )
  }

  if (is.na(toFile)) {
    stop(
      function_name, "(): toFile cannot be NA",
      call. = FALSE
    )
  }

  # Validate explicitExt
  if (!is.character(explicitExt)) {
    stop(
      function_name, "(): explicitExt must be a character string, not ",
      class(explicitExt)[1],
      call. = FALSE
    )
  }

  if (length(explicitExt) != 1) {
    stop(
      function_name, "(): explicitExt must be a single character string, not a vector of length ",
      length(explicitExt),
      call. = FALSE
    )
  }

  if (is.na(explicitExt)) {
    stop(
      function_name, "(): explicitExt cannot be NA",
      call. = FALSE
    )
  }

  if (nchar(explicitExt) == 0) {
    stop(
      function_name, "(): explicitExt cannot be empty",
      call. = FALSE
    )
  }

  # Validate explicitExt format (alphanumeric + underscore + hyphen only)
  if (!grepl("^[a-zA-Z0-9_-]+$", explicitExt)) {
    stop(
      function_name, "(): explicitExt must contain only letters, numbers, hyphens, and underscores. Got: '",
      explicitExt, "'",
      call. = FALSE
    )
  }

  # Validate outputDirectory
  if (!is.null(outputDirectory)) {
    if (!is.character(outputDirectory)) {
      stop(
        function_name, "(): outputDirectory must be a character string or NULL, not ",
        class(outputDirectory)[1],
        call. = FALSE
      )
    }

    if (length(outputDirectory) != 1) {
      stop(
        function_name, "(): outputDirectory must be a single character string, not a vector of length ",
        length(outputDirectory),
        call. = FALSE
      )
    }

    if (is.na(outputDirectory)) {
      stop(
        function_name, "(): outputDirectory cannot be NA (use NULL for default behavior)",
        call. = FALSE
      )
    }

    if (nchar(outputDirectory) == 0) {
      stop(
        function_name, "(): outputDirectory cannot be empty (use NULL for default behavior)",
        call. = FALSE
      )
    }

    # Check directory exists
    if (!dir.exists(outputDirectory)) {
      stop(
        function_name, "(): outputDirectory does not exist: ", outputDirectory,
        "\nCreate the directory first or use NULL to save in the same directory as input files.",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

#' Validate time window parameters
#'
#' Validates `beginTime` and `endTime` parameters used for time windowing
#' in DSP functions.
#'
#' @param beginTime Numeric value or vector for start time(s) in seconds
#' @param endTime Numeric value or vector for end time(s) in seconds
#' @param n_files Integer number of files being processed
#' @param function_name Character string with the calling function name (for error messages)
#'
#' @return Invisible TRUE if all validations pass
#' @keywords internal
#' @noRd
validate_time_window <- function(beginTime, endTime, n_files,
                                  function_name = "unknown function") {

  # Validate beginTime
  if (!is.numeric(beginTime)) {
    stop(
      function_name, "(): beginTime must be numeric, not ",
      class(beginTime)[1],
      call. = FALSE
    )
  }

  if (any(is.na(beginTime))) {
    stop(
      function_name, "(): beginTime cannot contain NA values",
      call. = FALSE
    )
  }

  if (any(beginTime < 0)) {
    stop(
      function_name, "(): beginTime cannot be negative. Got: ",
      paste(beginTime[beginTime < 0], collapse = ", "),
      call. = FALSE
    )
  }

  # Validate endTime
  if (!is.numeric(endTime)) {
    stop(
      function_name, "(): endTime must be numeric, not ",
      class(endTime)[1],
      call. = FALSE
    )
  }

  if (any(is.na(endTime))) {
    stop(
      function_name, "(): endTime cannot contain NA values",
      call. = FALSE
    )
  }

  if (any(endTime < 0)) {
    stop(
      function_name, "(): endTime cannot be negative. Got: ",
      paste(endTime[endTime < 0], collapse = ", "),
      call. = FALSE
    )
  }

  # Validate length consistency
  len_begin <- length(beginTime)
  len_end <- length(endTime)

  if (len_begin != 1 && len_begin != n_files) {
    stop(
      function_name, "(): beginTime must be either length 1 or length ", n_files,
      " (matching number of files). Got length ", len_begin,
      call. = FALSE
    )
  }

  if (len_end != 1 && len_end != n_files) {
    stop(
      function_name, "(): endTime must be either length 1 or length ", n_files,
      " (matching number of files). Got length ", len_end,
      call. = FALSE
    )
  }

  # Validate time ranges (endTime >= beginTime when both specified)
  # Note: endTime = 0 means "use full file duration"
  if (len_begin == len_end) {
    invalid_ranges <- (endTime > 0) & (endTime <= beginTime)
    if (any(invalid_ranges)) {
      stop(
        function_name, "(): endTime must be greater than beginTime (or 0 for full duration). ",
        "Invalid ranges at indices: ", paste(which(invalid_ranges), collapse = ", "),
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

#' Validate file paths
#'
#' Validates that file paths exist and are readable.
#'
#' @param file_paths Character vector of file paths
#' @param function_name Character string with the calling function name (for error messages)
#' @param allow_missing Logical indicating whether missing files should produce warning instead of error
#'
#' @return Invisible TRUE if all validations pass (or warning if allow_missing=TRUE)
#' @keywords internal
#' @noRd
validate_file_paths <- function(file_paths, function_name = "unknown function",
                                allow_missing = FALSE) {

  if (!is.character(file_paths)) {
    stop(
      function_name, "(): file paths must be character strings, not ",
      class(file_paths)[1],
      call. = FALSE
    )
  }

  if (length(file_paths) == 0) {
    stop(
      function_name, "(): no files provided",
      call. = FALSE
    )
  }

  if (any(is.na(file_paths))) {
    stop(
      function_name, "(): file paths cannot contain NA values",
      call. = FALSE
    )
  }

  # Check file existence
  files_exist <- file.exists(file_paths)

  if (!all(files_exist)) {
    missing_files <- file_paths[!files_exist]

    msg <- paste0(
      function_name, "(): unable to find ", length(missing_files), " file(s):\n",
      paste("  -", missing_files, collapse = "\n")
    )

    if (allow_missing) {
      warning(msg, call. = FALSE)
    } else {
      stop(msg, call. = FALSE)
    }
  }

  invisible(TRUE)
}

#' Validate sample rate
#'
#' Validates that a sample rate is positive and reasonable.
#'
#' @param sample_rate Numeric sample rate in Hz
#' @param function_name Character string with the calling function name (for error messages)
#' @param min_rate Minimum acceptable sample rate (default: 1000 Hz)
#' @param max_rate Maximum acceptable sample rate (default: 192000 Hz)
#'
#' @return Invisible TRUE if validation passes
#' @keywords internal
#' @noRd
validate_sample_rate <- function(sample_rate, function_name = "unknown function",
                                 min_rate = 1000, max_rate = 192000) {

  if (!is.numeric(sample_rate)) {
    stop(
      function_name, "(): sample_rate must be numeric, not ",
      class(sample_rate)[1],
      call. = FALSE
    )
  }

  if (length(sample_rate) != 1) {
    stop(
      function_name, "(): sample_rate must be a single value, not a vector of length ",
      length(sample_rate),
      call. = FALSE
    )
  }

  if (is.na(sample_rate)) {
    stop(
      function_name, "(): sample_rate cannot be NA",
      call. = FALSE
    )
  }

  if (sample_rate <= 0) {
    stop(
      function_name, "(): sample_rate must be positive. Got: ", sample_rate,
      call. = FALSE
    )
  }

  if (sample_rate < min_rate) {
    warning(
      function_name, "(): sample_rate ", sample_rate, " Hz is unusually low (< ",
      min_rate, " Hz). This may produce unexpected results.",
      call. = FALSE
    )
  }

  if (sample_rate > max_rate) {
    warning(
      function_name, "(): sample_rate ", sample_rate, " Hz is unusually high (> ",
      max_rate, " Hz). This may produce unexpected results.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}
