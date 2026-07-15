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
    cli::cli_abort("{function_name}(): toFile must be logical (TRUE/FALSE), not {.cls {class(toFile)[1]}}")
  }

  if (length(toFile) != 1) {
    cli::cli_abort("{function_name}(): toFile must be a single logical value, not a vector of length {length(toFile)}")
  }

  if (is.na(toFile)) {
    cli::cli_abort("{function_name}(): toFile cannot be NA")
  }

  # Validate explicitExt
  if (!is.character(explicitExt)) {
    cli::cli_abort("{function_name}(): explicitExt must be a character string, not {.cls {class(explicitExt)[1]}}")
  }

  if (length(explicitExt) != 1) {
    cli::cli_abort("{function_name}(): explicitExt must be a single character string, not a vector of length {length(explicitExt)}")
  }

  if (is.na(explicitExt)) {
    cli::cli_abort("{function_name}(): explicitExt cannot be NA")
  }

  if (nchar(explicitExt) == 0) {
    cli::cli_abort("{function_name}(): explicitExt cannot be empty")
  }

  # Validate explicitExt format (alphanumeric + underscore + hyphen only)
  if (!grepl("^[a-zA-Z0-9_-]+$", explicitExt)) {
    cli::cli_abort("{function_name}(): explicitExt must contain only letters, numbers, hyphens, and underscores. Got: {.val {explicitExt}}")
  }

  # Validate outputDirectory
  if (!is.null(outputDirectory)) {
    if (!is.character(outputDirectory)) {
      cli::cli_abort("{function_name}(): outputDirectory must be a character string or NULL, not {.cls {class(outputDirectory)[1]}}")
    }

    if (length(outputDirectory) != 1) {
      cli::cli_abort("{function_name}(): outputDirectory must be a single character string, not a vector of length {length(outputDirectory)}")
    }

    if (is.na(outputDirectory)) {
      cli::cli_abort("{function_name}(): outputDirectory cannot be NA (use NULL for default behavior)")
    }

    if (nchar(outputDirectory) == 0) {
      cli::cli_abort("{function_name}(): outputDirectory cannot be empty (use NULL for default behavior)")
    }

    # Check directory exists
    if (!dir.exists(outputDirectory)) {
      cli::cli_abort(c(
        "{function_name}(): outputDirectory does not exist: {.path {outputDirectory}}",
        "i" = "Create the directory first or use NULL to save in the same directory as input files."
      ))
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
    cli::cli_abort("{function_name}(): beginTime must be numeric, not {.cls {class(beginTime)[1]}}")
  }

  if (any(is.na(beginTime))) {
    cli::cli_abort("{function_name}(): beginTime cannot contain NA values")
  }

  if (any(beginTime < 0)) {
    cli::cli_abort("{function_name}(): beginTime cannot be negative. Got: {.val {beginTime[beginTime < 0]}}")
  }

  # Validate endTime
  if (!is.numeric(endTime)) {
    cli::cli_abort("{function_name}(): endTime must be numeric, not {.cls {class(endTime)[1]}}")
  }

  if (any(is.na(endTime))) {
    cli::cli_abort("{function_name}(): endTime cannot contain NA values")
  }

  if (any(endTime < 0)) {
    cli::cli_abort("{function_name}(): endTime cannot be negative. Got: {.val {endTime[endTime < 0]}}")
  }

  # Validate length consistency
  len_begin <- length(beginTime)
  len_end <- length(endTime)

  if (len_begin != 1 && len_begin != n_files) {
    cli::cli_abort("{function_name}(): beginTime must be either length 1 or length {n_files} (matching number of files). Got length {len_begin}")
  }

  if (len_end != 1 && len_end != n_files) {
    cli::cli_abort("{function_name}(): endTime must be either length 1 or length {n_files} (matching number of files). Got length {len_end}")
  }

  # Validate time ranges (endTime >= beginTime when both specified)
  # Note: endTime = 0 means "use full file duration"
  if (len_begin == len_end) {
    invalid_ranges <- (endTime > 0) & (endTime <= beginTime)
    if (any(invalid_ranges)) {
      cli::cli_abort("{function_name}(): endTime must be greater than beginTime (or 0 for full duration). Invalid ranges at indices: {.val {which(invalid_ranges)}}")
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
    cli::cli_abort("{function_name}(): file paths must be character strings, not {.cls {class(file_paths)[1]}}")
  }

  if (length(file_paths) == 0) {
    cli::cli_abort("{function_name}(): no files provided")
  }

  if (any(is.na(file_paths))) {
    cli::cli_abort("{function_name}(): file paths cannot contain NA values")
  }

  # Check file existence
  files_exist <- file.exists(file_paths)

  if (!all(files_exist)) {
    missing_files <- file_paths[!files_exist]
    bullets <- stats::setNames(missing_files, rep("*", length(missing_files)))
    msg <- c(
      "{function_name}(): unable to find {length(missing_files)} file{?s}:",
      bullets
    )
    if (allow_missing) {
      cli::cli_warn(msg)
    } else {
      cli::cli_abort(msg)
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
    cli::cli_abort("{function_name}(): sample_rate must be numeric, not {.cls {class(sample_rate)[1]}}")
  }

  if (length(sample_rate) != 1) {
    cli::cli_abort("{function_name}(): sample_rate must be a single value, not a vector of length {length(sample_rate)}")
  }

  if (is.na(sample_rate)) {
    cli::cli_abort("{function_name}(): sample_rate cannot be NA")
  }

  if (sample_rate <= 0) {
    cli::cli_abort("{function_name}(): sample_rate must be positive. Got: {.val {sample_rate}}")
  }

  if (sample_rate < min_rate) {
    cli::cli_warn("{function_name}(): sample_rate {.val {sample_rate}} Hz is unusually low (< {min_rate} Hz). This may produce unexpected results.")
  }

  if (sample_rate > max_rate) {
    cli::cli_warn("{function_name}(): sample_rate {.val {sample_rate}} Hz is unusually high (> {max_rate} Hz). This may produce unexpected results.")
  }

  invisible(TRUE)
}
