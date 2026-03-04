#' JsonTrackObj — JSON Track Format Object
#'
#' A list-based S3 class representing a JSTF (JSON Speech Track Format) file
#' in memory. Produced by `lst_*` functions with `toFile = FALSE` and read
#' back by `read_json_track()`.
#'
#' @name JsonTrackObj
#' @aliases JsonTrackObj
NULL

#' JSON Track Object Core Functions
#'
#' Core infrastructure for JSON Track Format (JSTF).
#'
#' @name json_track_core
#' @keywords internal
NULL

#' Create a JsonTrackObj
#'
#' @param results List of results from lst_* function
#' @param function_name Name of the DSP function
#' @param file_path Original audio file path
#' @param sample_rate Audio sample rate in Hz
#' @param audio_duration Total audio duration in seconds
#' @param beginTime Start time in seconds
#' @param endTime End time in seconds
#' @param parameters Function parameters used
#'
#' @return Object of class JsonTrackObj
create_json_track_obj <- function(results,
                                  function_name,
                                  file_path,
                                  sample_rate = NULL,
                                  audio_duration = NULL,
                                  beginTime = 0.0,
                                  endTime = 0.0,
                                  parameters = list()) {
  
  # Infer field schema from results
  field_schema <- infer_field_schema(results)
  
  # Convert named vector to list to preserve names in JSON
  field_schema <- as.list(field_schema)
  
  # Create slice
  slice <- list(
    begin_time = beginTime,
    end_time = if (endTime > 0) endTime else audio_duration,
    values = extract_values_from_results(results, names(field_schema))
  )
  
  # Build JsonTrackObj
  obj <- structure(
    list(
      format = "JSTF",
      version = "1.0",
      created = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
      function_name = function_name,
      file_path = file_path,
      sample_rate = sample_rate,
      audio_duration = audio_duration,
      metadata = list(
        function_version = utils::packageVersion("superassp"),
        parameters = parameters
      ),
      field_schema = field_schema,
      slices = list(slice)
    ),
    class = c("JsonTrackObj", "list")
  )
  
  return(obj)
}

#' Infer field schema from results
#'
#' @param results List of results
#' @return Named character vector of field types
#' @keywords internal
infer_field_schema <- function(results) {
  
  # Handle data frame results
  if (is.data.frame(results)) {
    return(sapply(results, function(col) {
      if (is.numeric(col)) "numeric"
      else if (is.integer(col)) "integer"
      else if (is.character(col)) "character"
      else if (is.logical(col)) "logical"
      else "list"
    }))
  }
  
  # Handle list results
  if (is.list(results)) {
    schema <- sapply(results, function(val) {
      if (is.numeric(val) && length(val) == 1) "numeric"
      else if (is.integer(val) && length(val) == 1) "integer"
      else if (is.character(val) && length(val) == 1) "character"
      else if (is.logical(val) && length(val) == 1) "logical"
      else if (is.numeric(val) && length(val) > 1) "numeric_vector"
      else if (is.matrix(val)) "matrix"
      else "list"
    })
    return(schema)
  }
  
  stop("Cannot infer schema from results type: ", class(results)[1])
}

#' Extract values from results matching field schema
#'
#' @param results Results object
#' @param field_names Names of fields to extract
#' @return List of values
#' @keywords internal
extract_values_from_results <- function(results, field_names) {
  
  if (is.data.frame(results)) {
    # Extract first row as list
    return(as.list(results[1, field_names]))
  }
  
  if (is.list(results)) {
    return(results[field_names])
  }
  
  stop("Cannot extract values from results type: ", class(results)[1])
}

#' Append a slice to JsonTrackObj
#'
#' @param obj JsonTrackObj
#' @param results New results to append
#' @param beginTime Start time
#' @param endTime End time
#' @return Updated JsonTrackObj
append_json_track_slice <- function(obj, results, beginTime, endTime) {
  
  stopifnot(inherits(obj, "JsonTrackObj"))
  
  # Extract values
  values <- extract_values_from_results(results, names(obj$field_schema))
  
  # Create new slice
  new_slice <- list(
    begin_time = beginTime,
    end_time = endTime,
    values = values
  )
  
  # Append to slices
  obj$slices <- c(obj$slices, list(new_slice))
  
  return(obj)
}

#' Validate JsonTrackObj
#'
#' @return TRUE if valid, otherwise throws error
validate_json_track <- function(obj) {
  
  # Check class
  if (!inherits(obj, "JsonTrackObj")) {
    stop("Object is not a JsonTrackObj")
  }
  
  # Check format
  if (obj$format != "JSTF") {
    stop("Invalid format: expected 'JSTF', got '", obj$format, "'")
  }
  
  # Check version
  if (!obj$version %in% c("1.0")) {
    warning("Unknown version: ", obj$version)
  }
  
  # Check required fields
  required <- c("format", "version", "function_name", "field_schema", "slices")
  missing <- setdiff(required, names(obj))
  if (length(missing) > 0) {
    stop("Missing required fields: ", paste(missing, collapse = ", "))
  }
  
  # Validate slices
  for (i in seq_along(obj$slices)) {
    slice <- obj$slices[[i]]
    
    # Check time range
    if (slice$begin_time >= slice$end_time) {
      stop("Slice ", i, ": begin_time must be < end_time")
    }
    
    # Check values length
    if (length(slice$values) != length(obj$field_schema)) {
      stop("Slice ", i, ": values length (", length(slice$values),
           ") doesn't match field_schema length (", length(obj$field_schema), ")")
    }
  }
  
  # Check for overlapping slices (warning only)
  if (length(obj$slices) > 1) {
    for (i in 1:(length(obj$slices) - 1)) {
      if (obj$slices[[i]]$end_time > obj$slices[[i + 1]]$begin_time) {
        warning("Overlapping slices detected at index ", i, " and ", i + 1)
      }
    }
  }
  
  return(TRUE)
}

#' @describeIn JsonTrackObj Print a compact summary of a JsonTrackObj.
#' @export
print.JsonTrackObj <- function(x, ...) {
  cat("JsonTrackObj (", x$format, " v", x$version, ")\n", sep = "")
  cat("Function:     ", x$function_name, "\n", sep = "")
  cat("File:         ", basename(x$file_path), "\n", sep = "")
  cat("Duration:     ", x$audio_duration, " s\n", sep = "")
  cat("Sample rate:  ", x$sample_rate, " Hz\n", sep = "")
  cat("Fields:       ", length(x$field_schema), " (", 
      paste(head(names(x$field_schema), 3), collapse = ", "),
      ifelse(length(x$field_schema) > 3, ", ...", ""), ")\n", sep = "")
  cat("Slices:       ", length(x$slices), "\n", sep = "")
  
  if (length(x$slices) > 0) {
    cat("\nTime ranges:\n")
    for (i in seq_len(min(5, length(x$slices)))) {
      cat(sprintf("  [%6.2f, %6.2f]", 
                  x$slices[[i]]$begin_time, 
                  x$slices[[i]]$end_time), "\n")
    }
    if (length(x$slices) > 5) {
      cat("  ... and", length(x$slices) - 5, "more slices\n")
    }
  }
  
  invisible(x)
}
