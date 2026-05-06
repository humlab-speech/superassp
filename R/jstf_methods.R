#' JSON Track Conversion Methods
#' 
#' Methods for converting JsonTrackObj to data.frame and tibble formats,
#' making them compatible with AsspDataObj workflows.
#'
#' @name json_track_methods
NULL

#' @describeIn JsonTrackObj Convert to a data.frame; each slice becomes a row.
#' @param x JsonTrackObj
#' @param row.names NULL or character vector of row names
#' @param optional Logical, if TRUE column names are checked for syntactic validity
#' @param ... Additional arguments (ignored)
#' @export
as.data.frame.JsonTrackObj <- function(x, row.names = NULL, optional = FALSE, ...) {
  
  # Validate input
  validate_json_track(x)
  
  # Extract field names
  field_names <- names(x$field_schema)
  n_fields <- length(field_names)
  n_slices <- length(x$slices)
  
  # Initialize data.frame
  df <- data.frame(
    begin_time = numeric(n_slices),
    end_time = numeric(n_slices),
    stringsAsFactors = FALSE
  )
  
  # Add columns for each field
  for (field in field_names) {
    df[[field]] <- vector("list", n_slices)
  }
  
  # Fill in data
  for (i in seq_along(x$slices)) {
    slice <- x$slices[[i]]
    df$begin_time[i] <- slice$begin_time
    df$end_time[i] <- slice$end_time
    
    # Fill field values (access by name, not index)
    for (field in field_names) {
      value <- slice$values[[field]]
      
      # Get field type (field_schema is now a list)
      field_type <- x$field_schema[[field]]
      
      # Handle different types
      if (field_type == "numeric" || 
          field_type == "integer" ||
          field_type == "character" ||
          field_type == "logical") {
        df[[field]][i] <- value
      } else {
        # Keep as list for complex types
        df[[field]][i] <- list(value)
      }
    }
  }
  
  # Convert list columns to vectors where appropriate
  for (field in field_names) {
    field_type <- x$field_schema[[field]]
    if (field_type %in% c("numeric", "integer", "character", "logical")) {
      df[[field]] <- unlist(df[[field]])
    }
  }
  
  # Add metadata as attributes
  attr(df, "function") <- x$function_name
  attr(df, "file_path") <- x$file_path
  attr(df, "sample_rate") <- x$sample_rate
  attr(df, "audio_duration") <- x$audio_duration
  attr(df, "jstf_version") <- x$version
  
  # Set row names if provided
  if (!is.null(row.names)) {
    rownames(df) <- row.names
  }
  
  return(df)
}

#' @describeIn JsonTrackObj Convert to a tibble with typed columns. Requires the tibble package.
#' @export
as_tibble.JsonTrackObj <- function(x, ...) {
  
  # Check if tibble is available
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("tibble package required. Install with: install.packages('tibble')")
  }
  
  # Convert to data.frame first
  df <- as.data.frame(x)
  
  # Convert to tibble
  tbl <- tibble::as_tibble(df, ...)
  
  # Preserve metadata attributes
  attr(tbl, "function") <- attr(df, "function")
  attr(tbl, "file_path") <- attr(df, "file_path")
  attr(tbl, "sample_rate") <- attr(df, "sample_rate")
  attr(tbl, "audio_duration") <- attr(df, "audio_duration")
  attr(tbl, "jstf_version") <- attr(df, "jstf_version")
  
  return(tbl)
}

#' Subset JsonTrackObj
#'
#' Extract slices within a time range
#'
#' @param x JsonTrackObj
#' @param start_time Start time in seconds
#' @param end_time End time in seconds
#'
#' @return JsonTrackObj with filtered slices
subset_json_track <- function(x, start_time = NULL, end_time = NULL) {

  stopifnot(inherits(x, "JsonTrackObj"))

  # Filter slices - keep only slices fully within the range
  keep <- rep(TRUE, length(x$slices))

  if (!is.null(start_time)) {
    keep <- keep & sapply(x$slices, function(s) s$begin_time >= start_time)
  }

  if (!is.null(end_time)) {
    keep <- keep & sapply(x$slices, function(s) s$end_time <= end_time)
  }

  x$slices <- x$slices[keep]

  return(x)
}

#' Merge multiple JsonTrackObj files
#'
#' Combines multiple JsonTrackObj objects from the same audio file.
#' Useful when processing different time slices separately.
#'
#' @param ... JsonTrackObj objects to merge
#'
#' @return Merged JsonTrackObj
merge_json_tracks <- function(...) {
  
  objs <- list(...)
  
  if (length(objs) == 0) {
    stop("No objects provided")
  }
  
  # Check all are JsonTrackObj
  if (!all(sapply(objs, inherits, "JsonTrackObj"))) {
    stop("All inputs must be JsonTrackObj")
  }
  
  # Check all from same function
  funcs <- sapply(objs, function(x) x$function_name)
  if (length(unique(funcs)) > 1) {
    stop("Cannot merge tracks from different functions: ", 
         paste(unique(funcs), collapse = ", "))
  }
  
  # Use first object as template
  merged <- objs[[1]]
  
  # Append slices from other objects
  for (i in 2:length(objs)) {
    merged$slices <- c(merged$slices, objs[[i]]$slices)
  }
  
  # Sort slices by begin_time
  merged$slices <- merged$slices[order(sapply(merged$slices, function(s) s$begin_time))]
  
  return(merged)
}

#' @describeIn JsonTrackObj Print a summary of a JsonTrackObj to the console.
#' @param object JsonTrackObj
#' @export
summary.JsonTrackObj <- function(object, ...) {
  
  cat("JSON Track Object Summary\n")
  cat("=========================\n\n")
  
  cat("Format:       ", object$format, " v", object$version, "\n", sep = "")
  cat("Function:     ", object$function_name, "\n", sep = "")
  cat("File:         ", object$file_path, "\n", sep = "")
  cat("Created:      ", object$created, "\n", sep = "")
  
  if (!is.null(object$sample_rate)) {
    cat("Sample rate:  ", object$sample_rate, " Hz\n", sep = "")
  }
  
  if (!is.null(object$audio_duration)) {
    cat("Duration:     ", object$audio_duration, " s\n", sep = "")
  }
  
  cat("\nFields (", length(object$field_schema), "):\n", sep = "")
  for (i in seq_along(object$field_schema)) {
    cat(sprintf("  %-30s %s\n", 
                names(object$field_schema)[i], 
                object$field_schema[[i]]))
  }
  
  cat("\nSlices (", length(object$slices), "):\n", sep = "")
  if (length(object$slices) > 0) {
    total_duration <- sum(sapply(object$slices, function(s) s$end_time - s$begin_time))
    cat("  Total duration: ", total_duration, " s\n", sep = "")
    cat("  Time range:     [", object$slices[[1]]$begin_time, ", ", 
        object$slices[[length(object$slices)]]$end_time, "]\n", sep = "")
  }
  
  invisible(object)
}
