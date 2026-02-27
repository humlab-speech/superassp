#' AsspDataObj — ASSP Data Object
#'
#' S3 class for in-memory ASSP/SSFF signal data. Produced by `read_ssff()`,
#' `read_audio()`, and all `trk_*` functions with `toFile = FALSE`.
#'
#' @name AsspDataObj
#' @aliases AsspDataObj
NULL

#' @describeIn AsspDataObj Convert to a data.frame with template expansion, optional clean names and unit assignment.
#' @param x AsspDataObj. Object to convert.
#' @param ... Additional arguments (currently unused).
#' @param convert_units Logical. Assign units to columns based on unit suffix. Default: TRUE.
#' @param clean_names Logical. Convert bracket notation to underscore notation. Default: TRUE.
#' @param na.zeros Logical. Convert zeros to NA. Default: FALSE.
#' @export
as.data.frame.AsspDataObj <- function(x, ...,
                                      convert_units = TRUE,
                                      clean_names = TRUE,
                                      na.zeros = FALSE) {

  # Get basic metadata
  sample_rate <- attr(x, "sampleRate")
  start_time <- attr(x, "startTime")
  if (is.null(start_time)) start_time <- 0.0

  # Get track names (template notation)
  # Priority: 1) Object attribute, 2) Registry lookup, 3) Fallback to list names
  track_names <- attr(x, "tracks")

  if (is.null(track_names)) {
    # Try registry lookup if we know the function name
    func_name <- attr(x, "func")  # Some objects may have this
    if (!is.null(func_name) && exists("wrasspOutputInfos")) {
      registry <- get("wrasspOutputInfos")
      if (func_name %in% names(registry)) {
        track_names <- registry[[func_name]]$tracks
      }
    }
  }

  if (is.null(track_names)) {
    track_names <- names(x)  # Fallback to list names
  }

  # Expand templates and collect columns
  expanded_cols <- list()

  for (i in seq_along(names(x))) {
    track_name <- names(x)[i]
    track_data <- x[[track_name]]

    # Get template (may be NULL for unnamed tracks)
    if (i <= length(track_names)) {
      template <- track_names[i]
    } else {
      template <- track_name  # Fallback
    }

    if (is.null(template) || is.na(template)) {
      template <- track_name
    }

    # Check if matrix (multi-column track)
    if (is.matrix(track_data) && ncol(track_data) > 1) {
      # EXPAND TEMPLATE
      col_names <- .expand_track_template(template, ncol(track_data))

      # Split matrix into individual columns
      for (j in seq_len(ncol(track_data))) {
        col_data <- track_data[, j]

        # Convert zeros to NA if requested
        if (na.zeros) {
          col_data[col_data == 0] <- NA
        }

        expanded_cols[[col_names[j]]] <- col_data
      }
    } else {
      # Single column or vector
      col_data <- as.vector(track_data)

      # Convert zeros to NA if requested
      if (na.zeros) {
        col_data[col_data == 0] <- NA
      }

      # For single-column templates, still expand (Fi → F1)
      if (.has_placeholder(template)) {
        col_name <- .expand_track_template(template, 1)[1]
      } else {
        col_name <- template
      }

      expanded_cols[[col_name]] <- col_data
    }
  }

  # Determine number of frames
  if (length(expanded_cols) > 0) {
    n_frames <- length(expanded_cols[[1]])
  } else {
    n_frames <- 0
  }

  # Create time vector
  if (!is.null(sample_rate) && sample_rate > 0) {
    frame_time <- start_time + (seq_len(n_frames) - 1) / sample_rate
  } else {
    # Fallback: assume 1 sample per second
    frame_time <- start_time + seq_len(n_frames) - 1
  }

  # Combine into data frame
  df <- data.frame(
    frame_time = frame_time,
    stringsAsFactors = FALSE
  )

  # Add expanded columns
  for (col_name in names(expanded_cols)) {
    df[[col_name]] <- expanded_cols[[col_name]]
  }

  # Clean names if requested
  if (clean_names) {
    # Store original names before cleaning for attribute mapping
    original_names <- names(df)

    names(df) <- .clean_track_names(names(df))

    # Update expanded_cols keys to match cleaned names
    names_map <- setNames(names(df), original_names)
  }

  # Assign units if requested
  if (convert_units) {
    df <- .assign_track_units(df)
  }

  # Generate and store label mappings as attributes
  attr(df, "track_labels") <- .generate_track_labels(names(df))
  attr(df, "track_descriptions") <- .generate_track_descriptions(names(df))

  # Store original AsspDataObj metadata
  attr(df, "sampleRate") <- sample_rate
  attr(df, "startTime") <- start_time

  df
}



#' Get track label for plotting
#'
#' Extracts the display label for a track column, either from stored attributes
#' or by intelligent inference.
#'
#' @param df data.frame or tibble. Data frame from `as.data.frame.AsspDataObj()`.
#' @param col Character. Column name.
#' @param full Logical. If TRUE, return full descriptive label. If FALSE
#'   (default), return short label suitable for plot axes.
#'
#' @return Character. Display label for the track.
#'
#' @details
#' This function first checks the `track_labels` or `track_descriptions`
#' attributes stored by `as.data.frame.AsspDataObj()`. If not found, it
#' intelligently infers the label from the column name.
#'
#' **Short labels** (full = FALSE):
#' - "fo \[Hz\]", "F1 \[Hz\]", "H1-H2c \[dB\]"
#' - Suitable for plot axes
#'
#' **Full labels** (full = TRUE):
#' - "Frequency of oscillation \[Hz\]"
#' - "First formant frequency \[Hz\]"
#' - "H1-H2 corrected for formants \[dB\]"
#' - Suitable for documentation, papers
#'
#' @examples
#' \dontrun{
#' fms <- wrassp::forest("audio.wav", toFile = FALSE)
#' df <- as.data.frame(fms)
#'
#' # Short label
#' get_track_label(df, "F1_Hz")
#' # [1] "F1 [Hz]"
#'
#' # Full descriptive label
#' get_track_label(df, "F1_Hz", full = TRUE)
#' # [1] "First formant frequency [Hz]"
#'
#' # Works even without attributes (inference)
#' df2 <- data.frame(frame_time = 1:10, fo_Hz = rnorm(10, 120, 10))
#' get_track_label(df2, "fo_Hz")
#' # [1] "fo [Hz]"
#' }
get_track_label <- function(df, col, full = FALSE) {
  labels_attr <- if (full) "track_descriptions" else "track_labels"
  labels <- attr(df, labels_attr)

  if (!is.null(labels) && !is.null(labels[[col]])) {
    return(labels[[col]])
  }

  # Fallback: infer from column name
  .infer_track_label(col, full = full)
}
