# AsspDataObj S3 Methods
#
# S3 methods for converting AsspDataObj to data.frame and tibble with
# automatic track name cleaning, unit assignment, and label generation.

#' Convert AsspDataObj to data.frame
#'
#' Converts an AsspDataObj (from wrassp or superassp DSP functions) to a
#' data.frame with clean column names, optional unit assignment, and automatic
#' label generation for plotting.
#'
#' @param x AsspDataObj. Object to convert.
#' @param ... Additional arguments (currently unused).
#' @param convert_units Logical. If TRUE (default), assign R units package
#'   units to columns based on their unit suffix. Requires 'units' package.
#' @param clean_names Logical. If TRUE (default), convert SSFF-style bracket
#'   notation to R-friendly underscore notation:
#'   - `fo[Hz]` → `fo_Hz`
#'   - `H1-H2c[dB]` → `H1_H2c_dB`
#' @param na.zeros Logical. If TRUE, convert zero values to NA. Useful for
#'   pitch tracks where zeros indicate unvoiced frames. Default: FALSE.
#'
#' @return data.frame with:
#'   - `frame_time` column (time in seconds)
#'   - Track columns with clean names (if `clean_names = TRUE`)
#'   - Units assigned (if `convert_units = TRUE`)
#'   - Attributes: `track_labels`, `track_descriptions` for plotting
#'
#' @details
#' This method implements the three-layer naming strategy:
#'
#' **Layer 1 (SSFF/AsspDataObj)**: Scientific notation
#' - Uses brackets: `fo[Hz]`, `Fi[Hz]`, `Bi[Hz]`
#' - Titze 2015 / Nylén 2024 compliant
#'
#' **Layer 2 (data.frame)**: R-friendly notation
#' - Uses underscores: `fo_Hz`, `F1_Hz`, `B1_Hz`
#' - No backticks needed for column access
#'
#' **Layer 3 (plotting)**: Display labels
#' - Short: "fo [Hz]", "F1 [Hz]"
#' - Full: "Frequency of oscillation [Hz]"
#'
#' **Template Expansion**:
#'
#' For multi-column tracks (formants, LP coefficients), templates with
#' placeholder 'i' are expanded based on actual matrix dimensions:
#' - `Fi[Hz]` → `F1[Hz]`, `F2[Hz]`, `F3[Hz]`, ... (runtime)
#' - `LPCi` → `LPC1`, `LPC2`, `LPC3`, ... (runtime)
#'
#' **Units Assignment**:
#'
#' When `convert_units = TRUE`, columns are assigned units based on suffix:
#' - `_Hz` → Hz
#' - `_dB` → dB (as attribute, not standard unit)
#' - `_pct` → percent
#' - `_us` → microseconds
#' - `_s` → seconds
#'
#' @seealso
#' - [as_tibble.AsspDataObj()] for tibble conversion
#' - [get_track_label()] for extracting plot labels
#' - [ggtrack()] for auto-labeled ggplot2 plots
#'
#' @export
#' @examples
#' \dontrun{
#' # Get formants with default settings
#' fms <- wrassp::forest("audio.wav", numFormants = 4, toFile = FALSE)
#'
#' # Convert to data.frame with clean names and units
#' df <- as.data.frame(fms)
#' names(df)
#' # [1] "frame_time" "F1_Hz" "F2_Hz" "F3_Hz" "F4_Hz"
#' # [6] "B1_Hz" "B2_Hz" "B3_Hz" "B4_Hz"
#'
#' class(df$F1_Hz)  # "units"
#'
#' # Without clean names (keep brackets)
#' df_brackets <- as.data.frame(fms, clean_names = FALSE)
#' names(df_brackets)
#' # [1] "frame_time" "F1[Hz]" "F2[Hz]" ... (requires backticks)
#'
#' # Without units
#' df_no_units <- as.data.frame(fms, convert_units = FALSE)
#' class(df_no_units$F1_Hz)  # "numeric"
#' }
as.data.frame.AsspDataObj <- function(x, ...,
                                      convert_units = TRUE,
                                      clean_names = TRUE,
                                      na.zeros = FALSE) {

  # Get basic metadata
  sample_rate <- attr(x, "sampleRate")
  start_time <- attr(x, "startTime")
  if (is.null(start_time)) start_time <- 0.0

  # Get track names (template notation)
  track_names <- attr(x, "tracks")
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


#' Convert AsspDataObj to tibble
#'
#' Converts an AsspDataObj to a tibble (tidyverse data frame) with the same
#' features as `as.data.frame.AsspDataObj()`.
#'
#' @param x AsspDataObj. Object to convert.
#' @param ... Additional arguments passed to `as.data.frame.AsspDataObj()`.
#' @param convert_units Logical. If TRUE (default), assign units.
#' @param clean_names Logical. If TRUE (default), use underscore notation.
#' @param na.zeros Logical. If TRUE, convert zeros to NA. Default: FALSE.
#'
#' @return tibble with clean column names, units (optional), and label attributes.
#'
#' @details
#' This is a thin wrapper around `as.data.frame.AsspDataObj()` that converts
#' the result to a tibble. All features (template expansion, clean names,
#' units) are identical.
#'
#' Requires the 'tibble' package.
#'
#' @seealso [as.data.frame.AsspDataObj()]
#'
#' @export
#' @examples
#' \dontrun{
#' library(tibble)
#'
#' fms <- wrassp::forest("audio.wav", numFormants = 4, toFile = FALSE)
#' tbl <- as_tibble(fms)
#'
#' print(tbl)
#' # # A tibble: 289 × 9
#' #    frame_time F1_Hz F2_Hz F3_Hz F4_Hz B1_Hz B2_Hz B3_Hz B4_Hz
#' #         <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#' #  1      0.000  800  1200  2400  3500   80   120   200   250
#' #  2      0.005  820  1220  2380  3480   85   125   210   260
#' #  ...
#' }
as_tibble.AsspDataObj <- function(x, ...,
                                  convert_units = TRUE,
                                  clean_names = TRUE,
                                  na.zeros = FALSE) {

  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package 'tibble' is required for as_tibble(). ",
         "Install it with: install.packages('tibble')",
         call. = FALSE)
  }

  # Convert to data.frame first
  df <- as.data.frame.AsspDataObj(x,
                                  convert_units = convert_units,
                                  clean_names = clean_names,
                                  na.zeros = na.zeros,
                                  ...)

  # Preserve attributes through tibble conversion
  labels <- attr(df, "track_labels")
  descriptions <- attr(df, "track_descriptions")
  sample_rate <- attr(df, "sampleRate")
  start_time <- attr(df, "startTime")

  # Convert to tibble
  tbl <- tibble::as_tibble(df)

  # Restore attributes
  attr(tbl, "track_labels") <- labels
  attr(tbl, "track_descriptions") <- descriptions
  attr(tbl, "sampleRate") <- sample_rate
  attr(tbl, "startTime") <- start_time

  tbl
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
#' - "fo [Hz]", "F1 [Hz]", "H1-H2c [dB]"
#' - Suitable for plot axes
#'
#' **Full labels** (full = TRUE):
#' - "Frequency of oscillation [Hz]"
#' - "First formant frequency [Hz]"
#' - "H1-H2 corrected for formants [dB]"
#' - Suitable for documentation, papers
#'
#' @export
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
