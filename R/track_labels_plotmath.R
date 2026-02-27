# Track Label Expression Generation
#
# Functions to generate plotmath expressions for scientific subscripts in
# ggplot2 axis labels.

#' Convert track name to plotmath expression
#'
#' Converts a cleaned track name to a plotmath expression with proper subscripts
#' for scientific notation.
#'
#' @param col Character. Column name (cleaned format, e.g., "fo_Hz", "F1_Hz").
#' @param full Logical. If TRUE, use full descriptive text. Default: FALSE.
#' @return Expression object suitable for ggplot2 labels.
#'
#' @details
#' Creates plotmath expressions with subscripts following Titze 2015 notation:
#' - `fo_Hz` → expression(f\[o\]~"\[Hz\]")  # f subscript o
#' - `F1_Hz` → expression(F\[1\]~"\[Hz\]")  # F subscript 1
#' - `H1_H2c_dB` → expression(H\[1\]-H\["2c"\]~"\[dB\]")  # H subscript 1 minus H subscript 2c
#' - `LPC12` → expression(LPC\[12\])  # LPC subscript 12
#'
#' The `~` operator adds a small space between the parameter and unit in plotmath.
#'
#' For full descriptive labels, returns a character string instead of expression,
#' as full text like "Frequency of oscillation \[Hz\]" doesn't need subscripts.
#'
#' @keywords internal
#' @examples
#' \dontrun{
#' # Short labels (expressions with subscripts)
#' .col_to_expression("fo_Hz")       # f\[o\] \[Hz\]
#' .col_to_expression("F1_Hz")       # F\[1\] \[Hz\]
#' .col_to_expression("H1_H2c_dB")   # H\[1\]-H\[2c\] \[dB\]
#'
#' # Full labels (text strings)
#' .col_to_expression("fo_Hz", full = TRUE)  # "Frequency of oscillation \[Hz\]"
#' }
.col_to_expression <- function(col, full = FALSE) {

  # For full labels, return descriptive text (no subscripts needed)
  if (full) {
    return(.infer_track_label(col, full = TRUE))
  }

  # Special case: frame_time
  if (col == "frame_time") {
    return(expression("Time [s]"))
  }

  # Pattern: fo (frequency of oscillation)
  if (grepl("^fo(_|$)", col)) {
    unit <- .parse_unit_from_colname(col)
    if (!is.na(unit)) {
      return(bquote(f[o]~"["*.(unit)*"]"))
    } else {
      return(expression(f[o]))
    }
  }

  # Pattern: F<number> (formants)
  if (grepl("^F([0-9]+)(_|$)", col)) {
    n <- as.integer(sub("^F([0-9]+).*", "\\1", col))
    unit <- .parse_unit_from_colname(col)
    if (!is.na(unit)) {
      return(bquote(F[.(n)]~"["*.(unit)*"]"))
    } else {
      return(bquote(F[.(n)]))
    }
  }

  # Pattern: B<number> (bandwidths)
  if (grepl("^B([0-9]+)(_|$)", col)) {
    n <- as.integer(sub("^B([0-9]+).*", "\\1", col))
    unit <- .parse_unit_from_colname(col)
    if (!is.na(unit)) {
      return(bquote(B[.(n)]~"["*.(unit)*"]"))
    } else {
      return(bquote(B[.(n)]))
    }
  }

  # Pattern: H<n1>_H<n2>[c|u] (harmonic differences) - CHECK FIRST!
  if (grepl("^H([0-9]+[a-z]*)_H([0-9]+[a-z]*)(_|$)", col)) {
    h1 <- sub("^H([0-9]+[a-z]*)_H([0-9]+[a-z]*).*", "\\1", col)
    h2 <- sub("^H([0-9]+[a-z]*)_H([0-9]+[a-z]*).*", "\\2", col)
    unit <- .parse_unit_from_colname(col)

    if (!is.na(unit)) {
      return(bquote(H[.(h1)]-H[.(h2)]~"["*.(unit)*"]"))
    } else {
      return(bquote(H[.(h1)]-H[.(h2)]))
    }
  }

  # Pattern: H<n>_A<n> (harmonic to amplitude differences) - CHECK SECOND!
  if (grepl("^H([0-9]+)_A([0-9]+[a-z]*)(_|$)", col)) {
    h <- sub("^H([0-9]+)_A([0-9]+[a-z]*).*", "\\1", col)
    a <- sub("^H([0-9]+)_A([0-9]+[a-z]*).*", "\\2", col)
    unit <- .parse_unit_from_colname(col)

    if (!is.na(unit)) {
      return(bquote(H[.(h)]-A[.(a)]~"["*.(unit)*"]"))
    } else {
      return(bquote(H[.(h)]-A[.(a)]))
    }
  }

  # Pattern: H<number> (harmonics) - CHECK AFTER DIFFERENCES!
  if (grepl("^H([0-9]+[a-z]*)(_|$)", col)) {
    subscript <- sub("^H([0-9]+[a-z]*).*", "\\1", col)
    unit <- .parse_unit_from_colname(col)
    if (!is.na(unit)) {
      return(bquote(H[.(subscript)]~"["*.(unit)*"]"))
    } else {
      return(bquote(H[.(subscript)]))
    }
  }

  # Pattern: A<number> (amplitudes at formants)
  if (grepl("^A([0-9]+)(_|$)", col)) {
    n <- as.integer(sub("^A([0-9]+).*", "\\1", col))
    unit <- .parse_unit_from_colname(col)
    if (!is.na(unit)) {
      return(bquote(A[.(n)]~"["*.(unit)*"]"))
    } else {
      return(bquote(A[.(n)]))
    }
  }

  # Pattern: LPC<number>, ARF<number>, etc. (LP coefficients)
  if (grepl("^(LPC|ARF|LAR|RFC)([0-9]+)$", col)) {
    prefix <- sub("^([A-Z]+)([0-9]+)$", "\\1", col)
    n <- as.integer(sub("^([A-Z]+)([0-9]+)$", "\\2", col))
    return(bquote(.(prefix)[.(n)]))
  }

  # Fallback: use plain text with bracket notation
  .infer_track_label(col, full = FALSE)
}


#' Generate expression labels for data frame columns
#'
#' Creates a named list of plotmath expressions for all columns in a data frame.
#'
#' @param col_names Character vector. Column names.
#' @param full Logical. If TRUE, use full descriptive labels. Default: FALSE.
#' @return Named list mapping column names to expressions.
#'
#' @keywords internal
.generate_track_expressions <- function(col_names, full = FALSE) {
  expressions <- list()

  for (col in col_names) {
    expressions[[col]] <- .col_to_expression(col, full = full)
  }

  expressions
}


#' Get track label as expression for plotting
#'
#' Extended version of `get_track_label()` that returns plotmath expressions
#' with subscripts when `use_subscripts = TRUE`.
#'
#' @param df data.frame or tibble. Data frame from `as.data.frame.AsspDataObj()`.
#' @param col Character. Column name.
#' @param full Logical. If TRUE, return full descriptive label. Default: FALSE.
#' @param use_subscripts Logical. If TRUE, return plotmath expression with
#'   subscripts. If FALSE, return plain text. Default: TRUE.
#'
#' @return Expression object (if use_subscripts = TRUE) or character string.
#'
#' @details
#' This function extends `get_track_label()` with the option to generate
#' plotmath expressions for scientific subscripts in ggplot2.
#'
#' **With subscripts** (use_subscripts = TRUE):
#' - fo → f\[o\] (f subscript o)
#' - F1 → F\[1\] (F subscript 1)
#' - H1-H2c → H\[1\]-H\[2c\] (proper subscripts)
#'
#' **Without subscripts** (use_subscripts = FALSE):
#' - Same as `get_track_label()`
#' - Returns: "fo \[Hz\]", "F1 \[Hz\]", "H1-H2c \[dB\]"
#'
#' For full descriptive labels (full = TRUE), subscripts are not used as the
#' full text doesn't need them (e.g., "Frequency of oscillation \[Hz\]").
#'
#' @examples
#' \dontrun{
#' fms <- wrassp::forest("audio.wav", toFile = FALSE)
#' df <- as.data.frame(fms)
#'
#' # Get expression with subscript
#' expr <- get_track_label_expr(df, "F1_Hz")
#' # Returns: expression(F\[1\]~"\[Hz\]")
#'
#' # Use in ggplot
#' library(ggplot2)
#' ggplot(df, aes(x = frame_time, y = F1_Hz)) +
#'   geom_line() +
#'   labs(y = get_track_label_expr(df, "F1_Hz"))
#' # Y-axis shows: F₁ \[Hz\] (with subscript)
#'
#' # Without subscripts
#' label <- get_track_label_expr(df, "F1_Hz", use_subscripts = FALSE)
#' # Returns: "F1 \[Hz\]"
#' }
get_track_label_expr <- function(df, col, full = FALSE, use_subscripts = TRUE) {

  if (!use_subscripts || full) {
    # Use plain text labels
    return(get_track_label(df, col, full = full))
  }

  # Check if stored in attributes as expression
  expr_attr <- if (full) "track_expression_descriptions" else "track_expressions"
  expressions <- attr(df, expr_attr)

  if (!is.null(expressions) && !is.null(expressions[[col]])) {
    return(expressions[[col]])
  }

  # Generate expression
  .col_to_expression(col, full = full)
}
