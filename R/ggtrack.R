# ggtrack: Auto-labeled ggplot2 for AsspDataObj Track Data
#
# Convenience function for creating ggplot2 plots with automatic axis labels
# from AsspDataObj track metadata.

#' Create ggplot with automatic track labels
#'
#' Creates a ggplot2 plot with automatic axis labels derived from track names.
#' This is a convenience wrapper around `ggplot()` that automatically extracts
#' and applies appropriate labels for acoustic track data.
#' @param data data.frame, tibble, `AsspDataObj` or `JsonTrackObj`. Track data,
#'   either as a table from [as.data.frame.AsspDataObj()] or as the object
#'   itself.
#' @param mapping ggplot2::aes() specification. Aesthetic mappings.
#' @param ... Additional ggplot2 layers to add (geoms, scales, themes, etc.).
#' @param full_labels Logical. If TRUE, use full descriptive labels. If FALSE
#'   (default), use short labels suitable for plot axes.
#' @param use_subscripts Logical. If TRUE (default), use plotmath expressions
#'   with subscripts (fo and F1 rendered with subscript digits). If FALSE,
#'   use plain text.
#' @return A ggplot object with automatic axis labels.
#' @details
#' This function simplifies plotting of acoustic track data by automatically
#' generating appropriate axis labels based on column names. It extracts the
#' x and y variables from the aesthetic mapping and applies labels using
#' `get_track_label()`.
#' When no mapping is given, the axes are labelled for the layers this package
#' provides: `x` becomes "Time [s]" for a table with a `frame_time` column, and
#' `y` becomes the track label when the data holds exactly one track
#' (`ggtrack(rms) + geom_track()`).
#' **Short labels** (full_labels = FALSE, default):
#' - "fo \[Hz\]", "F1 \[Hz\]", "H1-H2c \[dB\]"
#' - Concise, suitable for most plots
#' **Full labels** (full_labels = TRUE):
#' - "Frequency of oscillation \[Hz\]"
#' - "First formant frequency \[Hz\]"
#' - "H1-H2 corrected for formants \[dB\]"
#' - Descriptive, suitable for publications
#' **Additional layers** can be added using the `...` argument or by adding
#' to the returned ggplot object with `+`.
#' @seealso
#' - [geom_track()] and [geom_spectrogram()] for the layers to combine with
#' - [as.data.frame.AsspDataObj()] for the wide track table
#' - [get_track_label()] for label extraction
#' - [ggplot2::ggplot()] for the underlying plotting function
#' @examples
#' wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(ggplot2)
#'
#'   # tracks of the object, labels taken from the track metadata
#'   f0 <- trk_pitch_ksv(wav, toFile = FALSE, verbose = FALSE)
#'   ggtrack(f0) + geom_track()
#'
#'   # y-axis: "fo [Hz]"
#'   # a spectrogram of a spectrum track
#'   dft <- trk_dft_spectrum(wav, toFile = FALSE, verbose = FALSE)
#'   ggtrack(dft) + geom_spectrogram() + labs(y = "Frequency [Hz]")
#'
#'   # the wide table with a column-wise mapping. convert_units = FALSE keeps
#'   # the columns plain numeric; unit-assigned columns are "units" objects,
#'   # which ggplot2 can only scale with the units package attached.
#'   fms <- trk_formant_forest(wav, numFormants = 3, toFile = FALSE, verbose = FALSE)
#'   df_fms <- as.data.frame(fms, convert_units = FALSE)
#'   ggtrack(df_fms, aes(x = F2_Hz, y = F1_Hz)) +
#'     geom_point(alpha = 0.5) +
#'     scale_x_reverse() +
#'     scale_y_reverse() +
#'     theme_minimal()
#'   # Automatic labels: "F1 [Hz]", "F2 [Hz]"
#' }
#' }
#' @export
ggtrack <- function(data, mapping = ggplot2::aes(), ...,
                   full_labels = FALSE,
                   use_subscripts = TRUE) {

  .require_ggplot2("ggtrack")

  # Labels are looked up on the wide table, which carries the label attributes
  # written by as.data.frame.AsspDataObj(); the plot itself takes `data` as
  # given (ggplot2 coerces it with fortify()).
  label_data <- if (inherits(data, "AsspDataObj") || inherits(data, "JsonTrackObj")) {
    as.data.frame(data, convert_units = FALSE)
  } else {
    data
  }

  # Extract x and y variables from mapping (NULL when the aesthetic is not a
  # plain column name, e.g. y = log(F1_Hz): ggplot2's own label applies then)
  x_var <- .assp_aes_name(mapping$x)
  y_var <- .assp_aes_name(mapping$y)

  # Get labels (with or without subscripts)
  if (use_subscripts && !full_labels) {
    # Use plotmath expressions with subscripts
    x_label <- if (!is.null(x_var)) {
      get_track_label_expr(label_data, x_var, full = FALSE, use_subscripts = TRUE)
    } else {
      .default_time_label(label_data)
    }

    y_label <- if (!is.null(y_var)) {
      get_track_label_expr(label_data, y_var, full = FALSE, use_subscripts = TRUE)
    } else {
      .default_track_label(label_data, use_subscripts = TRUE)
    }
  } else {
    # Use plain text labels
    x_label <- if (!is.null(x_var)) {
      get_track_label(label_data, x_var, full = full_labels)
    } else {
      .default_time_label(label_data)
    }

    y_label <- if (!is.null(y_var)) {
      get_track_label(label_data, y_var, full = full_labels)
    } else {
      .default_track_label(label_data, full = full_labels,
                           use_subscripts = use_subscripts)
    }
  }

  # Create base plot
  p <- ggplot2::ggplot(data, mapping)

  # Add labels
  p <- p + ggplot2::labs(x = x_label, y = y_label)

  # Add additional layers from ...
  dots <- list(...)
  if (length(dots) > 0) {
    for (layer in dots) {
      p <- p + layer
    }
  }

  p
}

#' Column name behind an aesthetic mapping
#'
#' @param aes Aesthetic from a [ggplot2::aes()] mapping (a quosure), or `NULL`.
#' @return Character. The column name when the aesthetic is a single column,
#'   `NULL` for anything else (constants, calls like `log(F1_Hz)`), where
#'   ggplot2's own label is the right one.
#' @keywords internal
.assp_aes_name <- function(aes) {
  expr <- if (rlang::is_quosure(aes)) rlang::quo_get_expr(aes) else aes
  if (rlang::is_symbol(expr)) rlang::as_string(expr) else NULL
}

#' Default x-axis label for a track table
#'
#' @param data data.frame or NULL. Wide track table.
#' @return "Time \[s\]" for a table with a `frame_time` column, a waiver
#'   otherwise.
#' @keywords internal
.default_time_label <- function(data) {
  if (is.data.frame(data) && "frame_time" %in% names(data)) {
    return("Time [s]")
  }
  ggplot2::waiver()
}

#' Default y-axis label for a track table
#'
#' Only a table holding exactly one track is unambiguous; anything else keeps
#' ggplot2's own default label. The sampled waveform track ("audio") is
#' amplitude, not a measured quantity with a label of its own.
#'
#' @param data data.frame or NULL. Wide track table.
#' @param full Logical. Use the full descriptive label.
#' @param use_subscripts Logical. Use plotmath subscripts in the short label.
#' @return The label of the single track, a waiver otherwise.
#' @keywords internal
.default_track_label <- function(data, full = FALSE, use_subscripts = TRUE) {
  if (is.data.frame(data)) {
    tracks <- setdiff(names(data), c("frame_time", "begin_time", "end_time"))
    if (length(tracks) == 1L) {
      track <- tracks[1L]
      if (identical(track, "audio")) return("Amplitude")
      return(if (use_subscripts && !full) {
        get_track_label_expr(data, track, use_subscripts = TRUE)
      } else {
        get_track_label(data, track, full = full)
      })
    }
  }
  ggplot2::waiver()
}
