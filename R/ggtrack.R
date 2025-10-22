# ggtrack: Auto-labeled ggplot2 for AsspDataObj Track Data
#
# Convenience function for creating ggplot2 plots with automatic axis labels
# from AsspDataObj track metadata.

#' Create ggplot with automatic track labels
#'
#' Creates a ggplot2 plot with automatic axis labels derived from track names.
#' This is a convenience wrapper around `ggplot()` that automatically extracts
#' and applies appropriate labels for acoustic track data.
#'
#' @param data data.frame or tibble. Data from `as.data.frame.AsspDataObj()`.
#' @param mapping ggplot2::aes() specification. Aesthetic mappings.
#' @param ... Additional ggplot2 layers to add (geoms, scales, themes, etc.).
#' @param full_labels Logical. If TRUE, use full descriptive labels. If FALSE
#'   (default), use short labels suitable for plot axes.
#'
#' @return A ggplot object with automatic axis labels.
#'
#' @details
#' This function simplifies plotting of acoustic track data by automatically
#' generating appropriate axis labels based on column names. It extracts the
#' x and y variables from the aesthetic mapping and applies labels using
#' `get_track_label()`.
#'
#' **Short labels** (full_labels = FALSE, default):
#' - "fo [Hz]", "F1 [Hz]", "H1-H2c [dB]"
#' - Concise, suitable for most plots
#'
#' **Full labels** (full_labels = TRUE):
#' - "Frequency of oscillation [Hz]"
#' - "First formant frequency [Hz]"
#' - "H1-H2 corrected for formants [dB]"
#' - Descriptive, suitable for publications
#'
#' **Additional layers** can be added using the `...` argument or by adding
#' to the returned ggplot object with `+`.
#'
#' @seealso
#' - [as.data.frame.AsspDataObj()] for creating the data frame
#' - [get_track_label()] for label extraction
#' - [ggplot2::ggplot()] for the underlying plotting function
#'
#' @export
#' @examples
#' \dontrun{
#' library(ggplot2)
#'
#' # Get pitch track
#' pitch <- wrassp::ksvF0("audio.wav", toFile = FALSE)
#' df <- as.data.frame(pitch)
#'
#' # Quick plot with automatic labels
#' ggtrack(df, aes(x = frame_time, y = fo_Hz)) +
#'   geom_line()
#' # Y-axis automatically labeled: "fo [Hz]"
#'
#' # With full descriptive labels
#' ggtrack(df, aes(x = frame_time, y = fo_Hz), full_labels = TRUE) +
#'   geom_line() +
#'   theme_minimal()
#' # Y-axis: "Frequency of oscillation [Hz]"
#'
#' # Adding layers via ...
#' fms <- wrassp::forest("vowel.wav", toFile = FALSE)
#' df_fms <- as.data.frame(fms)
#'
#' ggtrack(df_fms, aes(x = F2_Hz, y = F1_Hz),
#'         geom_point(alpha = 0.5),
#'         scale_x_reverse(),
#'         scale_y_reverse(),
#'         theme_minimal())
#' # Automatic labels: "F1 [Hz]", "F2 [Hz]"
#'
#' # Can also add layers with +
#' ggtrack(df, aes(x = frame_time, y = fo_Hz)) +
#'   geom_line(color = "steelblue") +
#'   geom_smooth(method = "loess", se = FALSE, color = "red") +
#'   theme_bw()
#' }
ggtrack <- function(data, mapping = ggplot2::aes(), ..., full_labels = FALSE) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for ggtrack(). ",
         "Install it with: install.packages('ggplot2')",
         call. = FALSE)
  }

  if (!requireNamespace("rlang", quietly = TRUE)) {
    stop("Package 'rlang' is required for ggtrack(). ",
         "Install it with: install.packages('rlang')",
         call. = FALSE)
  }

  # Extract x and y variables from mapping
  x_var <- if (!is.null(mapping$x)) {
    rlang::as_name(mapping$x)
  } else {
    NULL
  }

  y_var <- if (!is.null(mapping$y)) {
    rlang::as_name(mapping$y)
  } else {
    NULL
  }

  # Get labels
  x_label <- if (!is.null(x_var)) {
    get_track_label(data, x_var, full = full_labels)
  } else {
    ggplot2::waiver()
  }

  y_label <- if (!is.null(y_var)) {
    get_track_label(data, y_var, full = full_labels)
  } else {
    ggplot2::waiver()
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
