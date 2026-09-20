# ggtrack_geoms: ggplot2 geoms for the data objects this package produces
#
# `geom_track()` draws the time-aligned tracks of an `AsspDataObj` (one line per
# track), `geom_spectrogram()` draws a multi-column spectral track as a
# time x frequency raster. Both are ordinary ggplot2 layers, so scales, facets,
# coordinates, legends and themes behave as usual.
#
# ggplot2 hands a layer only the *evaluated aesthetics* of its data (>= 4.0), so
# track selection and the wide -> long rewrite cannot happen in setup_data():
# the `geom_*()` functions prepare the data themselves, and fall back to a data
# function (applied to the plot data during the build) when data is inherited.

#' Check for a suggested ggplot2
#'
#' @param fn Character. Name of the calling function, used in the error message.
#' @return `TRUE`, invisibly. Aborts if ggplot2 is not installed.
#' @keywords internal
.require_ggplot2 <- function(fn) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort(c("Package {.pkg ggplot2} is required for {.fn {fn}}.",
                     "i" = "Install it with: {.code install.packages('ggplot2')}"))
  }
  invisible(TRUE)
}

#' Coerce an AsspDataObj to the table ggplot2 plots from
#'
#' Units columns are not assigned by default: they need the units package
#' attached to scale and the axis labels already carry the unit, so the
#' fortified table stays plain numeric. Pass `convert_units = TRUE` to keep
#' them.
#'
#' @param model AsspDataObj.
#' @param data Ignored; present for method compatibility.
#' @param ... Passed to [as.data.frame.AsspDataObj()].
#' @rdname AsspDataObj
#' @keywords internal
fortify.AsspDataObj <- function(model, data, ...) {
  args <- list(...)
  if (is.null(args$convert_units)) args$convert_units <- FALSE
  do.call(as.data.frame, c(list(model), args))
}

#' Build a geom class from a ggplot2 parent
#'
#' ggplot2 is a suggested dependency, so ggproto classes cannot be created at
#' build time; they are constructed per layer instead.
#'
#' @param name Character. Name of the new ggproto class.
#' @param parent Character. Name of the exported ggplot2 parent class.
#' @return A ggproto object.
#' @keywords internal
.assp_geom <- function(name, parent) {
  ggplot2::ggproto(name, getExportedValue("ggplot2", parent))
}

#' Locate the time column of a track table
#'
#' @param data data.frame. Track table.
#' @param time Character. Preferred name of the time column.
#' @return Character. Name of the column holding time in seconds.
#' @keywords internal
.assp_time_col <- function(data, time = "frame_time") {
  candidates <- intersect(c(time, "frame_time", "time"), names(data))
  if (length(candidates) == 0L) {
    cli::cli_abort(c("No time column in {.arg data}.",
                     "i" = "Expected a column named {.val frame_time}, as written by {.fn as.data.frame.AsspDataObj}.",
                     "i" = "Rename the column, pass it as {.arg time}, or use the {.cls AsspDataObj} itself."))
  }
  candidates[1L]
}

#' Split expanded column names back into coefficient bands
#'
#' `as.data.frame.AsspDataObj()` expands a matrix track into one column per
#' coefficient (`DFT_dB_1`, `DFT_dB_2`, ...). Grouping those names back lets a
#' spectral track keep a single band identity with a coefficient index.
#'
#' @param cols Character vector. Column names of the track table.
#' @return list with `band` (character) and `bin` (integer, `NA` outside a
#'   band).
#' @keywords internal
.assp_band_index <- function(cols) {
  has_index <- grepl("_[0-9]+$", cols)
  bin <- rep(NA_integer_, length(cols))
  bin[has_index] <- as.integer(sub("^.*_([0-9]+)$", "\\1", cols[has_index]))
  band <- ifelse(has_index, sub("_[0-9]+$", "", cols), cols)

  # A prefix only counts as a band when its indices run 1..k without gaps;
  # anything else is an ordinary column that merely ends in a number.
  grouped <- split(bin, band)
  is_band <- vapply(grouped, function(idx) {
    length(idx) > 1L && !anyNA(idx) && identical(sort(idx), seq_along(idx))
  }, logical(1))
  for (stray in names(grouped)[!is_band]) {
    band[band == stray] <- cols[band == stray]
  }

  list(band = band, bin = bin)
}

#' Select track columns by name
#'
#' Accepts the cleaned column name (`"RMS_dB"`), the track template
#' (`"RMS[dB]"`) or the object's track name.
#'
#' @param col,band Character vectors. Candidate column and band names.
#' @param tracks Character vector. Requested names.
#' @param available Character vector. All candidate names, for error messages.
#' @return Logical vector selecting `col`.
#' @keywords internal
.assp_select_tracks <- function(col, band, tracks, available) {
  keep <- col %in% tracks | band %in% tracks
  if (!any(keep)) {
    cli::cli_abort(c("{.arg tracks} does not match any track in {.arg data}.",
                     "i" = "Requested: {.val {tracks}}.",
                     "i" = "Available: {.val {unique(available)}}."))
  }
  keep
}

#' Assemble the long track table
#'
#' @param frame_time Numeric vector. Time of each record, in seconds.
#' @param values List of numeric vectors, one per track column.
#' @param track,band Character vectors. Column name and band name per track.
#' @param bin Integer vector. Coefficient index per track (`NA` if not indexed).
#' @param na.zeros Logical. Convert stored zeros to `NA`.
#' @param origFreq Numeric or NULL. Sample rate the tracks were computed from.
#' @return A data.frame with columns `frame_time`, `value`, `track`, `band` and
#'   `bin`.
#' @keywords internal
.assp_long_frame <- function(frame_time, values, track, band, bin, na.zeros,
                             origFreq = NULL) {
  n <- length(frame_time)
  k <- length(values)
  if (k == 0L) {
    cli::cli_abort(c("No tracks to plot in {.arg data}.",
                     "i" = "Expected numeric track columns, or an {.cls AsspDataObj} holding track data."))
  }
  if (any(lengths(values) != n)) {
    cli::cli_abort("Tracks in {.arg data} do not share a common number of records.")
  }

  value <- unlist(values, use.names = FALSE)
  if (na.zeros) {
    value[value == 0] <- NA_real_
  }

  out <- data.frame(
    frame_time = rep(frame_time, times = k),
    value = value,
    track = rep(track, each = n),
    band = rep(band, each = n),
    bin = rep(bin, each = n),
    stringsAsFactors = FALSE
  )
  attr(out, "origFreq") <- origFreq
  out
}

#' Long track table from an AsspDataObj
#'
#' @param obj AsspDataObj.
#' @param tracks Character or NULL. Tracks to keep.
#' @param na.zeros Logical. Convert stored zeros to `NA`.
#' @return A data.frame as described in [.assp_long_frame()].
#' @keywords internal
.assp_obj_to_long <- function(obj, tracks = NULL, na.zeros = FALSE) {
  sample_rate <- attr(obj, "sampleRate")
  start_time <- attr(obj, "startTime")
  if (is.null(start_time) || !is.finite(start_time)) start_time <- 0
  if (is.null(sample_rate) || !is.finite(sample_rate) || sample_rate <= 0) {
    cli::cli_abort(c("{.cls AsspDataObj} has no usable {.field sampleRate}.",
                     "i" = "Record times cannot be reconstructed without the record rate."))
  }

  # Track templates: the object attribute when set, the list names otherwise
  # (mirrors as.data.frame.AsspDataObj()).
  templates <- attr(obj, "tracks")
  obj_names <- names(obj)
  blocks <- list()

  for (i in seq_along(obj_names)) {
    track_data <- obj[[i]]
    template <- if (!is.null(templates) && length(templates) >= i &&
                    !is.na(templates[i]) && nzchar(templates[i])) {
      templates[i]
    } else {
      obj_names[i]
    }

    if (is.matrix(track_data) && ncol(track_data) > 1L) {
      cols <- .clean_track_names(.expand_track_template(template, ncol(track_data)))
      for (j in seq_len(ncol(track_data))) {
        blocks[[length(blocks) + 1L]] <- list(col = cols[j], band = template,
                                              bin = j, values = track_data[, j])
      }
    } else {
      col <- .clean_track_names(
        if (.has_placeholder(template)) .expand_track_template(template, 1L)[1L] else template
      )
      blocks[[length(blocks) + 1L]] <- list(col = col, band = template,
                                            bin = NA_integer_,
                                            values = as.vector(track_data))
    }
  }

  if (!is.null(tracks)) {
    keep <- .assp_select_tracks(
      col = vapply(blocks, `[[`, character(1), "col"),
      band = vapply(blocks, `[[`, character(1), "band"),
      tracks = tracks,
      available = vapply(blocks, function(b) c(b$col, b$band), character(2))
    )
    blocks <- blocks[keep]
  }

  n_records <- length(blocks[[1L]]$values)
  .assp_long_frame(
    frame_time = start_time + (seq_len(n_records) - 1) / sample_rate,
    values = lapply(blocks, `[[`, "values"),
    track = vapply(blocks, `[[`, character(1), "col"),
    band = vapply(blocks, `[[`, character(1), "band"),
    bin = vapply(blocks, function(b) as.integer(b$bin), integer(1)),
    na.zeros = na.zeros,
    origFreq = attr(obj, "origFreq")
  )
}

#' Long track table from a wide track table
#'
#' @param data data.frame. Wide table, as written by
#'   [as.data.frame.AsspDataObj()].
#' @param tracks Character or NULL. Tracks to keep.
#' @param time Character. Preferred name of the time column.
#' @param na.zeros Logical. Convert stored zeros to `NA`.
#' @return A data.frame as described in [.assp_long_frame()].
#' @keywords internal
.assp_wide_to_long <- function(data, tracks = NULL, time = "frame_time",
                               na.zeros = FALSE) {
  time_col <- .assp_time_col(data, time)
  cols <- names(data)[vapply(data, is.numeric, logical(1))]
  cols <- setdiff(cols, c(time_col, "PANEL", "group", "value", "track", "band", "bin"))
  if (length(cols) == 0L) {
    cli::cli_abort(c("No numeric track columns in {.arg data}.",
                     "i" = "Expected the table written by {.fn as.data.frame.AsspDataObj}."))
  }

  index <- .assp_band_index(cols)
  if (!is.null(tracks)) {
    keep <- .assp_select_tracks(cols, index$band, tracks, c(cols, index$band))
    cols <- cols[keep]
    index$band <- index$band[keep]
    index$bin <- index$bin[keep]
  }

  .assp_long_frame(
    frame_time = data[[time_col]],
    values = unname(data[cols]),
    track = cols,
    band = index$band,
    bin = index$bin,
    na.zeros = na.zeros,
    origFreq = attr(data, "origFreq")
  )
}

#' Long track table from any supported input
#'
#' @param data AsspDataObj, data.frame, or an already long track table.
#' @param tracks Character or NULL. Tracks to keep.
#' @param time Character. Preferred name of the time column.
#' @param na.zeros Logical. Convert stored zeros to `NA`.
#' @return A data.frame with columns `frame_time`, `value`, `track`, `band` and
#'   `bin`, carrying `origFreq` as an attribute when it is known.
#' @keywords internal
.assp_long_data <- function(data, tracks = NULL, time = "frame_time",
                            na.zeros = FALSE) {
  if (inherits(data, "AsspDataObj")) {
    return(.assp_obj_to_long(data, tracks = tracks, na.zeros = na.zeros))
  }
  if (inherits(data, "JsonTrackObj")) {
    cli::cli_abort(c("{.cls JsonTrackObj} holds time-bounded slices, not equally spaced tracks.",
                     "i" = "Plot slices with {.fn ggplot2::geom_segment}, for example",
                     "i" = "{.code aes(x = begin_time, xend = end_time, y = <field>, yend = <field>)}.",
                     "i" = "Tracks that follow the signal come from the {.code trk_*} functions."))
  }
  if (!is.data.frame(data)) {
    cli::cli_abort("{.arg data} must be an {.cls AsspDataObj} or a data.frame, not {.obj_type_friendly {data}}.")
  }

  if (all(c("value", "track") %in% names(data))) {
    # Already long: normalise the abscissa name so the default mapping applies.
    time_col <- .assp_time_col(data, time)
    if (time_col != "frame_time") {
      names(data)[names(data) == time_col] <- "frame_time"
    }
    if (!"band" %in% names(data)) data$band <- data$track
    if (!"bin" %in% names(data)) data$bin <- NA_integer_
    if (na.zeros) data$value[data$value == 0] <- NA_real_
    return(data)
  }

  .assp_wide_to_long(data, tracks = tracks, time = time, na.zeros = na.zeros)
}

#' Warn when a track plot would draw an unreadable number of lines
#'
#' @param long data.frame. Long track table.
#' @return `long`, invisibly.
#' @keywords internal
.assp_warn_many_tracks <- function(long) {
  n_tracks <- length(unique(long$track))
  if (n_tracks > 64L) {
    cli::cli_warn(c("{n_tracks} tracks selected, one line each.",
                    "i" = "A multi-column spectral track is usually drawn with {.fn geom_spectrogram}.",
                    "i" = "Use {.arg tracks} to draw a subset."))
  }
  invisible(long)
}

#' Long time x frequency table from a spectral track
#'
#' @param data AsspDataObj, data.frame, or long track table.
#' @param tracks Character or NULL. Spectral track (band) to use.
#' @param time Character. Preferred name of the time column.
#' @param freq_hz_per_bin Numeric or NULL. Spacing of the coefficients in Hz,
#'   overriding the value derived from `origFreq`.
#' @param na.zeros Logical. Convert stored zeros to `NA`.
#' @return A data.frame with columns `frame_time`, `freq` and `value`.
#' @keywords internal
.assp_spec_data <- function(data, tracks = NULL, time = "frame_time",
                            freq_hz_per_bin = NULL, na.zeros = FALSE) {
  long <- .assp_long_data(data, tracks = tracks, time = time, na.zeros = na.zeros)

  bands <- unique(long$band[!is.na(long$bin)])
  if (length(bands) == 0L) {
    cli::cli_abort(c("No multi-column (spectral) track in {.arg data}.",
                     "i" = "Spectra come from {.fn trk_dft_spectrum}, {.fn trk_lps_spectrum}",
                     "i" = "or {.fn trk_css_spectrum}.",
                     "i" = "Single-value tracks are drawn with {.fn geom_track}."))
  }
  if (length(bands) > 1L) {
    cli::cli_abort(c("{length(bands)} spectral tracks in {.arg data}.",
                     "i" = "Select one with {.arg tracks}, for example {.code tracks = {.val {bands[1]}}}."))
  }

  band <- bands[1L]
  values <- long[long$band == band, , drop = FALSE]
  n_bins <- max(values$bin)

  if (is.null(freq_hz_per_bin)) {
    orig_freq <- attr(long, "origFreq")
    if (is.null(orig_freq) || !is.finite(orig_freq) || orig_freq <= 0) {
      cli::cli_abort(c("Cannot derive the frequency axis of track {.val {band}}.",
                       "i" = "{.arg data} carries no {.val origFreq}, the sample rate the spectrum was computed from.",
                       "i" = "Pass the {.cls AsspDataObj} itself, or set {.arg freq_hz_per_bin}."))
    }
    # SSFF stores a spectrum from 0 Hz to the Nyquist rate; verified against
    # trk_dft_spectrum()/trk_lps_spectrum()/trk_css_spectrum() on a 1 kHz tone.
    freq_hz_per_bin <- orig_freq / (2 * (n_bins - 1))
  }

  data.frame(
    frame_time = values$frame_time,
    freq = (values$bin - 1) * freq_hz_per_bin,
    value = values$value,
    stringsAsFactors = FALSE
  )
}

#' @rdname JsonTrackObj
#' @param model JsonTrackObj.
#' @param data Ignored; present for method compatibility.
#' @keywords internal
fortify.JsonTrackObj <- function(model, data, ...) {
  as.data.frame(model, ...)
}

#' Plot time-aligned tracks of an AsspDataObj
#'
#' A ggplot2 layer that draws the tracks of an `AsspDataObj` (or of the wide
#' table written by [as.data.frame.AsspDataObj()]) as one line per track,
#' against time in seconds. Scales, facets, coordinates, legends and themes
#' apply as with any other layer.
#' @param mapping Set of aesthetic mappings created by [ggplot2::aes()],
#'   applied to the track table the layer prepares (`frame_time`, `value`,
#'   `track`, `band`, `bin`). When `NULL` (default) the tracks are drawn with
#'   `aes(x = frame_time, y = value, colour = track, group = track)`.
#' @param data The data to display: an `AsspDataObj` (from a `trk_*` function,
#'   [read_audio()], [read_ssff()] or [read_track()]), the wide table from
#'   [as.data.frame.AsspDataObj()], or an already long track table with
#'   `frame_time`, `value` and `track` columns. Inherited from the plot when
#'   `NULL`.
#' @param ... Other arguments passed on to [ggplot2::layer()].
#' @param tracks Character. Tracks to plot, matched against the column name
#'   (`"RMS_dB"`), the track template (`"RMS[dB]"`) or the object's track name.
#'   `NULL` (default) plots every track.
#' @param time Character. Preferred name of the time column of `data`. The
#'   plotted table always names it `frame_time`.
#' @param na.zeros Logical. Convert stored zeros to `NA`, so that undefined
#'   frames (unvoiced f0, absent formants) break the line instead of dropping to
#'   zero. Default: `FALSE`.
#' @param na.rm Logical. Remove `NA` values before drawing. Default: `FALSE`.
#' @param show.legend,inherit.aes,stat,position Passed on to
#'   [ggplot2::layer()], with ggplot2's usual defaults except
#'   `inherit.aes = FALSE` (see Details).
#' @return A ggplot2 layer.
#' @details
#' Track data is rewritten into a long table (`frame_time`, `value`, `track`,
#' `band`, `bin`) before the layer is built, because ggplot2 hands a layer only
#' its evaluated aesthetics. A `mapping` therefore refers to that table, and
#' `inherit.aes` defaults to `FALSE`: the mapping goes on the layer, as in
#' `geom_track(data = obj, mapping = aes(x = frame_time, y = value, group = track))`.
#' `tracks` also matches bands, so `tracks = "Fi[Hz]"` draws F1..Fn while
#' `tracks = "F1_Hz"` draws a single formant.
#' For a column-wise plot of the wide table (one column against another), use
#' [ggplot2::geom_line()] or [ggplot2::geom_point()] on
#' [as.data.frame.AsspDataObj()] instead: `ggtrack(df, aes(x = F2_Hz, y = F1_Hz)) + geom_point()`.
#' Layers on different objects compose, so a waveform can be drawn over a
#' spectrogram by adding one [geom_spectrogram()] and one `geom_track()`.
#' The layer carries one row per record and track (a waveform has one row per
#' sample), so window a long recording at read time with
#' `read_audio(wav, begin = …, end = …)`.
#' @seealso [geom_spectrogram()] for multi-column spectral tracks, [ggtrack()]
#'   for automatic track labels, [as.data.frame.AsspDataObj()] for the wide
#'   table.
#' @examples
#' wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(ggplot2)
#'
#'   # every track of the object, one line each
#'   rms <- trk_rms(wav, toFile = FALSE, verbose = FALSE)
#'   ggplot() + geom_track(data = rms)
#'
#'   # a subset, labelled from the track metadata
#'   fms <- trk_formant_forest(wav, numFormants = 3, toFile = FALSE, verbose = FALSE)
#'   ggtrack(fms) + geom_track(tracks = c("F1_Hz", "F2_Hz", "F3_Hz"))
#'
#'   # the raw waveform
#'   audio <- read_audio(wav, samples = TRUE)
#'   ggplot() + geom_track(data = audio)
#'
#'   # a waveform drawn over the spectrogram of the same file
#'   dft <- trk_dft_spectrum(wav, toFile = FALSE, verbose = FALSE)
#'   ggplot() +
#'     geom_spectrogram(data = dft) +
#'     geom_track(data = audio,
#'                mapping = aes(x = frame_time, y = value / 8000 + 3000),
#'                colour = "white", linewidth = 0.2) +
#'     coord_cartesian(ylim = c(0, 6000))
#' }
#' }
#' @export
geom_track <- function(mapping = NULL, data = NULL, ...,
                       tracks = NULL, time = "frame_time", na.zeros = FALSE,
                       na.rm = FALSE, show.legend = NA, inherit.aes = FALSE,
                       stat = "identity", position = "identity") {
  .require_ggplot2("geom_track")

  prep <- function(plot_data) {
    long <- .assp_long_data(plot_data, tracks = tracks, time = time,
                            na.zeros = na.zeros)
    if (is.null(tracks)) .assp_warn_many_tracks(long)
    long
  }
  # A data function is used when the data is inherited: ggplot2 applies it to
  # the plot data during the build, which is where the tracks are visible.
  data <- if (is.null(data)) prep else prep(data)

  if (is.null(mapping)) {
    mapping <- ggplot2::aes(x = frame_time, y = value, colour = track, group = track)
  }

  ggplot2::layer(
    geom = .assp_geom("GeomTrack", "GeomLine"),
    mapping = mapping, data = data, stat = stat, position = position,
    show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

#' Plot a spectral track as a spectrogram
#'
#' A ggplot2 layer that draws a multi-column spectral track of an `AsspDataObj`
#' as a time x frequency raster, with the track values (dB for the package's
#' spectra) on the fill scale.
#' @param mapping Set of aesthetic mappings created by [ggplot2::aes()],
#'   applied to the time x frequency table the layer prepares (`frame_time`,
#'   `freq`, `value`). When `NULL` (default) the spectrum is drawn with
#'   `aes(x = frame_time, y = freq, fill = value)`.
#' @param data The data to display: an `AsspDataObj` holding a spectrum (from
#'   [trk_dft_spectrum()], [trk_lps_spectrum()], [trk_css_spectrum()] or any
#'   other multi-column track), the wide table from
#'   [as.data.frame.AsspDataObj()], or an already long table. Inherited from the
#'   plot when `NULL`.
#' @param ... Other arguments passed on to [ggplot2::layer()].
#' @param tracks Character. The spectral track to draw when `data` holds more
#'   than one multi-column track. `NULL` (default) uses the only one.
#' @param time Character. Preferred name of the time column of `data`. The
#'   plotted table always names it `frame_time`.
#' @param freq_hz_per_bin Numeric. Spacing between the coefficients in Hz,
#'   overriding the spacing derived from the `origFreq` of the object. Needed
#'   when `data` carries no sample rate.
#' @param na.zeros Logical. Convert stored zeros to `NA`, leaving those raster
#'   cells unfilled. Default: `FALSE`.
#' @param na.rm Logical. Remove `NA` values before drawing. Default: `FALSE`.
#' @param interpolate Logical. Interpolate the raster for a smoother image.
#'   Default: `FALSE`.
#' @param show.legend,inherit.aes Passed on to [ggplot2::layer()], with
#'   ggplot2's usual defaults except `inherit.aes = FALSE` (see Details).
#' @return A ggplot2 layer.
#' @details
#' SSFF stores a spectrum from 0 Hz up to the Nyquist rate, so the coefficients
#' of a track with `n` columns sit at `0, bin_hz, ..., (n - 1) * bin_hz` with
#' `bin_hz = origFreq / (2 * (n - 1))`: a 2048-point spectrum of 44.1 kHz audio
#' has 1025 bins of ~21.5 Hz. The track must be a frequency-domain track and
#' `data` must hold exactly one multi-column track, or `tracks` must pick it;
#' multi-column tracks of anything else (formants, LPC coefficients) belong to
#' [geom_track()].
#' The raster is drawn by [ggplot2::geom_raster()]'s engine, which needs evenly
#' spaced times and frequencies (all spectra in this package are evenly spaced).
#' Zoom with [ggplot2::coord_cartesian()] rather than scale limits to keep the
#' raster intact.
#' @seealso [geom_track()] for single-value tracks and waveforms, [ggtrack()]
#'   for automatic track labels.
#' @examples
#' wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(ggplot2)
#'
#'   dft <- trk_dft_spectrum(wav, toFile = FALSE, verbose = FALSE)
#'   ggplot() +
#'     geom_spectrogram(data = dft) +
#'     coord_cartesian(ylim = c(0, 5000))
#'
#'   # an object inherited from the plot works the same way
#'   lps <- trk_lps_spectrum(wav, toFile = FALSE, verbose = FALSE)
#'   ggtrack(lps) +
#'     geom_spectrogram() +
#'     labs(y = "Frequency [Hz]") +
#'     scale_fill_viridis_c()
#' }
#' }
#' @export
geom_spectrogram <- function(mapping = NULL, data = NULL, ...,
                             tracks = NULL, time = "frame_time",
                             freq_hz_per_bin = NULL, na.zeros = FALSE,
                             na.rm = FALSE, interpolate = FALSE,
                             show.legend = NA, inherit.aes = FALSE) {
  .require_ggplot2("geom_spectrogram")

  prep <- function(plot_data) {
    .assp_spec_data(plot_data, tracks = tracks, time = time,
                    freq_hz_per_bin = freq_hz_per_bin, na.zeros = na.zeros)
  }
  data <- if (is.null(data)) prep else prep(data)

  if (is.null(mapping)) {
    mapping <- ggplot2::aes(x = frame_time, y = freq, fill = value)
  }

  ggplot2::layer(
    geom = .assp_geom("GeomSpectrogram", "GeomRaster"),
    mapping = mapping, data = data, stat = "identity", position = "identity",
    show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, interpolate = interpolate, ...)
  )
}
