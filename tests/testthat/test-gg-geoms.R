# ggplot2 geoms for AsspDataObj track data: geom_track(), geom_spectrogram(),
# the fortify() methods behind them, and the labels ggtrack() derives.

# The DSP objects are built once per file: each trk_* call runs the analysis.
sample_objects <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
      cache <<- if (!nzchar(wav)) {
        list(wav = "")
      } else {
        list(
          wav = wav,
          rms = trk_rms(wav, toFile = FALSE, verbose = FALSE),
          fms = trk_formant_forest(wav, numFormants = 4, toFile = FALSE, verbose = FALSE),
          dft = trk_dft_spectrum(wav, toFile = FALSE, verbose = FALSE, windowShift = 20),
          f0 = trk_pitch_ksv(wav, toFile = FALSE, verbose = FALSE)
        )
      }
    }
    cache
  }
})

skip_if_no_gg <- function() {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("superassp")
  if (!nzchar(sample_objects()$wav)) testthat::skip("Test file not found")
}

# The layer data is the table the geom prepared, before scales are applied.
prepared <- function(layer) layer$data

test_that("geom_track draws every track against time in seconds", {
  skip_if_no_gg()
  obj <- sample_objects()$rms
  df <- as.data.frame(obj, convert_units = FALSE)

  built <- ggplot2::ggplot_build(ggplot2::ggplot() + geom_track(data = obj))
  panel <- built@data[[1]]

  expect_equal(nrow(panel), nrow(df))
  expect_equal(panel$x, df$frame_time)
  expect_equal(panel$y, df$RMS_dB)
  expect_length(unique(panel$group), 1)
  expect_length(unique(panel$colour), 1)
})

test_that("geom_track works on data inherited from the plot", {
  skip_if_no_gg()
  obj <- sample_objects()$fms

  supplied <- ggplot2::ggplot_build(ggplot2::ggplot() + geom_track(data = obj))@data[[1]]
  inherited <- ggplot2::ggplot_build(ggplot2::ggplot(obj) + geom_track())@data[[1]]

  expect_equal(inherited$x, supplied$x)
  expect_equal(inherited$y, supplied$y)
  expect_length(unique(inherited$group), 8)  # F1-F4 and B1-B4
})

test_that("tracks are selected by column name, band template or object name", {
  skip_if_no_gg()
  obj <- sample_objects()$fms
  fields <- function(layer) unique(prepared(layer)$track)

  expect_equal(fields(geom_track(data = obj, tracks = "F1_Hz")), "F1_Hz")
  expect_setequal(fields(geom_track(data = obj, tracks = "Fi[Hz]")),
                  c("F1_Hz", "F2_Hz", "F3_Hz", "F4_Hz"))
  expect_setequal(fields(geom_track(data = obj, tracks = c("F2_Hz", "Bi[Hz]"))),
                  c("F2_Hz", "B1_Hz", "B2_Hz", "B3_Hz", "B4_Hz"))

  expect_error(geom_track(data = obj, tracks = "nope_Hz"),
               "does not match any track")
})

test_that("geom_track plots a wide or long track table", {
  skip_if_no_gg()
  obj <- sample_objects()$fms
  wide <- as.data.frame(obj, convert_units = FALSE)

  from_wide <- ggplot2::ggplot_build(
    ggplot2::ggplot(wide) + geom_track(tracks = c("F1_Hz", "F2_Hz"))
  )@data[[1]]
  from_object <- ggplot2::ggplot_build(
    ggplot2::ggplot() + geom_track(data = obj, tracks = c("F1_Hz", "F2_Hz"))
  )@data[[1]]

  expect_equal(from_wide$y, from_object$y)

  reused_long <- prepared(geom_track(data = obj, tracks = "F1_Hz"))
  reused <- ggplot2::ggplot_build(ggplot2::ggplot(reused_long) + geom_track())@data[[1]]
  expect_equal(reused$y, reused_long$value)
  expect_equal(reused$x, reused_long$frame_time)
})

test_that("a layer mapping refers to the prepared track table", {
  skip_if_no_gg()
  obj <- sample_objects()$fms

  # the documented long columns: frame_time, value, track (and band, bin)
  long <- prepared(geom_track(data = obj))
  expect_setequal(names(long), c("frame_time", "value", "track", "band", "bin"))

  spelled_out <- ggplot2::ggplot_build(
    ggplot2::ggplot() + geom_track(
      data = obj,
      mapping = ggplot2::aes(x = frame_time, y = value, colour = track, group = track)
    )
  )@data[[1]]
  defaulted <- ggplot2::ggplot_build(ggplot2::ggplot() + geom_track(data = obj))@data[[1]]

  expect_equal(spelled_out$y, defaulted$y)
  expect_equal(spelled_out$colour, defaulted$colour)

  # a column of the wide table is not part of that table
  expect_error(
    ggplot2::ggplot_build(
      ggplot2::ggplot() + geom_track(data = obj, mapping = ggplot2::aes(y = F1_Hz))
    )
  )
})

test_that("each track gets a legend key", {
  skip_if_no_gg()
  obj <- sample_objects()$fms

  plot <- ggplot2::ggplot() + geom_track(data = obj, tracks = "Fi[Hz]")
  guide <- ggplot2::get_guide_data(plot, "colour")

  expect_setequal(guide$.label, c("F1_Hz", "F2_Hz", "F3_Hz", "F4_Hz"))
})

test_that("stored zeros turn into line breaks on request", {
  skip_if_no_gg()
  obj <- sample_objects()$fms
  values <- prepared(geom_track(data = obj, tracks = "F1_Hz"))$value

  with_zeros <- prepared(geom_track(data = obj, tracks = "F1_Hz", na.zeros = TRUE))$value

  expect_gt(sum(values == 0), 0)
  expect_equal(sum(is.na(with_zeros)), sum(values == 0))
  expect_equal(with_zeros[values != 0], values[values != 0])
})

test_that("geom_spectrogram maps coefficients onto the Nyquist grid", {
  skip_if_no_gg()
  obj <- sample_objects()$dft
  values <- obj[[1]]
  n_bins <- ncol(values)
  bin_hz <- attr(obj, "origFreq") / (2 * (n_bins - 1))

  spec <- prepared(geom_spectrogram(data = obj))

  expect_equal(nrow(spec), nrow(values) * n_bins)
  expect_equal(range(spec$freq), c(0, attr(obj, "origFreq") / 2))
  expect_equal(sort(unique(spec$freq))[2] - sort(unique(spec$freq))[1], bin_hz,
               tolerance = 1e-9)

  # every row carries the matrix value of its own frame and coefficient
  frame <- round((spec$frame_time - attr(obj, "startTime")) * attr(obj, "sampleRate")) + 1
  bin <- round(spec$freq / bin_hz) + 1
  expect_equal(max(frame), nrow(values))
  expect_equal(max(bin), n_bins)
  expect_equal(spec$value, values[cbind(frame, bin)])

  # an explicit coefficient spacing overrides the one derived from the object
  coarse <- prepared(geom_spectrogram(data = obj, freq_hz_per_bin = 100))
  expect_equal(coarse$freq, (bin - 1) * 100)
})

test_that("geom_spectrogram works on the wide table and the raster it draws", {
  skip_if_no_gg()
  obj <- sample_objects()$dft
  wide <- as.data.frame(obj, convert_units = FALSE)

  spec <- prepared(geom_spectrogram(data = obj))
  from_wide <- prepared(geom_spectrogram(data = wide))
  expect_equal(from_wide$freq, spec$freq)
  expect_equal(from_wide$value, spec$value)

  panel <- ggplot2::ggplot_build(
    ggplot2::ggplot() + geom_spectrogram(data = obj) +
      ggplot2::coord_cartesian(ylim = c(0, 4000))
  )@data[[1]]
  expect_length(unique(panel$PANEL), 1)
  expect_equal(nrow(panel), nrow(spec))
  expect_equal(sort(unique(panel$y)), sort(unique(spec$freq)))
})

test_that("geom_spectrogram rejects tracks that carry no frequency axis", {
  skip_if_no_gg()

  expect_error(geom_spectrogram(data = sample_objects()$rms),
               "No multi-column \\(spectral\\) track")
  expect_error(geom_spectrogram(data = sample_objects()$fms),
               "2 spectral tracks")
  expect_error(geom_spectrogram(data = sample_objects()$fms, tracks = "nope"),
               "does not match any track")

  plain <- data.frame(frame_time = seq(0, 0.3, by = 0.1), SPC_dB_1 = 1:4, SPC_dB_2 = 5:8)
  expect_error(geom_spectrogram(data = plain), "freq_hz_per_bin")
})

test_that("fortify makes the data objects usable as plot data", {
  skip_if_no_gg()
  obj <- sample_objects()$fms

  fortified <- ggplot2::fortify(obj)
  expect_s3_class(fortified, "data.frame")
  expect_true(all(c("frame_time", "F1_Hz") %in% names(fortified)))
  expect_true(is.numeric(fortified$F1_Hz))

  # column-wise plotting keeps working, with no rewrite by the layer
  panel <- ggplot2::ggplot_build(
    ggplot2::ggplot(obj, ggplot2::aes(x = frame_time, y = F1_Hz)) + ggplot2::geom_point()
  )@data[[1]]
  expect_equal(nrow(panel), nrow(fortified))

  # JSTF slices are a table of bounded intervals, not tracks
  jstf <- lst_voice_report(sample_objects()$wav, toFile = FALSE, verbose = FALSE,
                           return_jstf = TRUE)
  expect_true("begin_time" %in% names(ggplot2::fortify(jstf)))
  expect_error(geom_track(data = jstf), "time-bounded slices")
})

test_that("ggtrack labels the axes from the track metadata", {
  skip_if_no_gg()
  obj <- sample_objects()$f0
  wide <- as.data.frame(obj, convert_units = FALSE)

  # an explicit mapping labels both axes from the track names
  mapped <- ggplot2::ggplot_build(
    ggtrack(wide, ggplot2::aes(x = frame_time, y = fo_Hz)) + ggplot2::geom_line()
  )@plot@labels
  expect_equal(mapped$x, get_track_label_expr(wide, "frame_time"))
  expect_equal(mapped$y, get_track_label_expr(wide, "fo_Hz"))

  # without a mapping the single track of the object supplies the y label
  single <- ggplot2::ggplot_build(ggtrack(obj) + geom_track())@plot@labels
  expect_equal(single$x, "Time [s]")
  expect_equal(single$y, get_track_label_expr(wide, "fo_Hz"))
})

test_that("ggtrack labels both axes of a plain column-wise plot", {
  skip_if_no_gg()
  wide <- as.data.frame(sample_objects()$fms, convert_units = FALSE)
  plotted <- function(...) {
    ggplot2::ggplot_build(
      ggtrack(wide, ggplot2::aes(x = F2_Hz, y = F1_Hz), ...) + ggplot2::geom_point()
    )@plot@labels
  }

  # subscripts by default
  labels <- plotted()
  expect_equal(labels$x, get_track_label_expr(wide, "F2_Hz"))
  expect_equal(labels$y, get_track_label_expr(wide, "F1_Hz"))

  # plain text on request
  plain <- plotted(use_subscripts = FALSE)
  expect_equal(plain$x, get_track_label(wide, "F2_Hz"))
  expect_equal(plain$y, get_track_label(wide, "F1_Hz"))
})

test_that("ggtrack leaves transformed aesthetics to ggplot2", {
  skip_if_no_gg()
  wide <- as.data.frame(sample_objects()$fms, convert_units = FALSE)

  labels <- ggplot2::ggplot_build(
    ggtrack(wide, ggplot2::aes(x = frame_time, y = log(F1_Hz))) + ggplot2::geom_point()
  )@plot@labels

  expect_equal(labels$x, get_track_label_expr(wide, "frame_time"))
  expect_equal(labels$y, "log(F1_Hz)")
})
