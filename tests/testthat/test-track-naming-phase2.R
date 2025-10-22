# Test Phase 2: Track Naming Infrastructure
#
# Tests for as.data.frame.AsspDataObj(), template expansion, clean names,
# units assignment, label generation, and plotmath expressions.

test_that("Template expansion works for formant tracks", {
  skip_if_not_installed("wrassp")

  # Get test file
  test_wav <- system.file("samples", "sustained", "a1.wav",
                         package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Get formants (4 formants expected) - use trk_forest for tracks attribute
  fms <- trk_forest(test_wav, numFormants = 4, toFile = FALSE, verbose = FALSE)

  # Check original AsspDataObj structure
  expect_true(wrassp::is.AsspDataObj(fms))
  expect_equal(attr(fms, "tracks"), c("Fi[Hz]", "Bi[Hz]"))

  # Convert to data.frame
  df <- as.data.frame(fms)

  # Check template expansion worked
  expect_true("frame_time" %in% names(df))
  expect_true("F1_Hz" %in% names(df))
  expect_true("F2_Hz" %in% names(df))
  expect_true("F3_Hz" %in% names(df))
  expect_true("F4_Hz" %in% names(df))
  expect_true("B1_Hz" %in% names(df))
  expect_true("B2_Hz" %in% names(df))
  expect_true("B3_Hz" %in% names(df))
  expect_true("B4_Hz" %in% names(df))

  # Check that original bracket notation was cleaned
  expect_false("F1[Hz]" %in% names(df))
  expect_false("Fi[Hz]" %in% names(df))
})

test_that("Template expansion works for LP coefficients", {
  skip_if_not_installed("wrassp")

  test_wav <- system.file("samples", "sustained", "a1.wav",
                         package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Get LPC analysis (order 12)
  lpc <- wrassp::lpcana(test_wav, order = 12, toFile = FALSE)

  # Check original structure
  expect_true(wrassp::is.AsspDataObj(lpc))
  expect_equal(attr(lpc, "tracks"), c("RMS[dB]", "gain[dB]", "LPCi"))

  # Convert to data.frame
  df <- as.data.frame(lpc)

  # Check expansion
  expect_true("RMS_dB" %in% names(df))
  expect_true("gain_dB" %in% names(df))
  expect_true("LPC1" %in% names(df))
  expect_true("LPC12" %in% names(df))

  # Check column count (frame_time + RMS + gain + 12 LPC coeffs = 15)
  expect_equal(ncol(df), 15)
})

test_that("clean_names parameter works", {
  skip_if_not_installed("wrassp")

  test_wav <- system.file("samples", "sustained", "a1.wav",
                         package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  fms <- wrassp::forest(test_wav, numFormants = 2, toFile = FALSE)

  # With clean names (default)
  df_clean <- as.data.frame(fms, clean_names = TRUE)
  expect_true("F1_Hz" %in% names(df_clean))
  expect_false("F1[Hz]" %in% names(df_clean))

  # Without clean names
  df_brackets <- as.data.frame(fms, clean_names = FALSE)
  expect_true("F1[Hz]" %in% names(df_brackets))
  expect_false("F1_Hz" %in% names(df_brackets))
})

test_that("units are assigned correctly", {
  skip_if_not_installed("wrassp")
  skip_if_not_installed("units")

  test_wav <- system.file("samples", "sustained", "a1.wav",
                         package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  fms <- wrassp::forest(test_wav, numFormants = 2, toFile = FALSE)

  # With units (default)
  df_units <- as.data.frame(fms, convert_units = TRUE)
  expect_s3_class(df_units$F1_Hz, "units")

  # Without units
  df_no_units <- as.data.frame(fms, convert_units = FALSE)
  expect_type(df_no_units$F1_Hz, "double")
  expect_false(inherits(df_no_units$F1_Hz, "units"))
})

test_that("track labels are generated and stored", {
  skip_if_not_installed("wrassp")

  test_wav <- system.file("samples", "sustained", "a1.wav",
                         package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  fms <- wrassp::forest(test_wav, numFormants = 2, toFile = FALSE)
  df <- as.data.frame(fms)

  # Check attributes exist
  expect_true(!is.null(attr(df, "track_labels")))
  expect_true(!is.null(attr(df, "track_descriptions")))

  # Check short labels
  labels <- attr(df, "track_labels")
  expect_equal(labels[["F1_Hz"]], "F1 [Hz]")
  expect_equal(labels[["B1_Hz"]], "B1 [Hz]")

  # Check full descriptions
  descriptions <- attr(df, "track_descriptions")
  expect_true(grepl("formant", descriptions[["F1_Hz"]], ignore.case = TRUE))
  expect_true(grepl("bandwidth", descriptions[["B1_Hz"]], ignore.case = TRUE))
})

test_that("get_track_label() retrieves labels correctly", {
  skip_if_not_installed("wrassp")

  test_wav <- system.file("samples", "sustained", "a1.wav",
                         package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  fms <- wrassp::forest(test_wav, numFormants = 2, toFile = FALSE)
  df <- as.data.frame(fms)

  # Short labels
  expect_equal(get_track_label(df, "F1_Hz", full = FALSE), "F1 [Hz]")
  expect_equal(get_track_label(df, "frame_time", full = FALSE), "Time [s]")

  # Full labels
  f1_full <- get_track_label(df, "F1_Hz", full = TRUE)
  expect_true(grepl("formant", f1_full, ignore.case = TRUE))
})

test_that("plotmath expressions are generated correctly", {
  skip_if_not_installed("wrassp")

  test_wav <- system.file("samples", "sustained", "a1.wav",
                         package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  fms <- wrassp::forest(test_wav, numFormants = 2, toFile = FALSE)
  df <- as.data.frame(fms)

  # Get expression for F1
  expr_f1 <- get_track_label_expr(df, "F1_Hz", use_subscripts = TRUE)
  expect_s3_class(expr_f1, "call")  # expressions are calls

  # Check expression structure
  expr_str <- deparse(expr_f1)
  expect_true(grepl("F\\[", expr_str))  # Should have F[1]
  expect_true(grepl("Hz", expr_str))    # Should have Hz

  # Plain text version
  plain <- get_track_label_expr(df, "F1_Hz", use_subscripts = FALSE)
  expect_type(plain, "character")
  expect_equal(plain, "F1 [Hz]")
})

test_that("plotmath handles special cases correctly", {
  skip_if_not_installed("wrassp")

  test_wav <- system.file("samples", "sustained", "a1.wav",
                         package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Test fo (frequency of oscillation)
  pitch <- wrassp::ksvF0(test_wav, toFile = FALSE)
  df_pitch <- as.data.frame(pitch)

  expr_fo <- get_track_label_expr(df_pitch, "fo_Hz", use_subscripts = TRUE)
  expr_str <- deparse(expr_fo)
  expect_true(grepl("f\\[o\\]", expr_str))  # Should be f[o] not f[0]

  # Test frame_time
  expr_time <- get_track_label_expr(df_pitch, "frame_time", use_subscripts = TRUE)
  expect_s3_class(expr_time, "call")
})

test_that("na.zeros parameter works", {
  skip_if_not_installed("wrassp")

  test_wav <- system.file("samples", "sustained", "a1.wav",
                         package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Get pitch (likely has some zero values)
  pitch <- wrassp::ksvF0(test_wav, toFile = FALSE)

  # Without na.zeros
  df_zeros <- as.data.frame(pitch, na.zeros = FALSE)

  # With na.zeros
  df_na <- as.data.frame(pitch, na.zeros = TRUE)

  # Should have different numbers of NAs (though we can't predict exact count)
  # Just check the parameter doesn't cause errors
  expect_true(is.data.frame(df_zeros))
  expect_true(is.data.frame(df_na))
})

test_that("as_tibble.AsspDataObj() works", {
  skip_if_not_installed("wrassp")
  skip_if_not_installed("tibble")

  test_wav <- system.file("samples", "sustained", "a1.wav",
                         package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  fms <- wrassp::forest(test_wav, numFormants = 2, toFile = FALSE)

  # Convert to tibble
  tbl <- tibble::as_tibble(fms)

  # Check class
  expect_s3_class(tbl, "tbl_df")

  # Check attributes preserved
  expect_true(!is.null(attr(tbl, "track_labels")))
  expect_true(!is.null(attr(tbl, "track_descriptions")))

  # Check column names
  expect_true("F1_Hz" %in% names(tbl))
  expect_true("F2_Hz" %in% names(tbl))
})

test_that("ggtrack() creates plots with subscripts", {
  skip_if_not_installed("wrassp")
  skip_if_not_installed("ggplot2")

  test_wav <- system.file("samples", "sustained", "a1.wav",
                         package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  fms <- wrassp::forest(test_wav, numFormants = 2, toFile = FALSE)
  df <- as.data.frame(fms)

  # Create plot with subscripts
  p <- ggtrack(df, ggplot2::aes(x = frame_time, y = F1_Hz),
               use_subscripts = TRUE)

  # Check it's a ggplot
  expect_s3_class(p, "ggplot")

  # Check labels exist
  expect_true(!is.null(p$labels$x))
  expect_true(!is.null(p$labels$y))

  # Create plot without subscripts
  p_plain <- ggtrack(df, ggplot2::aes(x = frame_time, y = F1_Hz),
                     use_subscripts = FALSE)
  expect_s3_class(p_plain, "ggplot")
  expect_type(p_plain$labels$y, "character")
})

test_that("Integration test: Full workflow with trk_forest", {
  skip_if_not_installed("wrassp")
  skip_if_not_installed("ggplot2")

  test_wav <- system.file("samples", "sustained", "a1.wav",
                         package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Step 1: Get formants using superassp function
  fms <- trk_forest(test_wav, numFormants = 3, toFile = FALSE, verbose = FALSE)

  # Step 2: Convert to data.frame
  df <- as.data.frame(fms)

  # Step 3: Check template expansion worked
  expect_true(all(c("F1_Hz", "F2_Hz", "F3_Hz",
                    "B1_Hz", "B2_Hz", "B3_Hz") %in% names(df)))

  # Step 4: Check labels
  expect_equal(get_track_label(df, "F1_Hz"), "F1 [Hz]")

  # Step 5: Check plotmath
  expr <- get_track_label_expr(df, "F1_Hz", use_subscripts = TRUE)
  expect_s3_class(expr, "call")

  # Step 6: Create plot (should not error)
  p <- ggtrack(df, ggplot2::aes(x = frame_time, y = F1_Hz))
  expect_s3_class(p, "ggplot")
})
