# Tests for automatic unit conversion in as.data.frame.AsspDataObj and as_tibble.AsspDataObj

test_that(".parse_unit_from_colname extracts units correctly", {
  skip_if_not_installed("superassp")

  # Test with valid unit labels
  expect_equal(.parse_unit_from_colname("fo[Hz]"), "Hz")
  expect_equal(.parse_unit_from_colname("intensity[dB]"), "dB")
  expect_equal(.parse_unit_from_colname("fm[kHz]"), "kHz")
  expect_equal(.parse_unit_from_colname("pitch[semitone]"), "semitone")

  # Test without unit labels
  expect_true(is.na(.parse_unit_from_colname("frame_time")))
  expect_true(is.na(.parse_unit_from_colname("fo")))
  expect_true(is.na(.parse_unit_from_colname("fm_1")))

  # Test with brackets not at end
  expect_true(is.na(.parse_unit_from_colname("fo[Hz]_1")))
})

test_that("as.data.frame.AsspDataObj converts units automatically", {
  skip_if_not_installed("superassp")
  skip_if_not_installed("units")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Load F0 data
  f0_obj <- fo(test_wav, toFile = FALSE, verbose = FALSE)

  # Convert to data frame with unit conversion
  df <- as.data.frame(f0_obj, convert_units = TRUE)

  # Check that the fo[Hz] column has units
  expect_true(inherits(df[["fo[Hz]"]], "units"))

  # Check that the unit is Hz
  unit_str <- as.character(units::deparse_unit(df[["fo[Hz]"]]))
  expect_equal(unit_str, "Hz")

  # Check that frame_time is NOT converted (no unit label)
  expect_false(inherits(df$frame_time, "units"))
})

test_that("as.data.frame.AsspDataObj can disable unit conversion", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Load F0 data
  f0_obj <- fo(test_wav, toFile = FALSE, verbose = FALSE)

  # Convert to data frame WITHOUT unit conversion
  df <- as.data.frame(f0_obj, convert_units = FALSE)

  # Check that the fo[Hz] column does NOT have units
  expect_false(inherits(df[["fo[Hz]"]], "units"))
  expect_true(is.numeric(df[["fo[Hz]"]]))
})

test_that("as_tibble.AsspDataObj converts units automatically", {
  skip_if_not_installed("superassp")
  skip_if_not_installed("units")
  skip_if_not_installed("tibble")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Load F0 data
  f0_obj <- fo(test_wav, toFile = FALSE, verbose = FALSE)

  # Convert to tibble with unit conversion
  tbl <- as_tibble(f0_obj, convert_units = TRUE)

  # Check that it's a tibble
  expect_s3_class(tbl, "tbl_df")

  # Check that the fo[Hz] column has units
  expect_true(inherits(tbl[["fo[Hz]"]], "units"))

  # Check that the unit is Hz
  unit_str <- as.character(units::deparse_unit(tbl[["fo[Hz]"]]))
  expect_equal(unit_str, "Hz")

  # Check that time columns are NOT converted (no unit labels)
  expect_false(inherits(tbl$times_orig, "units"))
  expect_false(inherits(tbl$times_rel, "units"))
  expect_false(inherits(tbl$times_norm, "units"))
})

test_that("as_tibble.AsspDataObj can disable unit conversion", {
  skip_if_not_installed("superassp")
  skip_if_not_installed("tibble")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Load F0 data
  f0_obj <- fo(test_wav, toFile = FALSE, verbose = FALSE)

  # Convert to tibble WITHOUT unit conversion
  tbl <- as_tibble(f0_obj, convert_units = FALSE)

  # Check that the fo[Hz] column does NOT have units
  expect_false(inherits(tbl[["fo[Hz]"]], "units"))
  expect_true(is.numeric(tbl[["fo[Hz]"]]))
})

test_that("unit conversion works with psychoacoustic units", {
  skip_if_not_installed("superassp")
  skip_if_not_installed("units")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Load F0 data and convert to Bark scale
  f0_obj <- fo(test_wav, toFile = FALSE, verbose = FALSE)
  df <- as.data.frame(f0_obj, convert_units = TRUE)

  # Extract Hz values
  hz_values <- as.numeric(df[["fo[Hz]"]])

  # Convert to Bark using our psychoacoustic scale functions
  bark_values <- hz_to_bark(hz_values[hz_values > 0], as_units = TRUE)

  # Check that conversion worked
  expect_true(inherits(bark_values, "units"))
  unit_str <- as.character(units::deparse_unit(bark_values))
  expect_equal(unit_str, "Bark")
})

test_that("unit conversion handles invalid units gracefully", {
  skip_if_not_installed("superassp")
  skip_if_not_installed("units")

  # Create a mock AsspDataObj with an invalid unit label
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  f0_obj <- fo(test_wav, toFile = FALSE, verbose = FALSE)

  # Rename track to use an invalid unit
  names(f0_obj)[1] <- "pitch[InvalidUnit]"

  # Try to convert - should warn but not error
  expect_warning(
    df <- as.data.frame(f0_obj, convert_units = TRUE),
    "Could not convert column"
  )

  # Check that column is still numeric (not converted)
  expect_false(inherits(df[["pitch[InvalidUnit]"]], "units"))
  expect_true(is.numeric(df[["pitch[InvalidUnit]"]]))
})

test_that("unit conversion preserves data values", {
  skip_if_not_installed("superassp")
  skip_if_not_installed("units")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Load F0 data
  f0_obj <- fo(test_wav, toFile = FALSE, verbose = FALSE)

  # Convert with and without units
  df_with_units <- as.data.frame(f0_obj, convert_units = TRUE)
  df_without_units <- as.data.frame(f0_obj, convert_units = FALSE)

  # Extract numeric values
  with_units_numeric <- as.numeric(df_with_units[["fo[Hz]"]])
  without_units_numeric <- df_without_units[["fo[Hz]"]]

  # Values should be identical
  expect_equal(with_units_numeric, without_units_numeric)
})

test_that("unit conversion works with multi-column tracks", {
  skip_if_not_installed("superassp")
  skip_if_not_installed("units")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Create a mock multi-column track by loading F0 and duplicating
  f0_obj <- fo(test_wav, toFile = FALSE, verbose = FALSE)

  # Manually add a second column to simulate multi-field tracks
  # (In real use, formant tracks have multiple columns like fm[Hz]_1, fm[Hz]_2, etc.)
  original_track <- f0_obj[[1]]
  two_col_track <- cbind(original_track, original_track)
  f0_obj[[1]] <- two_col_track

  # Convert to data frame with name separator
  df <- as.data.frame(f0_obj, name.separator = "_", convert_units = TRUE)

  # Multi-column tracks get numeric suffixes appended AFTER the unit label
  # So the pattern becomes "fo[Hz]_1" which does NOT match our regex (unit not at end)
  # This is expected behavior - the suffix prevents unit detection
  expect_false(inherits(df[["fo[Hz]_1"]], "units"))
  expect_false(inherits(df[["fo[Hz]_2"]], "units"))

  # They should still be numeric
  expect_true(is.numeric(df[["fo[Hz]_1"]]))
  expect_true(is.numeric(df[["fo[Hz]_2"]]))
})

test_that("as_tibble field selection works with unit conversion", {
  skip_if_not_installed("superassp")
  skip_if_not_installed("units")
  skip_if_not_installed("tibble")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Load F0 data
  f0_obj <- fo(test_wav, toFile = FALSE, verbose = FALSE)

  # Select only the first field
  tbl <- as_tibble(f0_obj, field = 1, convert_units = TRUE)

  # Check structure
  expect_s3_class(tbl, "tbl_df")
  expect_true("fo[Hz]" %in% names(tbl))

  # Check unit conversion
  expect_true(inherits(tbl[["fo[Hz]"]], "units"))
})

test_that("na.zeros works with unit conversion", {
  skip_if_not_installed("superassp")
  skip_if_not_installed("units")
  skip_if_not_installed("tibble")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Load F0 data
  f0_obj <- fo(test_wav, toFile = FALSE, verbose = FALSE)

  # Convert with na.zeros = TRUE (default)
  tbl_with_na <- as_tibble(f0_obj, convert_units = TRUE, na.zeros = TRUE)

  # Convert with na.zeros = FALSE
  tbl_without_na <- as_tibble(f0_obj, convert_units = TRUE, na.zeros = FALSE)

  # Count zeros vs NAs
  zeros_with_na <- sum(as.numeric(tbl_with_na[["fo[Hz]"]]) == 0, na.rm = TRUE)
  zeros_without_na <- sum(as.numeric(tbl_without_na[["fo[Hz]"]]) == 0, na.rm = TRUE)

  # Should have fewer zeros with na.zeros = TRUE
  expect_true(zeros_with_na <= zeros_without_na)
})
