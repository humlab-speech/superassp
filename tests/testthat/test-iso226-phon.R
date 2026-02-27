# Tests for ISO 226:2023 phon (loudness level) conversions

test_that("ucnv_db_and_hz_to_phon works at 1000 Hz reference", {
  skip_if_not_installed("superassp")

  # At 1000 Hz, phon = dB by definition
  expect_equal(ucnv_db_and_hz_to_phon(40, 1000), 40, tolerance = 0.1)
  expect_equal(ucnv_db_and_hz_to_phon(60, 1000), 60, tolerance = 0.1)
  expect_equal(ucnv_db_and_hz_to_phon(80, 1000), 80, tolerance = 0.1)
})

test_that("ucnv_phon_and_hz_to_db works at 1000 Hz reference", {
  skip_if_not_installed("superassp")

  # At 1000 Hz, dB = phon by definition
  expect_equal(ucnv_phon_and_hz_to_db(40, 1000), 40, tolerance = 0.1)
  expect_equal(ucnv_phon_and_hz_to_db(60, 1000), 60, tolerance = 0.1)
  expect_equal(ucnv_phon_and_hz_to_db(80, 1000), 80, tolerance = 0.1)
})

test_that("round-trip conversion preserves values at 1000 Hz", {
  skip_if_not_installed("superassp")

  test_spls <- c(30, 40, 50, 60, 70, 80)

  for (spl in test_spls) {
    phon <- ucnv_db_and_hz_to_phon(spl, 1000)
    spl_recovered <- ucnv_phon_and_hz_to_db(phon, 1000)
    expect_equal(spl_recovered, spl, tolerance = 0.1)
  }
})

test_that("round-trip conversion works across frequencies", {
  skip_if_not_installed("superassp")

  test_freqs <- c(100, 500, 1000, 2000, 4000, 8000)
  test_spls <- c(40, 50, 60, 70)

  for (freq in test_freqs) {
    for (spl in test_spls) {
      phon <- ucnv_db_and_hz_to_phon(spl, freq)
      spl_recovered <- ucnv_phon_and_hz_to_db(phon, freq)
      expect_equal(spl_recovered, spl, tolerance = 0.5)
    }
  }
})

test_that("low frequencies require higher SPL for same loudness", {
  skip_if_not_installed("superassp")

  # 40 phon should require more dB at low frequencies
  target_phon <- 40

  spl_100hz <- ucnv_phon_and_hz_to_db(target_phon, 100)
  spl_1000hz <- ucnv_phon_and_hz_to_db(target_phon, 1000)

  expect_true(spl_100hz > spl_1000hz)
  expect_equal(spl_1000hz, target_phon, tolerance = 0.1)  # Reference
})

test_that("high frequencies may require less SPL for same loudness", {
  skip_if_not_installed("superassp")

  # 40 phon at different frequencies
  target_phon <- 40

  spl_1000hz <- ucnv_phon_and_hz_to_db(target_phon, 1000)
  spl_4000hz <- ucnv_phon_and_hz_to_db(target_phon, 4000)

  # At moderate loudness levels, 4000 Hz is most sensitive
  expect_true(spl_4000hz <= spl_1000hz)
})

test_that("vectorized input works for ucnv_db_and_hz_to_phon", {
  skip_if_not_installed("superassp")

  spls <- c(40, 50, 60, 70)
  freqs <- c(100, 500, 1000, 4000)

  result <- ucnv_db_and_hz_to_phon(spls, freqs)

  expect_length(result, 4)
  expect_true(all(is.numeric(result)))
  expect_true(all(!is.na(result)))
})

test_that("vectorized input works for ucnv_phon_and_hz_to_db", {
  skip_if_not_installed("superassp")

  phons <- c(40, 50, 60, 70)
  freqs <- c(100, 500, 1000, 4000)

  result <- ucnv_phon_and_hz_to_db(phons, freqs)

  expect_length(result, 4)
  expect_true(all(is.numeric(result)))
  expect_true(all(!is.na(result)))
})

test_that("scalar recycling works", {
  skip_if_not_installed("superassp")

  # Single SPL, multiple frequencies
  result1 <- ucnv_db_and_hz_to_phon(60, c(100, 500, 1000, 4000))
  expect_length(result1, 4)

  # Multiple SPLs, single frequency
  result2 <- ucnv_db_and_hz_to_phon(c(40, 50, 60, 70), 1000)
  expect_length(result2, 4)

  # Single phon, multiple frequencies
  result3 <- ucnv_phon_and_hz_to_db(40, c(100, 500, 1000, 4000))
  expect_length(result3, 4)

  # Multiple phons, single frequency
  result4 <- ucnv_phon_and_hz_to_db(c(40, 50, 60, 70), 1000)
  expect_length(result4, 4)
})

test_that("standard 1/3-octave frequencies work", {
  skip_if_not_installed("superassp")

  # All standard frequencies from ISO 226 Table 1
  standard_freqs <- c(20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250,
                      315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500,
                      3150, 4000, 5000, 6300, 8000, 10000, 12500)

  result <- ucnv_phon_and_hz_to_db(40, standard_freqs)

  expect_length(result, length(standard_freqs))
  expect_true(all(is.numeric(result)))
  expect_true(all(!is.na(result)))
})

test_that("interpolation works for intermediate frequencies", {
  skip_if_not_installed("superassp")

  # Test frequency between 100 Hz and 125 Hz
  spl_100 <- ucnv_phon_and_hz_to_db(40, 100)
  spl_110 <- ucnv_phon_and_hz_to_db(40, 110)
  spl_125 <- ucnv_phon_and_hz_to_db(40, 125)

  # Interpolated value should be between the two endpoints
  expect_true(spl_110 >= min(spl_100, spl_125) && spl_110 <= max(spl_100, spl_125))
})

test_that("frequency range validation works", {
  skip_if_not_installed("superassp")

  # Below minimum frequency
  expect_error(ucnv_db_and_hz_to_phon(60, 10), "between 20 Hz and 12500 Hz")
  expect_error(ucnv_phon_and_hz_to_db(40, 10), "between 20 Hz and 12500 Hz")

  # Above maximum frequency
  expect_error(ucnv_db_and_hz_to_phon(60, 20000), "between 20 Hz and 12500 Hz")
  expect_error(ucnv_phon_and_hz_to_db(40, 20000), "between 20 Hz and 12500 Hz")
})

test_that("input validation works", {
  skip_if_not_installed("superassp")

  # Non-numeric input
  expect_error(ucnv_db_and_hz_to_phon("60", 1000), "must be numeric")
  expect_error(ucnv_db_and_hz_to_phon(60, "1000"), "must be numeric")
  expect_error(ucnv_phon_and_hz_to_db("40", 1000), "must be numeric")
  expect_error(ucnv_phon_and_hz_to_db(40, "1000"), "must be numeric")

  # Mismatched vector lengths
  expect_error(ucnv_db_and_hz_to_phon(c(40, 50), c(100, 500, 1000)),
               "same length")
  expect_error(ucnv_phon_and_hz_to_db(c(40, 50), c(100, 500, 1000)),
               "same length")
})

test_that("warning issued for low phon values", {
  skip_if_not_installed("superassp")

  # Below 20 phon should warn (near threshold)
  expect_warning(ucnv_db_and_hz_to_phon(10, 1000), "below 20 phon")
  expect_warning(ucnv_phon_and_hz_to_db(15, 1000), "below 20 phon")
})

test_that("warning issued for high phon values", {
  skip_if_not_installed("superassp")

  # Above 90 phon (20-4000 Hz) should warn
  expect_warning(ucnv_phon_and_hz_to_db(95, 1000), "exceeds reliable range")

  # Above 80 phon (5000-12500 Hz) should warn
  expect_warning(ucnv_phon_and_hz_to_db(85, 8000), "exceeds reliable range")
})

test_that("equal-loudness contours are monotonic in expected regions", {
  skip_if_not_installed("superassp")

  # For a fixed loudness level, SPL should generally increase
  # as frequency decreases below ~1000 Hz
  target_phon <- 60
  freqs_low <- c(100, 200, 500, 1000)

  spls <- ucnv_phon_and_hz_to_db(target_phon, freqs_low)

  # SPL should decrease as we go from 100 Hz to 1000 Hz
  expect_true(spls[1] > spls[4])  # 100 Hz > 1000 Hz
})

test_that("conversion matches ISO 226:2003 approximately", {
  skip_if_not_installed("superassp")

  # Test a few known values from ISO 226:2003
  # (2023 version differs by up to 0.6 dB)

  # 40 phon contour at various frequencies (approximate values)
  # These are ballpark figures for validation

  spl_100 <- ucnv_phon_and_hz_to_db(40, 100)
  expect_true(spl_100 > 50 && spl_100 < 65)  # Should be around 57 dB

  spl_1000 <- ucnv_phon_and_hz_to_db(40, 1000)
  expect_equal(spl_1000, 40, tolerance = 0.1)  # Exactly 40 dB

  spl_4000 <- ucnv_phon_and_hz_to_db(40, 4000)
  expect_true(spl_4000 > 30 && spl_4000 < 45)  # Should be around 37 dB
})

test_that("threshold region matches ISO 389-7", {
  skip_if_not_installed("superassp")

  # At very low phon levels, we approach the hearing threshold
  # ISO 226 aligns with ISO 389-7 for threshold values

  # The T_f values in the parameter table are the thresholds
  # At 1000 Hz, threshold is 2.4 dB
  # So very low phon should give SPL near 2.4 dB

  # Note: Below 20 phon is informative only
  suppressWarnings({
    spl_low <- ucnv_phon_and_hz_to_db(5, 1000)
  })

  # Should be close to threshold value (2.4 dB)
  expect_true(spl_low > 0 && spl_low < 15)
})

test_that("equal-loudness contour has expected U-shape", {
  skip_if_not_installed("superassp")

  # Equal-loudness contours typically have a U-shape
  # with minimum sensitivity around 3000-4000 Hz

  target_phon <- 60
  test_freqs <- c(100, 500, 1000, 2000, 4000, 8000)

  spls <- ucnv_phon_and_hz_to_db(target_phon, test_freqs)

  # Find minimum SPL (maximum sensitivity)
  min_idx <- which.min(spls)

  # Should be in the 2000-4000 Hz range
  expect_true(test_freqs[min_idx] >= 2000 && test_freqs[min_idx] <= 4000)
})

test_that("ISO 226 parameter table is correctly loaded", {
  skip_if_not_installed("superassp")

  # Check parameter table structure
  expect_true(exists(".iso226_params", where = "package:superassp", mode = "list"))

  params <- superassp:::.iso226_params

  expect_equal(nrow(params), 29)  # 29 standard frequencies
  expect_equal(ncol(params), 4)   # freq_hz, alpha_f, L_U, T_f

  # Check 1000 Hz reference values
  ref_row <- params[params$freq_hz == 1000, ]
  expect_equal(ref_row$alpha_f, 0.300)
  expect_equal(ref_row$L_U, 0.0)
  expect_equal(ref_row$T_f, 2.4)
})

test_that("parameter interpolation is reasonable", {
  skip_if_not_installed("superassp")

  # Get parameters for exact frequency
  params_100 <- superassp:::.get_iso226_params(100)

  # Get parameters for interpolated frequency
  params_110 <- superassp:::.get_iso226_params(110)

  # Get parameters for next exact frequency
  params_125 <- superassp:::.get_iso226_params(125)

  # Interpolated values should be between the two endpoints
  expect_true(params_110$alpha_f >= min(params_100$alpha_f, params_125$alpha_f))
  expect_true(params_110$alpha_f <= max(params_100$alpha_f, params_125$alpha_f))
})
