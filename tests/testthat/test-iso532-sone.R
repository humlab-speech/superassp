# Tests for ISO 532 sone (loudness) conversions

test_that("phon_to_sone works at reference value (40 phon = 1 sone)", {
  skip_if_not_installed("superassp")

  # Reference value: 1 sone = 40 phons by definition
  expect_equal(phon_to_sone(40), 1.0, tolerance = 0.01)
  expect_equal(phon_to_sone(40, method = "moore-glasberg"), 1.0, tolerance = 0.01)
})

test_that("sone_to_phon works at reference value (1 sone = 40 phon)", {
  skip_if_not_installed("superassp")

  # Reference value: 1 sone = 40 phons by definition
  expect_equal(sone_to_phon(1), 40, tolerance = 0.01)
  expect_equal(sone_to_phon(1, method = "moore-glasberg"), 40, tolerance = 0.01)
})

test_that("doubling property: +10 phon ≈ 2× sone", {
  skip_if_not_installed("superassp")

  # Each 10 phon increase should approximately double loudness
  sone_40 <- phon_to_sone(40)
  sone_50 <- phon_to_sone(50)
  sone_60 <- phon_to_sone(60)

  expect_equal(sone_50 / sone_40, 2.0, tolerance = 0.01)
  expect_equal(sone_60 / sone_50, 2.0, tolerance = 0.01)
})

test_that("doubling property: 2× sone ≈ +10 phon", {
  skip_if_not_installed("superassp")

  # Doubling loudness should add approximately 10 phon
  phon_1 <- sone_to_phon(1)
  phon_2 <- sone_to_phon(2)
  phon_4 <- sone_to_phon(4)

  expect_equal(phon_2 - phon_1, 10, tolerance = 0.01)
  expect_equal(phon_4 - phon_2, 10, tolerance = 0.01)
})

test_that("round-trip conversion preserves values (Zwicker)", {
  skip_if_not_installed("superassp")

  test_phons <- c(10, 20, 30, 40, 50, 60, 70, 80, 90)

  for (phon_orig in test_phons) {
    sone <- phon_to_sone(phon_orig, method = "zwicker")
    phon_recovered <- sone_to_phon(sone, method = "zwicker")
    expect_equal(phon_recovered, phon_orig, tolerance = 0.01)
  }
})

test_that("round-trip conversion preserves values (Moore-Glasberg)", {
  skip_if_not_installed("superassp")

  test_phons <- c(10, 20, 30, 40, 50, 60, 70, 80, 90)

  for (phon_orig in test_phons) {
    sone <- phon_to_sone(phon_orig, method = "moore-glasberg")
    phon_recovered <- sone_to_phon(sone, method = "moore-glasberg")
    expect_equal(phon_recovered, phon_orig, tolerance = 0.1)
  }
})

test_that("round-trip conversion from sone (Zwicker)", {
  skip_if_not_installed("superassp")

  test_sones <- c(0.1, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0)

  for (sone_orig in test_sones) {
    phon <- sone_to_phon(sone_orig, method = "zwicker")
    sone_recovered <- phon_to_sone(phon, method = "zwicker")
    expect_equal(sone_recovered, sone_orig, tolerance = 0.01)
  }
})

test_that("phon_to_sone is monotonically increasing", {
  skip_if_not_installed("superassp")

  phons <- seq(0, 100, by = 5)
  sones <- phon_to_sone(phons)

  # Check that sone values increase monotonically
  expect_true(all(diff(sones) > 0))
})

test_that("sone_to_phon is monotonically increasing", {
  skip_if_not_installed("superassp")

  sones <- 10^seq(-2, 2, by = 0.2)  # 0.01 to 100 sones
  phons <- sone_to_phon(sones)

  # Check that phon values increase monotonically
  expect_true(all(diff(phons) > 0))
})

test_that("vectorization works for phon_to_sone", {
  skip_if_not_installed("superassp")

  phons <- c(20, 40, 60, 80)
  result <- phon_to_sone(phons)

  expect_length(result, 4)
  expect_true(all(is.numeric(result)))
  expect_true(all(!is.na(result)))
  expect_true(all(result > 0))
})

test_that("vectorization works for sone_to_phon", {
  skip_if_not_installed("superassp")

  sones <- c(0.5, 1.0, 2.0, 4.0)
  result <- sone_to_phon(sones)

  expect_length(result, 4)
  expect_true(all(is.numeric(result)))
  expect_true(all(!is.na(result)))
  expect_true(all(result > 0))
})

test_that("phon_to_sone handles low values correctly", {
  skip_if_not_installed("superassp")

  # Below 40 phon uses different formula
  sone_0 <- phon_to_sone(0)
  sone_10 <- phon_to_sone(10)
  sone_20 <- phon_to_sone(20)
  sone_30 <- phon_to_sone(30)

  # All should be below 1 sone
  expect_true(sone_0 < 1)
  expect_true(sone_10 < 1)
  expect_true(sone_20 < 1)
  expect_true(sone_30 < 1)

  # Should be monotonically increasing
  expect_true(sone_10 > sone_0)
  expect_true(sone_20 > sone_10)
  expect_true(sone_30 > sone_20)
})

test_that("sone_to_phon handles low values correctly", {
  skip_if_not_installed("superassp")

  # Below 1 sone uses different formula
  phon_0.1 <- sone_to_phon(0.1)
  phon_0.3 <- sone_to_phon(0.3)
  phon_0.5 <- sone_to_phon(0.5)
  phon_0.8 <- sone_to_phon(0.8)

  # All should be below 40 phon
  expect_true(phon_0.1 < 40)
  expect_true(phon_0.3 < 40)
  expect_true(phon_0.5 < 40)
  expect_true(phon_0.8 < 40)

  # Should be monotonically increasing
  expect_true(phon_0.3 > phon_0.1)
  expect_true(phon_0.5 > phon_0.3)
  expect_true(phon_0.8 > phon_0.5)
})

test_that("input validation for phon_to_sone", {
  skip_if_not_installed("superassp")

  # Non-numeric input
  expect_error(phon_to_sone("40"), "must be numeric")

  # NA values
  expect_error(phon_to_sone(c(40, NA, 60)), "contains NA")

  # Negative values
  expect_error(phon_to_sone(-10), "must be non-negative")
})

test_that("input validation for sone_to_phon", {
  skip_if_not_installed("superassp")

  # Non-numeric input
  expect_error(sone_to_phon("1"), "must be numeric")

  # NA values
  expect_error(sone_to_phon(c(1, NA, 2)), "contains NA")

  # Zero or negative values
  expect_error(sone_to_phon(0), "must be positive")
  expect_error(sone_to_phon(-1), "must be positive")
})

test_that("warning for extreme phon values", {
  skip_if_not_installed("superassp")

  # Above 120 phon should warn
  expect_warning(phon_to_sone(130), "exceed 120")
})

test_that("warning for extreme sone values", {
  skip_if_not_installed("superassp")

  # Above ~340 sones should warn
  expect_warning(sone_to_phon(500), "exceed 340")
})

test_that("Zwicker and Moore-Glasberg methods give similar results", {
  skip_if_not_installed("superassp")

  # Test at several points
  test_phons <- c(20, 40, 60, 80, 100)

  for (phon in test_phons) {
    sone_zwicker <- phon_to_sone(phon, method = "zwicker")
    sone_mg <- phon_to_sone(phon, method = "moore-glasberg")

    # Should be similar (within ~30% at low levels, ~10% at higher levels)
    # The methods differ more at low loudness levels
    if (phon < 30) {
      expect_equal(sone_zwicker, sone_mg, tolerance = 0.3 * max(sone_zwicker, sone_mg))
    } else {
      expect_equal(sone_zwicker, sone_mg, tolerance = 0.1 * max(sone_zwicker, sone_mg))
    }
  }
})

test_that("Moore-Glasberg table values are correct", {
  skip_if_not_installed("superassp")

  # Check a few reference points from the table
  ref_table <- superassp:::.moore_glasberg_table

  # 0 phon = 0.001 sone
  idx_0 <- which(ref_table$phon == 0)
  expect_equal(ref_table$sone[idx_0], 0.001)

  # 40 phon = 1.00 sone (reference point, same as Zwicker)
  idx_40 <- which(ref_table$phon == 40)
  expect_equal(ref_table$sone[idx_40], 1.00)

  # 120 phon = 337.60 sone
  idx_120 <- which(ref_table$phon == 120)
  expect_equal(ref_table$sone[idx_120], 337.60)
})

test_that("Moore-Glasberg interpolation works", {
  skip_if_not_installed("superassp")

  # Test interpolation between table values
  # Table has 40 and 45 phon
  sone_40 <- phon_to_sone(40, method = "moore-glasberg")
  sone_42.5 <- phon_to_sone(42.5, method = "moore-glasberg")
  sone_45 <- phon_to_sone(45, method = "moore-glasberg")

  # Interpolated value should be between endpoints
  expect_true(sone_42.5 > sone_40)
  expect_true(sone_42.5 < sone_45)
})

test_that("db_and_hz_to_sone combines conversions correctly", {
  skip_if_not_installed("superassp")

  # At 1000 Hz reference
  sone_40db <- db_and_hz_to_sone(40, 1000)
  expect_equal(sone_40db, 1.0, tolerance = 0.01)

  # 50 dB at 1 kHz should be ~2 sones (+10 dB ≈ 2× loudness)
  sone_50db <- db_and_hz_to_sone(50, 1000)
  expect_equal(sone_50db, 2.0, tolerance = 0.1)

  # Compare to step-by-step conversion
  phon <- db_and_hz_to_phon(60, 2000)
  sone_manual <- phon_to_sone(phon)
  sone_direct <- db_and_hz_to_sone(60, 2000)
  expect_equal(sone_direct, sone_manual, tolerance = 0.01)
})

test_that("sone_and_hz_to_db combines conversions correctly", {
  skip_if_not_installed("superassp")

  # At 1000 Hz reference
  db_1sone <- sone_and_hz_to_db(1, 1000)
  expect_equal(db_1sone, 40, tolerance = 0.1)

  # 2 sones at 1 kHz should be ~50 dB (2× loudness ≈ +10 dB)
  db_2sone <- sone_and_hz_to_db(2, 1000)
  expect_equal(db_2sone, 50, tolerance = 0.5)

  # Compare to step-by-step conversion
  phon <- sone_to_phon(3)
  db_manual <- phon_and_hz_to_db(phon, 500)
  db_direct <- sone_and_hz_to_db(3, 500)
  expect_equal(db_direct, db_manual, tolerance = 0.01)
})

test_that("round-trip with frequency conversions", {
  skip_if_not_installed("superassp")

  # Test: dB/Hz → sone → dB/Hz
  freqs <- c(100, 500, 1000, 2000, 4000)
  spls <- c(60, 55, 50, 52, 48)

  for (i in seq_along(freqs)) {
    sone <- db_and_hz_to_sone(spls[i], freqs[i])
    spl_recovered <- sone_and_hz_to_db(sone, freqs[i])
    expect_equal(spl_recovered, spls[i], tolerance = 0.5)
  }
})

test_that("frequency affects loudness perception", {
  skip_if_not_installed("superassp")

  # Same SPL at different frequencies gives different loudness
  sone_100hz <- db_and_hz_to_sone(60, 100)
  sone_1000hz <- db_and_hz_to_sone(60, 1000)
  sone_4000hz <- db_and_hz_to_sone(60, 4000)

  # Low frequency (100 Hz) less loud than 1 kHz
  expect_true(sone_100hz < sone_1000hz)

  # High frequency (4 kHz) more loud than 1 kHz (peak sensitivity)
  expect_true(sone_4000hz > sone_1000hz)
})

test_that("same loudness requires different SPL at different frequencies", {
  skip_if_not_installed("superassp")

  # 2 sones at different frequencies
  target_sone <- 2.0

  db_100hz <- sone_and_hz_to_db(target_sone, 100)
  db_1000hz <- sone_and_hz_to_db(target_sone, 1000)
  db_4000hz <- sone_and_hz_to_db(target_sone, 4000)

  # Low frequency needs more dB
  expect_true(db_100hz > db_1000hz)

  # High frequency needs less dB (more sensitive)
  expect_true(db_4000hz < db_1000hz)
})

test_that("piecewise formula transition is smooth at 40 phon", {
  skip_if_not_installed("superassp")

  # Test continuity at transition point (40 phon = 1 sone)
  phons_near_40 <- c(39.9, 40.0, 40.1)
  sones <- phon_to_sone(phons_near_40)

  # Check smoothness: differences should be small
  expect_true(abs(sones[2] - sones[1]) < 0.01)
  expect_true(abs(sones[3] - sones[2]) < 0.01)
})

test_that("piecewise formula transition is smooth at 1 sone", {
  skip_if_not_installed("superassp")

  # Test continuity at transition point (1 sone = 40 phon)
  sones_near_1 <- c(0.99, 1.0, 1.01)
  phons <- sone_to_phon(sones_near_1)

  # Check smoothness: differences should be small
  expect_true(abs(phons[2] - phons[1]) < 0.5)
  expect_true(abs(phons[3] - phons[2]) < 0.5)
})

test_that("conversion respects Stevens' power law", {
  skip_if_not_installed("superassp")

  # Stevens' power law: loudness ∝ intensity^0.3
  # For phon ≥ 40: sone = 2^((phon-40)/10)
  # This gives ~2× loudness for 10× intensity (10 dB)

  phon_40 <- 40
  phon_50 <- 50  # 10 dB more = 10× intensity

  sone_40 <- phon_to_sone(phon_40)
  sone_50 <- phon_to_sone(phon_50)

  # 10× intensity should give 10^0.3 ≈ 2× loudness
  intensity_ratio <- 10
  expected_loudness_ratio <- intensity_ratio^0.3
  actual_loudness_ratio <- sone_50 / sone_40

  expect_equal(actual_loudness_ratio, expected_loudness_ratio, tolerance = 0.05)
})
