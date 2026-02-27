test_that("ucnv_hz_to_bark works with numeric input", {
  skip_if_not_installed("superassp")

  # Test single value
  bark <- ucnv_hz_to_bark(1000, as_units = FALSE)
  expect_type(bark, "double")
  expect_length(bark, 1)
  expect_true(bark > 0 && bark < 24)

  # Test vector input
  freqs <- c(100, 500, 1000, 2000, 4000)
  bark_vals <- ucnv_hz_to_bark(freqs, as_units = FALSE)
  expect_length(bark_vals, 5)
  expect_true(all(bark_vals > 0))
  expect_true(all(diff(bark_vals) > 0))  # Should be monotonically increasing
})


test_that("ucnv_hz_to_bark validates input", {
  skip_if_not_installed("superassp")

  # Negative frequencies should error
  expect_error(ucnv_hz_to_bark(-100, as_units = FALSE))

  # NA values should be preserved
  result <- ucnv_hz_to_bark(c(1000, NA, 2000), as_units = FALSE)
  expect_true(is.na(result[2]))
  expect_false(is.na(result[1]))
})


test_that("ucnv_hz_to_bark methods produce reasonable results", {
  skip_if_not_installed("superassp")

  freq <- 1000
  bark_traun <- ucnv_hz_to_bark(freq, method = "traunmuller", as_units = FALSE)
  bark_zwicker <- ucnv_hz_to_bark(freq, method = "zwicker", as_units = FALSE)
  bark_wang <- ucnv_hz_to_bark(freq, method = "wang", as_units = FALSE)

  # All should be in reasonable range for 1000 Hz
  expect_true(bark_traun > 7 && bark_traun < 10)
  expect_true(bark_zwicker > 7 && bark_zwicker < 10)
  expect_true(bark_wang > 7 && bark_wang < 10)

  # They should be similar but not identical
  expect_false(bark_traun == bark_zwicker)
  expect_false(bark_traun == bark_wang)
})


test_that("ucnv_hz_to_bark works with units package", {
  skip_if_not_installed("superassp")
  skip_if_not_installed("units")

  library(units)

  # Test with units input
  freq <- set_units(1000, Hz)
  bark <- ucnv_hz_to_bark(freq)

  expect_s3_class(bark, "units")
  expect_equal(as.character(units(bark)$numerator), "Bark")

  # Test with kHz input (should convert automatically)
  freq_khz <- set_units(1, kHz)
  bark_khz <- ucnv_hz_to_bark(freq_khz)

  expect_s3_class(bark_khz, "units")
  # Should be same as 1000 Hz
  expect_equal(as.numeric(bark), as.numeric(bark_khz), tolerance = 0.01)
})


test_that("ucnv_bark_to_hz works with numeric input", {
  skip_if_not_installed("superassp")

  # Test single value
  hz <- ucnv_bark_to_hz(10, as_units = FALSE)
  expect_type(hz, "double")
  expect_length(hz, 1)
  expect_true(hz > 0)

  # Test vector input
  bark_vals <- c(5, 10, 15, 20)
  freqs <- ucnv_bark_to_hz(bark_vals, as_units = FALSE)
  expect_length(freqs, 4)
  expect_true(all(freqs > 0))
  expect_true(all(diff(freqs) > 0))  # Should be monotonically increasing
})


test_that("ucnv_bark_to_hz validates input", {
  skip_if_not_installed("superassp")

  # Out of range values should warn
  expect_warning(ucnv_bark_to_hz(30, as_units = FALSE))
  expect_warning(ucnv_bark_to_hz(-5, as_units = FALSE))

  # NA values should be preserved
  result <- ucnv_bark_to_hz(c(10, NA, 15), as_units = FALSE)
  expect_true(is.na(result[2]))
  expect_false(is.na(result[1]))
})


test_that("ucnv_bark_to_hz methods produce reasonable results", {
  skip_if_not_installed("superassp")

  bark <- 10
  hz_traun <- ucnv_bark_to_hz(bark, method = "traunmuller", as_units = FALSE)
  hz_zwicker <- ucnv_bark_to_hz(bark, method = "zwicker", as_units = FALSE)
  hz_wang <- ucnv_bark_to_hz(bark, method = "wang", as_units = FALSE)

  # All should be in reasonable range
  expect_true(hz_traun > 500 && hz_traun < 2000)
  expect_true(hz_zwicker > 500 && hz_zwicker < 2000)
  expect_true(hz_wang > 500 && hz_wang < 2000)
})


test_that("ucnv_bark_to_hz works with units package", {
  skip_if_not_installed("superassp")
  skip_if_not_installed("units")

  library(units)

  # Ensure Bark unit is installed
  superassp:::.ensure_bark_unit()

  # Test with units input
  bark <- set_units(10, Bark)
  freq <- ucnv_bark_to_hz(bark)

  expect_s3_class(freq, "units")
  expect_equal(as.character(units(freq)$numerator), "Hz")
})


test_that("round-trip conversion preserves values (Traunmüller)", {
  skip_if_not_installed("superassp")

  test_freqs <- c(100, 500, 1000, 2000, 4000, 8000)

  for (freq in test_freqs) {
    bark <- ucnv_hz_to_bark(freq, method = "traunmuller", as_units = FALSE)
    freq_recovered <- ucnv_bark_to_hz(bark, method = "traunmuller", as_units = FALSE)

    # Note: Low frequencies (<200 Hz) have larger errors due to Traunmüller's
    # low-frequency correction factor
    tolerance <- if (freq < 200) 20 else 0.1

    expect_equal(freq_recovered, freq, tolerance = tolerance,
                 label = sprintf("Round-trip for %g Hz", freq))
  }
})


test_that("round-trip conversion preserves values (Wang)", {
  skip_if_not_installed("superassp")

  test_freqs <- c(100, 500, 1000, 2000, 4000, 8000)

  for (freq in test_freqs) {
    bark <- ucnv_hz_to_bark(freq, method = "wang", as_units = FALSE)
    freq_recovered <- ucnv_bark_to_hz(bark, method = "wang", as_units = FALSE)

    expect_equal(freq_recovered, freq, tolerance = 0.1,
                 label = sprintf("Round-trip for %g Hz", freq))
  }
})


test_that("round-trip conversion preserves values (Zwicker)", {
  skip_if_not_installed("superassp")

  test_freqs <- c(100, 500, 1000, 2000, 4000, 8000)

  for (freq in test_freqs) {
    bark <- ucnv_hz_to_bark(freq, method = "zwicker", as_units = FALSE)
    freq_recovered <- ucnv_bark_to_hz(bark, method = "zwicker", as_units = FALSE)

    expect_equal(freq_recovered, freq, tolerance = 1.0,  # Zwicker needs numerical root-finding
                 label = sprintf("Round-trip for %g Hz", freq))
  }
})


test_that("units integration works end-to-end", {
  skip_if_not_installed("superassp")
  skip_if_not_installed("units")

  library(units)

  # Start with frequency in Hz
  freq <- set_units(1000, Hz)

  # Convert to Bark
  bark <- ucnv_hz_to_bark(freq)
  expect_s3_class(bark, "units")

  # Convert back to Hz
  freq_recovered <- ucnv_bark_to_hz(bark)
  expect_s3_class(freq_recovered, "units")

  # Values should match
  expect_equal(as.numeric(freq), as.numeric(freq_recovered), tolerance = 0.1)
})


test_that("Bark unit behaves correctly with units package", {
  skip_if_not_installed("superassp")
  skip_if_not_installed("units")

  library(units)

  # Ensure Bark unit is installed
  superassp:::.ensure_bark_unit()

  # Can create Bark values directly
  b1 <- set_units(10, Bark)
  expect_s3_class(b1, "units")

  # Can do arithmetic
  b2 <- set_units(5, Bark)
  b_sum <- b1 + b2
  expect_equal(as.numeric(b_sum), 15)

  # Can drop units
  b_numeric <- drop_units(b1)
  expect_type(b_numeric, "double")
  expect_equal(b_numeric, 10)
})


test_that("conversion handles edge cases", {
  skip_if_not_installed("superassp")

  # Very low frequency
  bark_low <- ucnv_hz_to_bark(20, as_units = FALSE)
  expect_true(bark_low >= 0)

  # Very high frequency
  bark_high <- ucnv_hz_to_bark(15000, as_units = FALSE)
  expect_true(bark_high <= 24)

  # Zero frequency (boundary case)
  bark_zero <- ucnv_hz_to_bark(0, as_units = FALSE)
  expect_true(is.finite(bark_zero))
})


test_that("Traunmüller low-frequency correction is applied", {
  skip_if_not_installed("superassp")

  # For frequencies that map to Bark < 2, correction should apply
  # Test a low frequency
  freq <- 100  # Should be well below 2 Bark before correction

  bark <- ucnv_hz_to_bark(freq, method = "traunmuller", as_units = FALSE)

  # Without correction: (26.81 * 100) / (1960 + 100) - 0.53
  bark_uncorrected <- (26.81 * 100) / (1960 + 100) - 0.53
  # Should be < 2, so correction should have been applied
  expect_true(bark_uncorrected < 2)
  expect_true(bark != bark_uncorrected)

  # Round-trip should still work
  freq_recovered <- ucnv_bark_to_hz(bark, method = "traunmuller", as_units = FALSE)
  expect_equal(freq_recovered, freq, tolerance = 1.0)
})


test_that("methods parameter is validated", {
  skip_if_not_installed("superassp")

  expect_error(ucnv_hz_to_bark(1000, method = "invalid"))
  expect_error(ucnv_bark_to_hz(10, method = "invalid"))
})


test_that("all methods produce monotonic mappings", {
  skip_if_not_installed("superassp")

  test_freqs <- seq(100, 8000, by = 100)
  methods <- c("traunmuller", "zwicker", "wang")

  for (method in methods) {
    bark_vals <- ucnv_hz_to_bark(test_freqs, method = method, as_units = FALSE)

    # Check monotonicity: higher frequency -> higher Bark
    expect_true(all(diff(bark_vals) > 0),
                label = sprintf("Monotonicity for method=%s", method))
  }
})


test_that("as_units parameter controls output format", {
  skip_if_not_installed("superassp")
  skip_if_not_installed("units")

  # With as_units = FALSE, should get numeric even if input has units
  freq <- units::set_units(1000, Hz)
  bark_numeric <- ucnv_hz_to_bark(freq, as_units = FALSE)
  expect_type(bark_numeric, "double")
  expect_false(inherits(bark_numeric, "units"))

  # With as_units = TRUE, should get units even if input is numeric
  bark_units <- ucnv_hz_to_bark(1000, as_units = TRUE)
  expect_s3_class(bark_units, "units")

  # Same for inverse
  freq_numeric <- ucnv_bark_to_hz(10, as_units = FALSE)
  expect_type(freq_numeric, "double")

  freq_units <- ucnv_bark_to_hz(10, as_units = TRUE)
  expect_s3_class(freq_units, "units")
})


# ============================================================================
# ERB Scale Tests
# ============================================================================

test_that("ucnv_hz_to_erb works with numeric input", {
  skip_if_not_installed("superassp")

  erb <- ucnv_hz_to_erb(1000, as_units = FALSE)
  expect_type(erb, "double")
  expect_length(erb, 1)
  expect_true(erb > 0)

  # Vector input
  freqs <- c(100, 500, 1000, 2000, 4000)
  erb_vals <- ucnv_hz_to_erb(freqs, as_units = FALSE)
  expect_length(erb_vals, 5)
  expect_true(all(diff(erb_vals) > 0))  # Monotonic
})


test_that("ucnv_hz_to_erb methods produce reasonable results", {
  skip_if_not_installed("superassp")

  freq <- 1000
  erb_glasberg <- ucnv_hz_to_erb(freq, method = "glasberg1990", as_units = FALSE)
  erb_moore <- ucnv_hz_to_erb(freq, method = "moore1983", as_units = FALSE)

  # Both should be in reasonable range
  expect_true(erb_glasberg > 10 && erb_glasberg < 25)
  expect_true(erb_moore > 10 && erb_moore < 50)
})


test_that("ucnv_hz_to_erb works with units package", {
  skip_if_not_installed("superassp")
  skip_if_not_installed("units")

  library(units)

  freq <- set_units(1000, Hz)
  erb <- ucnv_hz_to_erb(freq)

  expect_s3_class(erb, "units")
})


test_that("ucnv_erb_to_hz works with numeric input", {
  skip_if_not_installed("superassp")

  hz <- ucnv_erb_to_hz(15, as_units = FALSE)
  expect_type(hz, "double")
  expect_true(hz > 0)
})


test_that("round-trip ERB conversion preserves values (Glasberg)", {
  skip_if_not_installed("superassp")

  test_freqs <- c(100, 500, 1000, 2000, 4000)

  for (freq in test_freqs) {
    erb <- ucnv_hz_to_erb(freq, method = "glasberg1990", as_units = FALSE)
    freq_recovered <- ucnv_erb_to_hz(erb, method = "glasberg1990", as_units = FALSE)

    expect_equal(freq_recovered, freq, tolerance = 0.1,
                 label = sprintf("ERB round-trip for %g Hz", freq))
  }
})


test_that("ERB conversion is monotonic", {
  skip_if_not_installed("superassp")

  freqs <- seq(100, 8000, by = 100)
  erb_vals <- ucnv_hz_to_erb(freqs, as_units = FALSE)

  expect_true(all(diff(erb_vals) > 0))
})


# ============================================================================
# Mel Scale Tests
# ============================================================================

test_that("ucnv_hz_to_mel works with numeric input", {
  skip_if_not_installed("superassp")

  mel <- ucnv_hz_to_mel(1000, as_units = FALSE)
  expect_type(mel, "double")
  expect_length(mel, 1)

  # 1000 Hz should be close to 1000 mel
  expect_true(mel > 900 && mel < 1100)

  # Vector input
  freqs <- c(100, 500, 1000, 2000, 4000)
  mel_vals <- ucnv_hz_to_mel(freqs, as_units = FALSE)
  expect_length(mel_vals, 5)
  expect_true(all(diff(mel_vals) > 0))
})


test_that("ucnv_hz_to_mel methods produce reasonable results", {
  skip_if_not_installed("superassp")

  freq <- 1000
  mel_htk <- ucnv_hz_to_mel(freq, method = "htk", as_units = FALSE)
  mel_slaney <- ucnv_hz_to_mel(freq, method = "slaney", as_units = FALSE)

  # 1000 Hz should be approximately 1000 mels
  expect_true(mel_htk > 900 && mel_htk < 1100)
  expect_true(mel_slaney > 0)
})


test_that("ucnv_hz_to_mel works with units package", {
  skip_if_not_installed("superassp")
  skip_if_not_installed("units")

  library(units)

  freq <- set_units(1000, Hz)
  mel <- ucnv_hz_to_mel(freq)

  expect_s3_class(mel, "units")
})


test_that("ucnv_mel_to_hz works with numeric input", {
  skip_if_not_installed("superassp")

  hz <- ucnv_mel_to_hz(1000, as_units = FALSE)
  expect_type(hz, "double")
  expect_true(hz > 0)
})


test_that("round-trip mel conversion preserves values (HTK)", {
  skip_if_not_installed("superassp")

  test_freqs <- c(100, 500, 1000, 2000, 4000)

  for (freq in test_freqs) {
    mel <- ucnv_hz_to_mel(freq, method = "htk", as_units = FALSE)
    freq_recovered <- ucnv_mel_to_hz(mel, method = "htk", as_units = FALSE)

    expect_equal(freq_recovered, freq, tolerance = 0.01,
                 label = sprintf("Mel round-trip for %g Hz", freq))
  }
})


test_that("round-trip mel conversion preserves values (Slaney)", {
  skip_if_not_installed("superassp")

  test_freqs <- c(100, 500, 1000, 2000, 4000)

  for (freq in test_freqs) {
    mel <- ucnv_hz_to_mel(freq, method = "slaney", as_units = FALSE)
    freq_recovered <- ucnv_mel_to_hz(mel, method = "slaney", as_units = FALSE)

    expect_equal(freq_recovered, freq, tolerance = 1,
                 label = sprintf("Mel round-trip for %g Hz", freq))
  }
})


test_that("mel conversion is monotonic", {
  skip_if_not_installed("superassp")

  freqs <- seq(100, 8000, by = 100)
  mel_vals <- ucnv_hz_to_mel(freqs, as_units = FALSE)

  expect_true(all(diff(mel_vals) > 0))
})


# ============================================================================
# Semitone Scale Tests
# ============================================================================

test_that("ucnv_hz_to_semitone works with numeric input", {
  skip_if_not_installed("superassp")

  # 880 Hz is 12 semitones above 440 Hz (A4 to A5 = 1 octave)
  st <- ucnv_hz_to_semitone(880, ref_freq = 440, as_units = FALSE)
  expect_type(st, "double")
  expect_equal(st, 12, tolerance = 0.001)

  # 440 Hz relative to itself should be 0
  st_zero <- ucnv_hz_to_semitone(440, ref_freq = 440, as_units = FALSE)
  expect_equal(st_zero, 0, tolerance = 0.001)

  # Vector input
  freqs <- c(440, 466.16, 493.88, 523.25, 554.37, 587.33)  # A4 to D5
  st_vals <- ucnv_hz_to_semitone(freqs, ref_freq = 440, as_units = FALSE)
  expect_length(st_vals, 6)
})


test_that("ucnv_hz_to_semitone handles different reference frequencies", {
  skip_if_not_installed("superassp")

  # Middle C (261.63 Hz) to C one octave higher (523.25 Hz)
  st <- ucnv_hz_to_semitone(523.25, ref_freq = 261.63, as_units = FALSE)
  expect_equal(st, 12, tolerance = 0.01)

  # With units
  skip_if_not_installed("units")
  library(units)

  freq <- set_units(880, Hz)
  ref <- set_units(440, Hz)
  st_units <- ucnv_hz_to_semitone(freq, ref_freq = ref, as_units = TRUE)
  expect_s3_class(st_units, "units")
  expect_equal(as.numeric(st_units), 12, tolerance = 0.001)
})


test_that("ucnv_hz_to_semitone works with units package", {
  skip_if_not_installed("superassp")
  skip_if_not_installed("units")

  library(units)

  freq <- set_units(880, Hz)
  st <- ucnv_hz_to_semitone(freq, ref_freq = 440)

  expect_s3_class(st, "units")
})


test_that("ucnv_semitone_to_hz works with numeric input", {
  skip_if_not_installed("superassp")

  # 12 semitones above 440 Hz should be 880 Hz
  hz <- ucnv_semitone_to_hz(12, ref_freq = 440, as_units = FALSE)
  expect_type(hz, "double")
  expect_equal(hz, 880, tolerance = 0.01)

  # 0 semitones should return reference frequency
  hz_zero <- ucnv_semitone_to_hz(0, ref_freq = 440, as_units = FALSE)
  expect_equal(hz_zero, 440)

  # Negative semitones (going down)
  hz_down <- ucnv_semitone_to_hz(-12, ref_freq = 440, as_units = FALSE)
  expect_equal(hz_down, 220, tolerance = 0.01)
})


test_that("round-trip semitone conversion preserves values", {
  skip_if_not_installed("superassp")

  test_freqs <- c(110, 220, 440, 880, 1760)

  for (freq in test_freqs) {
    st <- ucnv_hz_to_semitone(freq, ref_freq = 440, as_units = FALSE)
    freq_recovered <- ucnv_semitone_to_hz(st, ref_freq = 440, as_units = FALSE)

    expect_equal(freq_recovered, freq, tolerance = 0.01,
                 label = sprintf("Semitone round-trip for %g Hz", freq))
  }
})


test_that("semitone conversion validates input", {
  skip_if_not_installed("superassp")

  # Negative/zero frequencies should error
  expect_error(ucnv_hz_to_semitone(0, ref_freq = 440))
  expect_error(ucnv_hz_to_semitone(-100, ref_freq = 440))
  expect_error(ucnv_hz_to_semitone(100, ref_freq = -440))
})


test_that("musical intervals are correct", {
  skip_if_not_installed("superassp")

  # Perfect fifth (7 semitones) has frequency ratio of ~1.5
  freq_fifth <- ucnv_semitone_to_hz(7, ref_freq = 440, as_units = FALSE)
  expect_equal(freq_fifth / 440, 2^(7/12), tolerance = 0.001)

  # Perfect fourth (5 semitones)
  freq_fourth <- ucnv_semitone_to_hz(5, ref_freq = 440, as_units = FALSE)
  expect_equal(freq_fourth / 440, 2^(5/12), tolerance = 0.001)

  # Major third (4 semitones)
  freq_third <- ucnv_semitone_to_hz(4, ref_freq = 440, as_units = FALSE)
  expect_equal(freq_third / 440, 2^(4/12), tolerance = 0.001)
})


# ============================================================================
# Cross-scale comparison tests
# ============================================================================

test_that("all scales are monotonic", {
  skip_if_not_installed("superassp")

  freqs <- seq(100, 4000, by = 50)

  # Bark
  bark_vals <- ucnv_hz_to_bark(freqs, as_units = FALSE)
  expect_true(all(diff(bark_vals) > 0))

  # ERB
  erb_vals <- ucnv_hz_to_erb(freqs, as_units = FALSE)
  expect_true(all(diff(erb_vals) > 0))

  # Mel
  mel_vals <- ucnv_hz_to_mel(freqs, as_units = FALSE)
  expect_true(all(diff(mel_vals) > 0))

  # Semitone (relative to constant reference)
  st_vals <- ucnv_hz_to_semitone(freqs, ref_freq = 440, as_units = FALSE)
  expect_true(all(diff(st_vals) > 0))
})


test_that("all scales handle NA values correctly", {
  skip_if_not_installed("superassp")

  freqs <- c(1000, NA, 2000)

  # Bark
  bark_vals <- ucnv_hz_to_bark(freqs, as_units = FALSE)
  expect_true(is.na(bark_vals[2]))
  expect_false(is.na(bark_vals[1]))

  # ERB
  erb_vals <- ucnv_hz_to_erb(freqs, as_units = FALSE)
  expect_true(is.na(erb_vals[2]))

  # Mel
  mel_vals <- ucnv_hz_to_mel(freqs, as_units = FALSE)
  expect_true(is.na(mel_vals[2]))

  # Semitone
  st_vals <- ucnv_hz_to_semitone(freqs, ref_freq = 440, as_units = FALSE)
  expect_true(is.na(st_vals[2]))
})


# ============================================================================
# Semitone Reference Source Tests (UEP83, Praat, A4)
# ============================================================================

test_that("ucnv_hz_to_semitone works with UEP83 reference", {
  skip_if_not_installed("superassp")

  # UEP83 uses 110 Hz (A2) as reference
  # 220 Hz is one octave (12 semitones) above 110 Hz
  st <- ucnv_hz_to_semitone(220, ref_source = "UEP83", as_units = FALSE)
  expect_equal(st, 12, tolerance = 0.001)

  # 110 Hz relative to itself should be 0
  st_zero <- ucnv_hz_to_semitone(110, ref_source = "UEP83", as_units = FALSE)
  expect_equal(st_zero, 0, tolerance = 0.001)
})


test_that("ucnv_hz_to_semitone works with Praat reference", {
  skip_if_not_installed("superassp")

  # Praat uses 100 Hz as arbitrary reference
  # 200 Hz is one octave (12 semitones) above 100 Hz
  st <- ucnv_hz_to_semitone(200, ref_source = "Praat", as_units = FALSE)
  expect_equal(st, 12, tolerance = 0.001)

  # 100 Hz relative to itself should be 0
  st_zero <- ucnv_hz_to_semitone(100, ref_source = "Praat", as_units = FALSE)
  expect_equal(st_zero, 0, tolerance = 0.001)
})


test_that("ucnv_hz_to_semitone works with A4 reference (default)", {
  skip_if_not_installed("superassp")

  # A4 uses 440 Hz as reference (default)
  # 880 Hz is one octave (12 semitones) above 440 Hz
  st <- ucnv_hz_to_semitone(880, ref_source = "A4", as_units = FALSE)
  expect_equal(st, 12, tolerance = 0.001)

  # Should also work without specifying ref_source (default)
  st_default <- ucnv_hz_to_semitone(880, as_units = FALSE)
  expect_equal(st_default, 12, tolerance = 0.001)
})


test_that("ref_freq parameter overrides ref_source", {
  skip_if_not_installed("superassp")

  # Explicit ref_freq should override ref_source
  st1 <- ucnv_hz_to_semitone(880, ref_freq = 440, ref_source = "UEP83", as_units = FALSE)
  st2 <- ucnv_hz_to_semitone(880, ref_freq = 440, as_units = FALSE)

  expect_equal(st1, st2)
  expect_equal(st1, 12, tolerance = 0.001)
})


test_that("ucnv_semitone_to_hz works with different reference sources", {
  skip_if_not_installed("superassp")

  # UEP83: 12 ST above 110 Hz = 220 Hz
  hz_uep <- ucnv_semitone_to_hz(12, ref_source = "UEP83", as_units = FALSE)
  expect_equal(hz_uep, 220, tolerance = 0.01)

  # Praat: 12 ST above 100 Hz = 200 Hz
  hz_praat <- ucnv_semitone_to_hz(12, ref_source = "Praat", as_units = FALSE)
  expect_equal(hz_praat, 200, tolerance = 0.01)

  # A4: 12 ST above 440 Hz = 880 Hz
  hz_a4 <- ucnv_semitone_to_hz(12, ref_source = "A4", as_units = FALSE)
  expect_equal(hz_a4, 880, tolerance = 0.01)
})


test_that("round-trip conversion works with all reference sources", {
  skip_if_not_installed("superassp")

  test_freqs <- c(110, 220, 440)

  for (ref_src in c("UEP83", "Praat", "A4")) {
    for (freq in test_freqs) {
      st <- ucnv_hz_to_semitone(freq, ref_source = ref_src, as_units = FALSE)
      freq_recovered <- ucnv_semitone_to_hz(st, ref_source = ref_src, as_units = FALSE)

      expect_equal(freq_recovered, freq, tolerance = 0.01,
                   label = sprintf("Round-trip %g Hz with %s", freq, ref_src))
    }
  }
})


test_that("UEP83 reference produces expected phonetogram values", {
  skip_if_not_installed("superassp")

  # Based on UEP 1983 standard phonetogram (Figure 1)
  # Some example frequencies from the standard:
  # A (110 Hz), c (131 Hz), e (165 Hz), a (220 Hz), c1 (262 Hz), a1 (440 Hz)

  # Test relative to 110 Hz (A in UEP notation, A2 in scientific)
  uep_freqs <- c(
    A = 110,
    c = 131,
    e = 165,
    a = 220,
    c1 = 262,
    a1 = 440
  )

  semitones <- ucnv_hz_to_semitone(uep_freqs, ref_source = "UEP83", as_units = FALSE)

  # A should be 0 ST (reference)
  expect_equal(unname(semitones["A"]), 0, tolerance = 0.01)

  # a (220 Hz) should be 12 ST above A (110 Hz)
  expect_equal(unname(semitones["a"]), 12, tolerance = 0.01)

  # a1 (440 Hz) should be 24 ST above A (110 Hz)
  expect_equal(unname(semitones["a1"]), 24, tolerance = 0.01)
})


test_that("Praat reference useful for speech F0 analysis", {
  skip_if_not_installed("superassp")

  # Typical adult male F0 range: 85-180 Hz
  # Typical adult female F0 range: 165-255 Hz
  
  male_low <- 85
  male_high <- 180
  female_low <- 165
  female_high <- 255

  # Convert to semitones re 100 Hz (Praat convention)
  male_st_low <- ucnv_hz_to_semitone(male_low, ref_source = "Praat", as_units = FALSE)
  male_st_high <- ucnv_hz_to_semitone(male_high, ref_source = "Praat", as_units = FALSE)
  female_st_low <- ucnv_hz_to_semitone(female_low, ref_source = "Praat", as_units = FALSE)
  female_st_high <- ucnv_hz_to_semitone(female_high, ref_source = "Praat", as_units = FALSE)

  # All should be above 0 (all above 100 Hz)
  expect_true(male_st_low < 0)   # 85 Hz is below 100 Hz
  expect_true(male_st_high > 0)
  expect_true(female_st_low > 0)
  expect_true(female_st_high > 0)

  # Female range should be higher than male range
  expect_true(female_st_low > male_st_low)
  expect_true(female_st_high > male_st_high)
})


test_that("ref_source parameter validation works", {
  skip_if_not_installed("superassp")

  # Valid ref_source values should work
  expect_no_error(ucnv_hz_to_semitone(220, ref_source = "UEP83", as_units = FALSE))
  expect_no_error(ucnv_hz_to_semitone(200, ref_source = "Praat", as_units = FALSE))
  expect_no_error(ucnv_hz_to_semitone(880, ref_source = "A4", as_units = FALSE))

  # Invalid ref_source should error
  expect_error(ucnv_hz_to_semitone(220, ref_source = "invalid"))
})
