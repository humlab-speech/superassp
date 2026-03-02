test_that("lst_covarep_vq works with single file (basic)", {
  # C++ implementation — no Python dependency

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- lst_covarep_vq(test_wav, verbose = FALSE)

  expect_type(result, "list")
  expect_true("glottal_flow_max" %in% names(result))
  expect_true("glottal_flow_min" %in% names(result))
  expect_true("glottal_derivative_peak" %in% names(result))
  expect_true("NAQ" %in% names(result))
  expect_true("QOQ" %in% names(result))
  expect_true("H1_H2" %in% names(result))
  expect_true("HRF" %in% names(result))
  expect_true("PSP" %in% names(result))

  # Values should be numeric
  expect_type(result$glottal_flow_max, "double")
  expect_type(result$HRF, "double")
  expect_type(result$PSP, "double")

  # Without F0 and GCI, some params should be NA
  expect_true(is.na(result$NAQ))  # No GCI
  expect_true(is.na(result$QOQ))  # No GCI
  expect_true(is.na(result$H1_H2))  # No F0

  # But these should always be computed
  expect_false(is.na(result$glottal_flow_max))
  expect_false(is.na(result$HRF))
  expect_false(is.na(result$PSP))
})

test_that("lst_covarep_vq works with scalar F0", {
  # C++ implementation — no Python dependency

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Provide F0
  f0_value <- 150.0

  result <- lst_covarep_vq(test_wav, f0 = f0_value, verbose = FALSE)

  expect_type(result, "list")
  expect_true("H1_H2" %in% names(result))

  # H1-H2 should now be computed (not NA)
  expect_false(is.na(result$H1_H2))
  expect_type(result$H1_H2, "double")

  # NAQ/QOQ still NA (no GCI)
  expect_true(is.na(result$NAQ))
  expect_true(is.na(result$QOQ))
})

test_that("lst_covarep_vq works with F0 vector", {
  # C++ implementation — no Python dependency
  skip("COVAREP SRH requires Python")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Get F0 contour from SRH
  f0_data <- trk_covarep_srh(test_wav, toFile = FALSE, verbose = FALSE)
  f0_vector <- f0_data$`F0[Hz]`[, 1]

  result <- lst_covarep_vq(test_wav, f0 = f0_vector, verbose = FALSE)

  expect_type(result, "list")
  expect_false(is.na(result$H1_H2))  # Should be computed from median
})

test_that("lst_covarep_vq handles batch processing", {
  # C++ implementation — no Python dependency

  test_files <- list.files(
    system.file("samples", "sustained", package = "superassp"),
    pattern = "\\.wav$",
    full.names = TRUE
  )

  skip_if(length(test_files) < 2, "Need at least 2 test files")
  test_files <- test_files[1:2]

  results <- lst_covarep_vq(test_files, verbose = FALSE)

  expect_type(results, "list")
  expect_length(results, 2)
  expect_type(results[[1]], "list")
  expect_type(results[[2]], "list")

  # Both should have all parameters
  expect_true("HRF" %in% names(results[[1]]))
  expect_true("HRF" %in% names(results[[2]]))
})

test_that("lst_covarep_vq supports time windowing", {
  # C++ implementation — no Python dependency

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result_full <- lst_covarep_vq(test_wav, verbose = FALSE)
  result_window <- lst_covarep_vq(test_wav,
                                  beginTime = 0.5,
                                  endTime = 1.5,
                                  verbose = FALSE)

  expect_type(result_full, "list")
  expect_type(result_window, "list")

  # Both should have same parameters
  expect_equal(names(result_full), names(result_window))

  # Values may differ due to different analysis window
  expect_type(result_window$HRF, "double")
})

test_that("lst_covarep_vq handles missing files gracefully", {
  # C++ implementation — no Python dependency

  expect_warning(
    result <- lst_covarep_vq("nonexistent.wav", verbose = FALSE),
    "File not found"
  )

  expect_null(result)
})

test_that("lst_covarep_vq validates F0 parameter", {
  # C++ implementation — no Python dependency

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Invalid F0 (negative)
  expect_error(
    lst_covarep_vq(test_wav, f0 = -100, verbose = FALSE),
    "f0 values must be non-negative"
  )

  # Invalid F0 (non-numeric)
  expect_error(
    lst_covarep_vq(test_wav, f0 = "invalid", verbose = FALSE),
    "f0 must be numeric"
  )
})

test_that("lst_covarep_vq validates GCI parameter", {
  # C++ implementation — no Python dependency

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Invalid GCI (non-numeric)
  expect_error(
    lst_covarep_vq(test_wav, gci = "invalid", verbose = FALSE),
    "gci must be numeric"
  )
})

test_that("lst_covarep_vq provides consistent results", {
  # C++ implementation — no Python dependency

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result1 <- lst_covarep_vq(test_wav, verbose = FALSE)
  result2 <- lst_covarep_vq(test_wav, verbose = FALSE)

  expect_equal(result1$HRF, result2$HRF)
  expect_equal(result1$PSP, result2$PSP)
  expect_equal(result1$glottal_flow_max, result2$glottal_flow_max)
})

test_that("lst_covarep_vq returns correct structure", {
  # C++ implementation — no Python dependency

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- lst_covarep_vq(test_wav, verbose = FALSE)

  # Should have exactly 8 parameters
  expected_params <- c("glottal_flow_max", "glottal_flow_min",
                      "glottal_derivative_peak", "NAQ", "QOQ",
                      "H1_H2", "HRF", "PSP")

  expect_length(result, length(expected_params))
  expect_true(all(expected_params %in% names(result)))

  # All values should be numeric (including NA)
  for (param in names(result)) {
    expect_type(result[[param]], "double")
  }
})

test_that("lst_covarep_vq with F0 and GCI integration", {
  skip("COVAREP SRH requires Python")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Get F0 from SRH
  f0_data <- trk_covarep_srh(test_wav, toFile = FALSE, verbose = FALSE)
  f0_mean <- median(f0_data$`F0[Hz]`[f0_data$VUV == 1, 1])

  # Test with F0
  result_with_f0 <- lst_covarep_vq(test_wav, f0 = f0_mean, verbose = FALSE)
  result_without_f0 <- lst_covarep_vq(test_wav, verbose = FALSE)

  # H1-H2 should differ
  expect_false(is.na(result_with_f0$H1_H2))
  expect_true(is.na(result_without_f0$H1_H2))

  # Other always-computed params should be the same
  expect_equal(result_with_f0$HRF, result_without_f0$HRF, tolerance = 1e-10)
  expect_equal(result_with_f0$PSP, result_without_f0$PSP, tolerance = 1e-10)
})

test_that("lst_covarep_vq batch with F0", {
  # C++ implementation — no Python dependency

  test_files <- list.files(
    system.file("samples", "sustained", package = "superassp"),
    pattern = "\\.wav$",
    full.names = TRUE
  )

  skip_if(length(test_files) < 2, "Need at least 2 test files")
  test_files <- test_files[1:2]

  # With scalar F0
  results <- lst_covarep_vq(test_files, f0 = 150, verbose = FALSE)

  expect_length(results, 2)
  expect_false(is.na(results[[1]]$H1_H2))
  expect_false(is.na(results[[2]]$H1_H2))
})
