test_that("trk_covarep_iaif works with single file", {
  # C++ implementation — no Python dependency

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- trk_covarep_iaif(test_wav, toFile = FALSE, verbose = FALSE)

  expect_s3_class(result, "AsspDataObj")
  expect_true("glottal_flow" %in% names(result))
  expect_true("glottal_derivative" %in% names(result))

  # Check dimensions
  expect_true(nrow(result$glottal_flow) > 0)
  expect_equal(ncol(result$glottal_flow), 1)
  expect_equal(nrow(result$glottal_flow), nrow(result$glottal_derivative))

  # Check attributes
  expect_true(!is.null(attr(result, "sampleRate")))
  expect_true(!is.null(attr(result, "startTime")))
  expect_equal(attr(result, "startTime"), 0.0)
})

test_that("trk_covarep_iaif respects custom filter orders", {
  # C++ implementation — no Python dependency

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- trk_covarep_iaif(test_wav,
                             order_vt = 16,
                             order_gl = 3,
                             toFile = FALSE,
                             verbose = FALSE)

  expect_s3_class(result, "AsspDataObj")
  expect_true(nrow(result$glottal_flow) > 0)
})

test_that("trk_covarep_iaif works with auto filter orders", {
  # C++ implementation — no Python dependency

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # NULL should trigger automatic order selection
  result <- trk_covarep_iaif(test_wav,
                             order_vt = NULL,
                             order_gl = NULL,
                             toFile = FALSE,
                             verbose = FALSE)

  expect_s3_class(result, "AsspDataObj")
  expect_true(nrow(result$glottal_flow) > 0)
})

test_that("trk_covarep_iaif respects leaky coefficient", {
  # C++ implementation — no Python dependency

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result_095 <- trk_covarep_iaif(test_wav, leaky_coef = 0.95,
                                 toFile = FALSE, verbose = FALSE)
  result_099 <- trk_covarep_iaif(test_wav, leaky_coef = 0.99,
                                 toFile = FALSE, verbose = FALSE)

  expect_s3_class(result_095, "AsspDataObj")
  expect_s3_class(result_099, "AsspDataObj")

  # Results should differ
  expect_false(identical(result_095$glottal_flow, result_099$glottal_flow))
})

test_that("trk_covarep_iaif works with and without HPF", {
  # C++ implementation — no Python dependency

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result_hpf <- trk_covarep_iaif(test_wav, hpfilt = TRUE,
                                 toFile = FALSE, verbose = FALSE)
  result_nohpf <- trk_covarep_iaif(test_wav, hpfilt = FALSE,
                                   toFile = FALSE, verbose = FALSE)

  expect_s3_class(result_hpf, "AsspDataObj")
  expect_s3_class(result_nohpf, "AsspDataObj")

  # Results should differ
  expect_false(identical(result_hpf$glottal_flow, result_nohpf$glottal_flow))
})

test_that("trk_covarep_iaif supports time windowing", {
  # C++ implementation — no Python dependency

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result_full <- trk_covarep_iaif(test_wav, toFile = FALSE, verbose = FALSE)
  result_window <- trk_covarep_iaif(test_wav,
                                    beginTime = 0.5,
                                    endTime = 1.5,
                                    toFile = FALSE,
                                    verbose = FALSE)

  expect_s3_class(result_window, "AsspDataObj")
  expect_true(nrow(result_window$glottal_flow) < nrow(result_full$glottal_flow))
  expect_equal(attr(result_window, "startTime"), 0.5)
})

test_that("trk_covarep_iaif handles batch processing", {
  # C++ implementation — no Python dependency

  test_files <- list.files(
    system.file("samples", "sustained", package = "superassp"),
    pattern = "\\.wav$",
    full.names = TRUE
  )

  skip_if(length(test_files) < 2, "Need at least 2 test files")
  test_files <- test_files[1:2]

  results <- trk_covarep_iaif(test_files, toFile = FALSE, verbose = FALSE)

  expect_type(results, "list")
  expect_length(results, 2)
  expect_s3_class(results[[1]], "AsspDataObj")
  expect_s3_class(results[[2]], "AsspDataObj")
})

test_that("trk_covarep_iaif writes files when toFile=TRUE", {
  # C++ implementation — no Python dependency

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  tmp_dir <- tempdir()
  result <- trk_covarep_iaif(test_wav,
                             toFile = TRUE,
                             outputDirectory = tmp_dir,
                             verbose = FALSE)

  expect_type(result, "character")
  expect_true(file.exists(result))
  expect_match(result, "\\.glf$")

  # Clean up
  unlink(result)
})

test_that("trk_covarep_iaif handles missing files gracefully", {
  # C++ implementation — no Python dependency

  expect_warning(
    result <- trk_covarep_iaif("nonexistent.wav", toFile = FALSE, verbose = FALSE),
    "File not found"
  )

  expect_null(result)
})

test_that("trk_covarep_iaif validates parameters", {
  # C++ implementation — no Python dependency

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Invalid leaky coefficient
  expect_error(
    trk_covarep_iaif(test_wav, leaky_coef = 1.5, toFile = FALSE, verbose = FALSE),
    "leaky_coef must be between 0 and 1"
  )

  expect_error(
    trk_covarep_iaif(test_wav, leaky_coef = 0, toFile = FALSE, verbose = FALSE),
    "leaky_coef must be between 0 and 1"
  )
})

test_that("trk_covarep_iaif provides consistent results", {
  # C++ implementation — no Python dependency

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result1 <- trk_covarep_iaif(test_wav, toFile = FALSE, verbose = FALSE)
  result2 <- trk_covarep_iaif(test_wav, toFile = FALSE, verbose = FALSE)

  expect_equal(nrow(result1$glottal_flow), nrow(result2$glottal_flow))
  expect_equal(result1$glottal_flow, result2$glottal_flow)
  expect_equal(result1$glottal_derivative, result2$glottal_derivative)
})

test_that("trk_covarep_iaif output dimensions match input", {
  # C++ implementation — no Python dependency

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Load audio to check expected length
  audio_bin <- av::read_audio_bin(test_wav, channels = 1)
  expected_samples <- length(audio_bin)

  result <- trk_covarep_iaif(test_wav, toFile = FALSE, verbose = FALSE)

  # IAIF may trim edges slightly, but should be close
  expect_true(nrow(result$glottal_flow) > 0.9 * expected_samples)
  expect_true(nrow(result$glottal_flow) <= expected_samples)
})

test_that("trk_covarep_iaif function attributes are set", {
  expect_equal(attr(trk_covarep_iaif, "ext"), "glf")
  expect_equal(attr(trk_covarep_iaif, "tracks"), c("glottal_flow", "glottal_derivative"))
  expect_equal(attr(trk_covarep_iaif, "outputType"), "SSFF")
})
