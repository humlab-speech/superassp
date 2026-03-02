# GFM-IAIF is now C++ — always available
gfmiaif_available <- function() TRUE

test_that("trk_gfmiaif works with single file", {
  skip_if_not(gfmiaif_available(), "GFM-IAIF dependencies not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- trk_gfmiaif(test_wav, toFile = FALSE, verbose = FALSE)

  expect_s3_class(result, "AsspDataObj")

  # Check that vocal tract tracks exist (av_0 to av_48 by default)
  expect_true("av_0" %in% names(result))
  expect_true("av_48" %in% names(result))

  # Check that glottis tracks exist (ag_0 to ag_3 by default)
  expect_true("ag_0" %in% names(result))
  expect_true("ag_3" %in% names(result))

  # Check that lip radiation tracks exist (al_0, al_1)
  expect_true("al_0" %in% names(result))
  expect_true("al_1" %in% names(result))

  # Check dimensions
  expect_true(nrow(result$av_0) > 0)
  expect_equal(ncol(result$av_0), 1)  # Each track is single column

  # All tracks should have same number of frames
  expect_equal(nrow(result$av_0), nrow(result$ag_0))
  expect_equal(nrow(result$ag_0), nrow(result$al_0))

  # Check attributes
  expect_true(!is.null(attr(result, "sampleRate")))
  expect_true(!is.null(attr(result, "startTime")))
  expect_true(!is.null(attr(result, "origFreq")))
})

test_that("trk_gfmiaif respects custom vocal tract order", {
  skip_if_not(gfmiaif_available(), "GFM-IAIF dependencies not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- trk_gfmiaif(test_wav,
                       nv = 24,
                       toFile = FALSE,
                       verbose = FALSE)

  expect_s3_class(result, "AsspDataObj")

  # Should have av_0 to av_24
  expect_true("av_0" %in% names(result))
  expect_true("av_24" %in% names(result))
  expect_false("av_48" %in% names(result))  # Should not have higher order

  # Glottis and lip radiation should be unaffected
  expect_true("ag_3" %in% names(result))
  expect_true("al_1" %in% names(result))
})

test_that("trk_gfmiaif respects custom glottis order", {
  skip_if_not(gfmiaif_available(), "GFM-IAIF dependencies not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Note: This tests functionality, but ng=3 is highly recommended
  result <- trk_gfmiaif(test_wav,
                       ng = 5,
                       toFile = FALSE,
                       verbose = FALSE)

  expect_s3_class(result, "AsspDataObj")

  # Should have ag_0 to ag_5
  expect_true("ag_0" %in% names(result))
  expect_true("ag_5" %in% names(result))
})

test_that("trk_gfmiaif respects leaky integration coefficient", {
  skip_if_not(gfmiaif_available(), "GFM-IAIF dependencies not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result_095 <- trk_gfmiaif(test_wav, d = 0.95,
                           toFile = FALSE, verbose = FALSE)
  result_099 <- trk_gfmiaif(test_wav, d = 0.99,
                           toFile = FALSE, verbose = FALSE)

  expect_s3_class(result_095, "AsspDataObj")
  expect_s3_class(result_099, "AsspDataObj")

  # Lip radiation coefficients should differ
  expect_false(identical(result_095$al_1, result_099$al_1))

  # Vocal tract and glottis should also differ (affected by pre-emphasis)
  expect_false(identical(result_095$av_1, result_099$av_1))
})

test_that("trk_gfmiaif supports time windowing", {
  skip_if_not(gfmiaif_available(), "GFM-IAIF dependencies not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result_full <- trk_gfmiaif(test_wav, toFile = FALSE, verbose = FALSE)
  result_window <- trk_gfmiaif(test_wav,
                              beginTime = 0.1,
                              endTime = 0.3,
                              toFile = FALSE,
                              verbose = FALSE)

  expect_s3_class(result_full, "AsspDataObj")
  expect_s3_class(result_window, "AsspDataObj")

  # Windowed version should have fewer frames
  expect_true(nrow(result_window$av_0) < nrow(result_full$av_0))

  # Start time should reflect windowing (use tolerance for floating point)
  expect_true(attr(result_window, "startTime") >= 0.1 - 1e-3)
})

test_that("trk_gfmiaif supports different window types", {
  skip_if_not(gfmiaif_available(), "GFM-IAIF dependencies not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result_hann <- trk_gfmiaif(test_wav, window = "HANN",
                            toFile = FALSE, verbose = FALSE)
  result_hamming <- trk_gfmiaif(test_wav, window = "HAMMING",
                               toFile = FALSE, verbose = FALSE)
  result_blackman <- trk_gfmiaif(test_wav, window = "BLACKMAN",
                                toFile = FALSE, verbose = FALSE)

  expect_s3_class(result_hann, "AsspDataObj")
  expect_s3_class(result_hamming, "AsspDataObj")
  expect_s3_class(result_blackman, "AsspDataObj")

  # Results should differ slightly
  expect_false(identical(result_hann$av_1, result_hamming$av_1))
  expect_false(identical(result_hann$av_1, result_blackman$av_1))
})

test_that("trk_gfmiaif supports custom window parameters", {
  skip_if_not(gfmiaif_available(), "GFM-IAIF dependencies not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result_small <- trk_gfmiaif(test_wav,
                             windowShift = 5.0,
                             windowSize = 20.0,
                             toFile = FALSE,
                             verbose = FALSE)

  result_large <- trk_gfmiaif(test_wav,
                             windowShift = 20.0,
                             windowSize = 40.0,
                             toFile = FALSE,
                             verbose = FALSE)

  expect_s3_class(result_small, "AsspDataObj")
  expect_s3_class(result_large, "AsspDataObj")

  # Smaller shift = more frames
  expect_true(nrow(result_small$av_0) > nrow(result_large$av_0))

  # Check sample rates differ
  expect_true(abs(attr(result_small, "sampleRate") -
                 attr(result_large, "sampleRate")) > 1)
})

test_that("trk_gfmiaif validates parameter ranges", {
  skip_if_not(gfmiaif_available(), "GFM-IAIF dependencies not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # nv out of range
  expect_error(trk_gfmiaif(test_wav, nv = 0, toFile = FALSE, verbose = FALSE))
  expect_error(trk_gfmiaif(test_wav, nv = 200, toFile = FALSE, verbose = FALSE))

  # ng out of range
  expect_error(trk_gfmiaif(test_wav, ng = 0, toFile = FALSE, verbose = FALSE))
  expect_error(trk_gfmiaif(test_wav, ng = 20, toFile = FALSE, verbose = FALSE))

  # d out of range
  expect_error(trk_gfmiaif(test_wav, d = 0.5, toFile = FALSE, verbose = FALSE))
  expect_error(trk_gfmiaif(test_wav, d = 1.1, toFile = FALSE, verbose = FALSE))

  # Invalid window type
  expect_error(trk_gfmiaif(test_wav, window = "INVALID",
                          toFile = FALSE, verbose = FALSE))
})

test_that("trk_gfmiaif produces valid LP coefficients", {
  skip_if_not(gfmiaif_available(), "GFM-IAIF dependencies not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- trk_gfmiaif(test_wav, toFile = FALSE, verbose = FALSE)

  # First coefficient should always be 1.0 (LP polynomial convention)
  expect_equal(mean(abs(result$av_0 - 1.0)), 0, tolerance = 1e-10)
  expect_equal(mean(abs(result$ag_0 - 1.0)), 0, tolerance = 1e-10)
  expect_equal(mean(abs(result$al_0 - 1.0)), 0, tolerance = 1e-10)

  # Lip radiation second coefficient should be approximately -d
  # al = [1, -d]
  expect_true(all(result$al_1 < 0))  # Should be negative
  expect_true(all(abs(result$al_1) > 0.9))  # Should be close to -0.99
})

test_that("trk_gfmiaif can write to file", {
  skip_if_not(gfmiaif_available(), "GFM-IAIF dependencies not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Create temporary directory
  temp_dir <- tempdir()

  # Process with toFile = TRUE (use nv=12 to fit SSFF 1KB header limit)
  n_processed <- trk_gfmiaif(test_wav,
                            nv = 12,
                            outputDirectory = temp_dir,
                            explicitExt = "gfm",
                            toFile = TRUE,
                            verbose = FALSE)

  expect_equal(n_processed, 1)

  # Check output file exists
  output_file <- file.path(temp_dir, sub("\\.[^.]+$", ".gfm", basename(test_wav)))
  expect_true(file.exists(output_file))

  # Should be able to read it back
  result <- superassp:::read.AsspDataObj(output_file)
  expect_s3_class(result, "AsspDataObj")
  expect_true("av_0" %in% names(result))

  # Clean up
  unlink(output_file)
})

test_that("trk_gfmiaif batch processing works", {
  skip_if_not(gfmiaif_available(), "GFM-IAIF dependencies not installed")

  # Get multiple test files
  test_files <- list.files(
    system.file("samples", "sustained", package = "superassp"),
    pattern = "\\.wav$",
    full.names = TRUE
  )

  skip_if(length(test_files) < 2, "Not enough test files")

  # Take first 2 files
  test_files <- test_files[1:2]

  temp_dir <- tempdir()

  n_processed <- trk_gfmiaif(test_files,
                            nv = 12,
                            outputDirectory = temp_dir,
                            toFile = TRUE,
                            verbose = FALSE)

  expect_equal(n_processed, 2)

  # Check both files exist
  for (f in test_files) {
    output_file <- file.path(temp_dir, sub("\\.[^.]+$", ".gfm", basename(f)))
    expect_true(file.exists(output_file))
    unlink(output_file)
  }
})

test_that("trk_gfmiaif fails gracefully with invalid input", {
  skip_if_not(gfmiaif_available(), "GFM-IAIF dependencies not installed")

  # Non-existent file
  expect_error(trk_gfmiaif("nonexistent.wav", toFile = FALSE, verbose = FALSE))

  # Multiple files with toFile = FALSE
  test_files <- c("file1.wav", "file2.wav")
  expect_error(trk_gfmiaif(test_files, toFile = FALSE, verbose = FALSE),
               "toFile=FALSE only permitted for single files")
})

test_that("trk_gfmiaif centerTime parameter works", {
  skip_if_not(gfmiaif_available(), "GFM-IAIF dependencies not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result_center <- trk_gfmiaif(test_wav, centerTime = TRUE,
                              toFile = FALSE, verbose = FALSE)
  result_start <- trk_gfmiaif(test_wav, centerTime = FALSE,
                             toFile = FALSE, verbose = FALSE)

  expect_s3_class(result_center, "AsspDataObj")
  expect_s3_class(result_start, "AsspDataObj")

  # Start times should differ
  expect_true(attr(result_center, "startTime") != attr(result_start, "startTime"))
})
