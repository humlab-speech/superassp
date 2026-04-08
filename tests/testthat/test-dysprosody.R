# Tests for lst_dysprosody prosodic assessment

test_that("lst_dysprosody works with single file", {
  skip_if_not_installed("pladdrr")

  test_wav <- system.file("samples", "sustained", "a32b.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- lst_dysprosody(test_wav, verbose = FALSE)

  expect_type(result, "list")
  expect_true(length(result) > 0)

  expect_true("Duration" %in% names(result))
  expect_true("PitchMean" %in% names(result))
  expect_true("PitchRange" %in% names(result))
  expect_true("PitchKey" %in% names(result))
  expect_true("IntsIntLabels" %in% names(result))
  expect_true("UniqueIntsInt" %in% names(result))

  expect_equal(length(result), 193)

  expect_type(result$Duration, "double")
  expect_type(result$PitchMean, "double")

  expect_gt(result$Duration, 0)
  expect_lt(result$Duration, 10)
})

test_that("lst_dysprosody handles time windowing", {
  skip_if_not_installed("pladdrr")

  test_wav <- system.file("samples", "sustained", "a32b.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result_full <- lst_dysprosody(test_wav, verbose = FALSE)
  skip_if(is.null(result_full), "Full file too short")

  full_duration <- result_full$Duration

  result_windowed <- lst_dysprosody(
    test_wav,
    beginTime = 0.5,
    endTime = min(2.0, full_duration - 0.1),
    verbose = FALSE
  )

  if (!is.null(result_windowed)) {
    expect_type(result_windowed, "list")
    expect_true("Duration" %in% names(result_windowed))
    expect_lt(result_windowed$Duration, full_duration)
    expect_gt(result_windowed$Duration, 0)
  }
})

test_that("lst_dysprosody skips very short files", {
  skip_if_not_installed("pladdrr")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- lst_dysprosody(
    test_wav,
    beginTime = 0,
    endTime = 0.5,
    verbose = FALSE
  )

  expect_null(result)
})

test_that("lst_dysprosody handles batch processing (sequential)", {
  skip_if_not_installed("pladdrr")

  test_files <- list.files(
    system.file("samples", "sustained", package = "superassp"),
    pattern = "\\.wav$",
    full.names = TRUE
  )

  skip_if(length(test_files) < 2, "Not enough test files")

  files <- test_files[1:min(3, length(test_files))]

  results <- lst_dysprosody(
    files,
    verbose = FALSE,
    parallel = FALSE
  )

  expect_type(results, "list")
  expect_true(length(results) > 0)

  for (result in results) {
    if (!is.null(result)) {
      expect_type(result, "list")
      expect_true("Duration" %in% names(result))
      expect_true("PitchMean" %in% names(result))
    }
  }
})

test_that("lst_dysprosody handles missing files", {
  skip_if_not_installed("pladdrr")

  expect_error(
    lst_dysprosody("nonexistent_file.wav", verbose = FALSE)
  )
})

test_that("lst_dysprosody handles custom F0 range", {
  skip_if_not_installed("pladdrr")

  test_wav <- system.file("samples", "sustained", "a32b.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- lst_dysprosody(
    test_wav,
    minF = 100,
    maxF = 500,
    verbose = FALSE
  )

  if (!is.null(result)) {
    expect_type(result, "list")
    expect_true("PitchMean" %in% names(result))
  }
})

test_that("lst_dysprosody returns consistent feature names", {
  skip_if_not_installed("pladdrr")

  test_files <- list.files(
    system.file("samples", "sustained", package = "superassp"),
    pattern = "\\.wav$",
    full.names = TRUE
  )[1:2]

  skip_if(length(test_files) < 2, "Not enough test files")

  results <- lst_dysprosody(test_files, verbose = FALSE, parallel = FALSE)

  valid_results <- Filter(Negate(is.null), results)
  skip_if(length(valid_results) < 2, "Not enough valid results")

  names1 <- names(valid_results[[1]])
  names2 <- names(valid_results[[2]])

  expect_equal(length(names1), length(names2))
  expect_equal(sort(names1), sort(names2))
})
