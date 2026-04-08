# Tests for JSON Track Format (JSTF) functionality

test_that("create_json_track_obj works with list results", {
  skip_if_not_installed("superassp")
  
  results <- list(
    jitter = 85.3,
    shimmer = 4.2,
    hnr = 15.7
  )
  
  obj <- create_json_track_obj(
    results = results,
    function_name = "lst_test",
    file_path = "test.wav",
    sample_rate = 16000,
    audio_duration = 5.0,
    beginTime = 0.0,
    endTime = 5.0
  )
  
  expect_s3_class(obj, "JsonTrackObj")
  expect_equal(obj$format, "JSTF")
  expect_equal(obj$version, "1.0")
  expect_equal(obj$function_name, "lst_test")
  expect_equal(length(obj$field_schema), 3)
  expect_equal(length(obj$slices), 1)
  expect_equal(obj$slices[[1]]$begin_time, 0.0)
  expect_equal(obj$slices[[1]]$end_time, 5.0)
})

test_that("create_json_track_obj works with data.frame results", {
  skip_if_not_installed("superassp")
  
  results <- data.frame(
    measure1 = 123.4,
    measure2 = 567.8,
    measure3 = 901.2
  )
  
  obj <- create_json_track_obj(
    results = results,
    function_name = "lst_test",
    file_path = "test.wav",
    sample_rate = 16000,
    audio_duration = 3.0
  )
  
  expect_s3_class(obj, "JsonTrackObj")
  expect_equal(length(obj$field_schema), 3)
  expect_equal(names(obj$field_schema), c("measure1", "measure2", "measure3"))
})

test_that("validate_json_track catches invalid objects", {
  skip_if_not_installed("superassp")
  
  # Valid object
  valid_obj <- create_json_track_obj(
    list(x = 1), "test", "test.wav", 16000, 1.0, 0, 1
  )
  expect_true(validate_json_track(valid_obj))
  
  # Invalid: wrong format
  invalid1 <- valid_obj
  invalid1$format <- "INVALID"
  expect_error(validate_json_track(invalid1), "Invalid format")
  
  # Invalid: begin_time >= end_time
  invalid2 <- valid_obj
  invalid2$slices[[1]]$begin_time <- 2.0
  invalid2$slices[[1]]$end_time <- 1.0
  expect_error(validate_json_track(invalid2), "begin_time must be < end_time")
  
  # Invalid: values length mismatch
  invalid3 <- valid_obj
  invalid3$slices[[1]]$values <- c(1, 2, 3)  # schema has 1 field
  expect_error(validate_json_track(invalid3), "values length")
})

test_that("write_json_track and read_json_track round-trip", {
  skip_if_not_installed("superassp")
  skip_if_not_installed("jsonlite")
  
  # Create object
  obj1 <- create_json_track_obj(
    results = list(f0 = 150, intensity = 70),
    function_name = "lst_test",
    file_path = "audio.wav",
    sample_rate = 16000,
    audio_duration = 2.0,
    beginTime = 0.0,
    endTime = 2.0
  )
  
  # Write to temp file
  tmp_file <- tempfile(fileext = ".json")
  on.exit(unlink(tmp_file), add = TRUE)
  
  write_json_track(obj1, tmp_file)
  expect_true(file.exists(tmp_file))
  
  # Read back
  obj2 <- read_json_track(tmp_file)
  
  expect_s3_class(obj2, "JsonTrackObj")
  expect_equal(obj2$format, obj1$format)
  expect_equal(obj2$function_name, obj1$function_name)
  expect_equal(length(obj2$slices), 1)
})

test_that("as.data.frame.JsonTrackObj works correctly", {
  skip_if_not_installed("superassp")
  
  obj <- create_json_track_obj(
    results = list(measure1 = 100, measure2 = 200),
    function_name = "lst_test",
    file_path = "test.wav",
    sample_rate = 16000,
    audio_duration = 1.0,
    beginTime = 0.0,
    endTime = 1.0
  )
  
  # Add another slice
  obj <- append_json_track_slice(
    obj,
    results = list(measure1 = 110, measure2 = 210),
    beginTime = 1.0,
    endTime = 2.0
  )
  
  df <- as.data.frame(obj)
  
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 2)
  expect_equal(ncol(df), 4)  # begin_time, end_time, measure1, measure2
  expect_true("begin_time" %in% names(df))
  expect_true("end_time" %in% names(df))
  expect_equal(df$measure1, c(100, 110))
  expect_equal(df$measure2, c(200, 210))
  expect_type(df$measure1, "double")
  expect_type(df$measure2, "double")
})

test_that("as_tibble.JsonTrackObj works with tibble package", {
  skip_if_not_installed("superassp")
  skip_if_not_installed("tibble")
  
  obj <- create_json_track_obj(
    results = list(x = 1, y = 2),
    function_name = "lst_test",
    file_path = "test.wav",
    sample_rate = 16000,
    audio_duration = 1.0
  )
  
  tbl <- as_tibble(obj)
  
  expect_s3_class(tbl, "tbl_df")
  expect_equal(nrow(tbl), 1)
  expect_true("begin_time" %in% names(tbl))
})

test_that("read_track dispatches correctly", {
  skip_if_not_installed("superassp")
  skip_if_not_installed("jsonlite")
  
  # Create JSTF file
  obj <- create_json_track_obj(
    list(test = 123),
    "lst_test",
    "test.wav",
    16000,
    1.0
  )
  
  tmp_vat <- tempfile(fileext = ".vat")
  on.exit(unlink(tmp_vat), add = TRUE)
  
  write_json_track(obj, tmp_vat)
  
  # Read back using read_track
  result <- read_track(tmp_vat)
  
  expect_s3_class(result, "JsonTrackObj")
})

test_that("append_json_track_slice adds slices correctly", {
  skip_if_not_installed("superassp")
  
  obj <- create_json_track_obj(
    list(x = 1), "test", "test.wav", 16000, 3.0, 0, 1
  )
  
  expect_equal(length(obj$slices), 1)
  
  obj <- append_json_track_slice(obj, list(x = 2), 1.0, 2.0)
  expect_equal(length(obj$slices), 2)
  
  obj <- append_json_track_slice(obj, list(x = 3), 2.0, 3.0)
  expect_equal(length(obj$slices), 3)
  
  df <- as.data.frame(obj)
  expect_equal(nrow(df), 3)
  expect_equal(df$x, c(1, 2, 3))
})

test_that("merge_json_tracks combines multiple objects", {
  skip_if_not_installed("superassp")
  
  obj1 <- create_json_track_obj(
    list(x = 1), "test", "test.wav", 16000, 2.0, 0, 1
  )
  
  obj2 <- create_json_track_obj(
    list(x = 2), "test", "test.wav", 16000, 2.0, 1, 2
  )
  
  merged <- merge_json_tracks(obj1, obj2)
  
  expect_s3_class(merged, "JsonTrackObj")
  expect_equal(length(merged$slices), 2)
  
  df <- as.data.frame(merged)
  expect_equal(nrow(df), 2)
})

test_that("subset_json_track filters by time", {
  skip_if_not_installed("superassp")
  
  obj <- create_json_track_obj(
    list(x = 1), "test", "test.wav", 16000, 5.0, 0, 1
  )
  obj <- append_json_track_slice(obj, list(x = 2), 1, 2)
  obj <- append_json_track_slice(obj, list(x = 3), 2, 3)
  obj <- append_json_track_slice(obj, list(x = 4), 3, 4)
  
  # Filter to middle slices
  filtered <- subset_json_track(obj, start_time = 1.0, end_time = 3.0)
  
  expect_equal(length(filtered$slices), 2)
  df <- as.data.frame(filtered)
  expect_equal(df$x, c(2, 3))
})

test_that("get_jstf_extension returns correct extensions", {
  skip_if_not_installed("superassp")
  
  ext3 <- get_jstf_extension("lst_dysprosody")
  expect_equal(ext3, "dyp")
})

test_that("print.JsonTrackObj displays summary", {
  skip_if_not_installed("superassp")
  
  obj <- create_json_track_obj(
    list(x = 1, y = 2, z = 3),
    "lst_test",
    "test.wav",
    16000,
    5.0
  )
  
  # Should not error
  expect_output(print(obj), "JsonTrackObj")
  expect_output(print(obj), "lst_test")
  expect_output(print(obj), "JSTF")
})

test_that("summary.JsonTrackObj provides detailed info", {
  skip_if_not_installed("superassp")

  obj <- create_json_track_obj(
    list(measure1 = 100, measure2 = 200),
    "lst_test",
    "test.wav",
    16000,
    10.0
  )

  # Call method explicitly due to devtools::load_all() S3 dispatch issues
  expect_output(summary.JsonTrackObj(obj), "JSON Track Object Summary")
  expect_output(summary.JsonTrackObj(obj), "measure1")
  expect_output(summary.JsonTrackObj(obj), "measure2")
})
