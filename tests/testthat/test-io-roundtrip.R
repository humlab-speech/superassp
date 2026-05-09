# Round-trip tests for the symmetric I/O surface.
# Verify that read_<fmt> / write_<fmt> are true inverses for each track format,
# and that read_track() / write_track() dispatch correctly by extension/class.

test_that("write_jstf -> read_jstf round-trips a JsonTrackObj", {
  obj1 <- superassp:::create_json_track_obj(
    results = list(f0_mean = 150.0, f0_sd = 20.0),
    function_name = "lst_example",
    file_path = "audio.wav",
    sample_rate = 16000,
    audio_duration = 5.0
  )

  tmp <- tempfile(fileext = ".vat")
  on.exit(unlink(tmp), add = TRUE)
  write_jstf(obj1, tmp)
  expect_true(file.exists(tmp))

  obj2 <- read_jstf(tmp)
  expect_s3_class(obj2, "JsonTrackObj")
  expect_equal(obj2$function_name, "lst_example")
  expect_equal(obj2$sample_rate, 16000)
})

test_that("read_jstf accepts begin/end/samples for symmetry (currently no-op)", {
  obj <- superassp:::create_json_track_obj(
    results = list(f0_mean = 120.0),
    function_name = "lst_example",
    file_path = "audio.wav",
    sample_rate = 16000,
    audio_duration = 1.0
  )
  tmp <- tempfile(fileext = ".vat")
  on.exit(unlink(tmp), add = TRUE)
  write_jstf(obj, tmp)

  # Should not error when begin/end/samples are provided
  full     <- read_jstf(tmp)
  windowed <- read_jstf(tmp, begin = 0.5, end = 0.9, samples = FALSE)
  in_samp  <- read_jstf(tmp, begin = 100, end = 200, samples = TRUE)

  # All three return the same object today (args are reserved)
  expect_equal(full$function_name, windowed$function_name)
  expect_equal(full$function_name, in_samp$function_name)
})

test_that("read_track dispatches JSTF by extension", {
  obj <- superassp:::create_json_track_obj(
    results = list(f0_mean = 150.0),
    function_name = "lst_example",
    file_path = "audio.wav",
    sample_rate = 16000,
    audio_duration = 1.0
  )
  tmp <- tempfile(fileext = ".vat")
  on.exit(unlink(tmp), add = TRUE)
  write_jstf(obj, tmp)

  read_via_track <- read_track(tmp)
  read_via_jstf  <- read_jstf(tmp)
  expect_s3_class(read_via_track, "JsonTrackObj")
  expect_equal(read_via_track$function_name, read_via_jstf$function_name)
})

test_that("write_track dispatches JsonTrackObj to write_jstf", {
  obj <- superassp:::create_json_track_obj(
    results = list(f0_mean = 150.0),
    function_name = "lst_example",
    file_path = "audio.wav",
    sample_rate = 16000,
    audio_duration = 1.0
  )
  tmp <- tempfile(fileext = ".vat")
  on.exit(unlink(tmp), add = TRUE)
  write_track(obj, tmp)
  expect_true(file.exists(tmp))

  obj2 <- read_track(tmp)
  expect_s3_class(obj2, "JsonTrackObj")
})

test_that("write_track rejects unknown classes", {
  expect_error(write_track(list(foo = 1), tempfile()), "Cannot write object")
})
