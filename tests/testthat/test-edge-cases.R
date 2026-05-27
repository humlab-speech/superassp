# Edge case tests for trk_* functions and I/O helpers.

wav_file <- function() {
  f <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  if (f == "") skip("Test wav not found")
  f
}

# ---- File existence checks ----

test_that("trk_pitch_rapt errors on nonexistent file", {
  expect_error(
    trk_pitch_rapt("/no/such/file.wav", toFile = FALSE, verbose = FALSE),
    regexp = "do not exist"
  )
})

test_that("trk_d4c errors on nonexistent file", {
  expect_error(
    trk_d4c("/no/such/file.wav", toFile = FALSE, verbose = FALSE),
    regexp = "do not exist"
  )
})

# NOTE: trk_acf (C-ASSP path) does not have an upfront existence check — nonexistent
# files propagate as mclapply worker warnings. Not tested here to avoid hanging.

# ---- Time window validation ----

test_that("trk_acf errors when beginTime > endTime", {
  expect_error(
    trk_acf(wav_file(), beginTime = 2.0, endTime = 0.5, toFile = FALSE, verbose = FALSE),
    regexp = "time window|begin|end|start_time"
  )
})

test_that("trk_rms errors when beginTime > endTime", {
  expect_error(
    trk_rms(wav_file(), beginTime = 2.0, endTime = 0.5, toFile = FALSE, verbose = FALSE),
    regexp = "time window|begin|end|start_time"
  )
})

# ---- outputDirectory handling ----

test_that("trk_acf creates nonexistent outputDirectory", {
  tmpdir <- tempfile("superassp_test_")
  on.exit(unlink(tmpdir, recursive = TRUE))
  expect_false(dir.exists(tmpdir))
  # makeOutputDirectory is called even when toFile = FALSE
  trk_acf(wav_file(), toFile = FALSE, outputDirectory = tmpdir, verbose = FALSE)
  expect_true(dir.exists(tmpdir))
})

test_that("trk_acf errors when outputDirectory path is a file not a directory", {
  tmpfile <- tempfile()
  writeLines("x", tmpfile)
  on.exit(unlink(tmpfile))
  expect_error(
    trk_acf(wav_file(), toFile = FALSE, outputDirectory = tmpfile, verbose = FALSE),
    regexp = "not a directory"
  )
})

# ---- Malformed SSFF file ----

test_that("read_ssff errors on malformed file", {
  tmpf <- tempfile(fileext = ".f0")
  writeLines("this is not SSFF content", tmpf)
  on.exit(unlink(tmpf))
  expect_error(
    read_ssff(tmpf),
    regexp = "Unknown file format|not.*SSFF|invalid"
  )
})

test_that("read_track errors on malformed SSFF file", {
  tmpf <- tempfile(fileext = ".f0")
  writeLines("this is not SSFF content", tmpf)
  on.exit(unlink(tmpf))
  expect_error(
    read_track(tmpf),
    regexp = "Unknown file format|not.*SSFF|invalid"
  )
})

# ---- Batch processing ----

test_that("trk_pitch_rapt errors cleanly on empty listOfFiles", {
  expect_error(
    trk_pitch_rapt(character(0), toFile = FALSE, verbose = FALSE),
    regexp = "No input files|length|do not exist"
  )
})

test_that("trk_d4c errors cleanly on empty listOfFiles", {
  expect_error(
    trk_d4c(character(0), toFile = FALSE, verbose = FALSE),
    regexp = "No input files|length"
  )
})
