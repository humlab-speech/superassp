f0_path <- testthat::test_path("golden", "ssff", "f0_ksv.ssff")

test_that("default (snap = 'none') still errors on off-grid begin == end", {
  skip_if(!file.exists(f0_path), "fixture not found")
  expect_error(read_ssff(f0_path, begin = 0.01, end = 0.01))
})

test_that("snap = 'nearest' resolves an off-grid single time point to exactly 1 record", {
  skip_if(!file.exists(f0_path), "fixture not found")

  snapped <- read_ssff(f0_path, begin = 0.01, end = 0.01, snap = "nearest")
  expect_equal(n_records(snapped), 1L)

  rate  <- attr(snapped, "sampleRate")
  idx   <- round(0.01 * rate)
  exact <- read_ssff(f0_path, begin = idx / rate, end = idx / rate)
  expect_equal(snapped[[1]], exact[[1]])
})

test_that("snap = 'nearest' is a no-op on an exact frame boundary", {
  skip_if(!file.exists(f0_path), "fixture not found")

  full <- read_ssff(f0_path)
  rate <- attr(full, "sampleRate")
  t    <- 5 / rate

  plain   <- read_ssff(f0_path, begin = t, end = t)
  snapped <- read_ssff(f0_path, begin = t, end = t, snap = "nearest")
  expect_equal(plain[[1]], snapped[[1]])
})

test_that("begin == end == 0 still means whole file regardless of snap", {
  skip_if(!file.exists(f0_path), "fixture not found")

  full      <- read_ssff(f0_path)
  full_snap <- read_ssff(f0_path, snap = "nearest")
  expect_equal(n_records(full), n_records(full_snap))
})

test_that("read_track propagates snap for SSFF but ignores it for JSTF", {
  skip_if(!file.exists(f0_path), "fixture not found")

  expect_error(read_track(f0_path, begin = 0.01, end = 0.01))
  expect_no_error(read_track(f0_path, begin = 0.01, end = 0.01, snap = "nearest"))
})

test_that("default snap='none' error message mentions snap parameter", {
  skip_if(!file.exists(f0_path), "fixture not found")

  err <- tryCatch(read_ssff(f0_path, begin = 0.01, end = 0.01), error = function(e) e)
  expect_match(conditionMessage(err), "snap", ignore.case = TRUE)
})
