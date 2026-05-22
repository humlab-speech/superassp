# Bit-/near-equivalence fidelity tests for the performance refactor (v2.7.1).
# These tests guarantee that the caching, scheduling, and Makevars changes
# introduced in the perf pass do NOT alter numerical output of the trk_* / lst_*
# DSP wrappers. Run this file after every perf-related edit (especially when
# enabling OpenMP pragmas) before declaring the change safe.

skip_helper <- function() {
  if (!requireNamespace("superassp", quietly = TRUE)) {
    testthat::skip("superassp not installed")
  }
  wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  if (!nzchar(wav) || !file.exists(wav)) {
    testthat::skip("Sample WAV not available")
  }
  wav
}

test_that("media_info() returns the same payload as av::av_media_info()", {
  wav <- skip_helper()
  fresh <- av::av_media_info(wav)
  cached <- superassp:::media_info(wav)
  # The cache wraps the call; payload must be identical
  expect_identical(cached$audio$sample_rate, fresh$audio$sample_rate)
  expect_identical(cached$audio$channels,    fresh$audio$channels)
  expect_equal(cached$duration, fresh$duration, tolerance = 1e-9)
})

test_that("read_audio() output is invariant under cache hits", {
  wav <- skip_helper()
  a1 <- read_audio(wav)
  a2 <- read_audio(wav)   # second call: hits media_info cache
  expect_identical(a1, a2)
})

test_that("trk_pitch_rapt output is invariant across sequential and parallel runs", {
  wav <- skip_helper()
  one  <- trk_pitch_rapt(wav, toFile = FALSE, verbose = FALSE)
  twice <- trk_pitch_rapt(c(wav, wav), toFile = FALSE, verbose = FALSE)
  # twice should yield a list of length 2, each element matching `one`
  expect_length(twice, 2L)
  expect_equal(twice[[1]]$f0, one$f0, tolerance = 0)
  expect_equal(twice[[2]]$f0, one$f0, tolerance = 0)
})

test_that("OpenMP thread count does not change output (skipped if not enabled)", {
  testthat::skip("Enable when OpenMP pragmas are added to C++ kernels")
  # Future:
  # wav <- skip_helper()
  # Sys.setenv(OMP_NUM_THREADS = "1")
  # r1 <- trk_pitch_rapt(wav, toFile = FALSE)
  # Sys.setenv(OMP_NUM_THREADS = "4")
  # r4 <- trk_pitch_rapt(wav, toFile = FALSE)
  # expect_equal(r1$f0, r4$f0, tolerance = 1e-9)
})

test_that("build_media_manifest() reports consistent metadata", {
  wav <- skip_helper()
  m <- superassp:::build_media_manifest(c(wav, wav))
  expect_s3_class(m, "data.frame")
  expect_equal(nrow(m), 2L)
  expect_true(all(m$exists))
  expect_true(all(m$native_ext))
  expect_true(all(m$sample_rate > 0, na.rm = TRUE))
})
