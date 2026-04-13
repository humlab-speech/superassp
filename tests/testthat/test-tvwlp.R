# Low-level C++ tests

test_that("tvlp_l2_cpp returns correct shape", {
  x <- as.numeric(cos(seq_len(320) * pi / 50))
  aki <- tvlp_l2_cpp(x, p = 8L, q = 3L)
  expect_equal(dim(aki), c(4L, 8L))
})

test_that("tvwlp_l2_cpp returns correct shape", {
  x <- rnorm(320)
  w <- seq(0, 1, length.out = 320)
  aki <- tvwlp_l2_cpp(x, p = 8L, q = 3L, w = w)
  expect_equal(dim(aki), c(4L, 8L))
})

test_that("tvlptoformants_akitofi_cpp returns fi and bw", {
  x <- as.numeric(cos(seq_len(320) * pi / 50))
  aki <- tvlp_l2_cpp(x, p = 8L, q = 3L)
  out <- tvlptoformants_akitofi_cpp(aki, nx = 100L, npeaks = 3L, fs = 8000)
  expect_named(out, c("fi", "bw", "ak"))
  expect_equal(dim(out$fi), c(100L, 3L))
  expect_equal(dim(out$bw), c(100L, 3L))
  expect_true(all(is.finite(out$fi)))
})

test_that("pre_emphasis_cpp works", {
  x <- rnorm(100)
  pe <- pre_emphasis_cpp(x, 0.97)
  expect_equal(nrow(pe), 100L)
  expect_equal(ncol(pe), 1L)
  expect_equal(pe[1, 1], x[1])
})

test_that("get_lpc_residual_cpp preserves length", {
  x <- rnorm(512)
  r <- get_lpc_residual_cpp(x, l = 200L, shift = 40L, order = 12L)
  expect_length(r, 512L)
  expect_true(all(is.finite(r)))
})

# End-to-end wrapper tests

test_that("trk_formant_tvwlp() works with single WAV", {
  wfile <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(wfile == "", "Test audio not found")
  res <- trk_formant_tvwlp(wfile, toFile = FALSE, verbose = FALSE)
  expect_s3_class(res, "AsspDataObj")
  expect_true("fm" %in% names(res))
  expect_true("bw" %in% names(res))
  expect_equal(ncol(res$fm), 3L)
  expect_true(nrow(res$fm) > 0)
  expect_true(all(is.finite(res$fm)))
})

test_that("trk_formant_tvwlp() file output works", {
  wfile <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(wfile == "", "Test audio not found")
  n <- trk_formant_tvwlp(wfile, toFile = TRUE,
                          outputDirectory = tempdir(), verbose = FALSE)
  expect_equal(n, 1L)
  expect_true(file.exists(file.path(tempdir(), "a1.tvf")))
})

test_that("trk_formant_tvwlp() tvlp_l2 method works", {
  wfile <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(wfile == "", "Test audio not found")
  res <- trk_formant_tvwlp(wfile, lptype = "tvlp_l2",
                            toFile = FALSE, verbose = FALSE)
  expect_s3_class(res, "AsspDataObj")
  expect_equal(ncol(res$fm), 3L)
})

test_that("trk_formant_tvwlp() repeated calls don't crash", {
  wfile <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(wfile == "", "Test audio not found")
  for (i in 1:3) {
    res <- trk_formant_tvwlp(wfile, toFile = FALSE, verbose = FALSE)
  }
  expect_s3_class(res, "AsspDataObj")
})

test_that("trk_formant_tvwlp() function attributes are correct", {
  expect_equal(attr(trk_formant_tvwlp, "ext"), "tvf")
  expect_equal(attr(trk_formant_tvwlp, "tracks"), c("fm", "bw"))
  expect_equal(attr(trk_formant_tvwlp, "outputType"), "SSFF")
})
