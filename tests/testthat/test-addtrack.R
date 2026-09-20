# Container mechanics of addTrack()/delTrack().
#
# addTrack() used to discard the result of append() when it added a *new* track,
# so the track appeared in the data and in the tibble/data.frame views but not
# in trackFormats()/attr(., "trackFormats") -- and the object could not be
# written ("Not enough format specifiers for the data tracks."). Every caller in
# the package worked around it; the wrappers that build their result with
# addTrack() therefore assert the invariant below as well.

test_that("addTrack() registers the format of a new track", {
  obj <- list(a = matrix(1:6, ncol = 2))
  attr(obj, "trackFormats") <- "REAL32"
  attr(obj, "sampleRate")   <- 100
  attr(obj, "origFreq")     <- 44100
  attr(obj, "startTime")    <- 0
  attr(obj, "startRecord")  <- 1L
  attr(obj, "endRecord")    <- 3L
  attr(obj, "fileInfo")     <- c(20L, 2L)
  class(obj) <- "AsspDataObj"

  extended <- addTrack(obj, "b", matrix(1:3, ncol = 1), "INT16")
  expect_identical(names(extended), c("a", "b"))
  expect_identical(attr(extended, "trackFormats"), c("REAL32", "INT16"))
  expect_length(track_formats(extended), length(names(extended)))

  out <- tempfile(fileext = ".ssff")
  on.exit(unlink(out), add = TRUE)
  write_ssff(extended, out)
  back <- read_ssff(out)
  expect_identical(names(back), c("a", "b"))
  expect_identical(attr(back, "trackFormats"), c("REAL32", "INT16"))
  expect_equal(back$b[, 1], 1:3, tolerance = 0)

  # a second new track appends again, keeping the order of the calls
  twice <- addTrack(extended, "c", matrix(4:6, ncol = 1), "REAL64")
  expect_identical(attr(twice, "trackFormats"), c("REAL32", "INT16", "REAL64"))
})

test_that("addTrack() replaces a track in place and knows its edge cases", {
  obj <- list(a = matrix(1:6, ncol = 2))
  attr(obj, "trackFormats") <- "REAL32"
  class(obj) <- "AsspDataObj"

  # replacing keeps the number of formats (the attribute is written in place)
  replaced <- addTrack(obj, "a", matrix(7:12, ncol = 2), "REAL64", deleteExisting = TRUE)
  expect_identical(attr(replaced, "trackFormats"), "REAL64")
  expect_length(track_formats(replaced), length(names(replaced)))

  # replacing the *only* track relaxes the row check (single-track branch); it
  # still needs deleteExisting, like any other replacement
  only <- addTrack(obj, "a", matrix(1:2, ncol = 1), "INT16", deleteExisting = TRUE)
  expect_identical(attr(only, "trackFormats"), "INT16")
  expect_identical(dim(only$a), c(2L, 1L))

  # guarded edge cases
  expect_error(addTrack(obj, "a", matrix(1:2, ncol = 1), "REAL32"), "will not be deleted")
  expect_error(addTrack(obj, "b", matrix(1:2, ncol = 1), "REAL32"), "number of rows")
  expect_error(addTrack(obj, "b", "not numeric", "REAL32"), "numeric matrix")
  expect_error(addTrack(obj, c("b", "c"), matrix(1:3, ncol = 1), "REAL32"), "atomic string")
  expect_error(addTrack(list(a = 1), "b", matrix(1:2, ncol = 1), "REAL32"), "AsspDataObj")
})

test_that("delTrack() keeps trackFormats in step with the tracks", {
  obj <- list(a = matrix(1:6, ncol = 2), b = matrix(1:3, ncol = 1))
  attr(obj, "trackFormats") <- c("REAL32", "INT16")
  class(obj) <- "AsspDataObj"

  reduced <- delTrack(obj, "a")
  expect_identical(names(reduced), "b")
  expect_identical(attr(reduced, "trackFormats"), "INT16")
  expect_error(delTrack(obj, "nope"), "Invalid trackname")
})

test_that("harmonics() builds a consistent container", {
  # internal helper (R/ssff.R) that adds its track with addTrack()
  f0 <- list(f0 = matrix(seq(100, 200, length.out = 40), ncol = 1))
  attr(f0, "trackFormats") <- "REAL32"
  attr(f0, "sampleRate")   <- 100
  attr(f0, "origFreq")     <- 44100
  attr(f0, "startTime")    <- 0
  attr(f0, "startRecord")  <- 1L
  attr(f0, "endRecord")    <- 40L
  attr(f0, "fileInfo")     <- c(20L, 2L)
  class(f0) <- "AsspDataObj"

  path <- tempfile(fileext = ".f0")
  on.exit(unlink(path), add = TRUE)
  write_ssff(f0, path)

  har <- harmonics(path, column = "f0", n = 3, toFile = FALSE)
  expect_s3_class(har, "AsspDataObj")
  expect_identical(names(har), "har")
  expect_identical(attr(har, "trackFormats"), "INT16")
  expect_length(track_formats(har), length(names(har)))

  out <- tempfile(fileext = ".har")
  on.exit(unlink(out), add = TRUE)
  write_ssff(har, out)
  expect_identical(names(read_ssff(out)), "har")
})
