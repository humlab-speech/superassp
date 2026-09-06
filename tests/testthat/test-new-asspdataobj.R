# Faithfulness gate for new_asspdataobj(): the shared constructor must reproduce
# the hand-built `list(...) + attr() + class()` pattern exactly — identical
# attributes AND identical written SSFF bytes — for every attribute signature
# used across the create_*_asspobj helpers.

# Build an AsspDataObj the old hand-rolled way, given which optional attrs to set.
build_inline <- function(tracks, trackFormats, sampleRate, endRecord,
                         origFreq = NULL, fileInfo = NULL, filePath = NULL,
                         startTime = 0.0, startRecord = 1L) {
  obj <- tracks
  attr(obj, "trackFormats") <- trackFormats
  attr(obj, "sampleRate")   <- sampleRate
  if (!is.null(origFreq)) attr(obj, "origFreq") <- origFreq
  attr(obj, "startTime")    <- startTime
  attr(obj, "startRecord")  <- as.integer(startRecord)
  attr(obj, "endRecord")    <- as.integer(endRecord)
  if (!is.null(fileInfo)) attr(obj, "fileInfo") <- fileInfo
  if (!is.null(filePath)) attr(obj, "filePath") <- filePath
  class(obj) <- "AsspDataObj"
  obj
}

# Same attributes (regardless of set order) => identical SSFF output.
expect_same_attrs <- function(a, b) {
  aa <- attributes(a); bb <- attributes(b)
  expect_setequal(names(aa), names(bb))
  for (nm in names(aa)) expect_equal(aa[[nm]], bb[[nm]], info = nm)
}

make_tracks <- function(ncol = 1L, n = 40L) {
  list(F0 = matrix(as.double(seq_len(n * ncol)), nrow = n, ncol = ncol))
}

test_that("full signature (f0-style: origFreq + fileInfo) matches inline", {
  tr <- make_tracks()
  old <- build_inline(tr, "REAL32", 100.0, 40L, origFreq = 44100, fileInfo = c(20L, 2L))
  new <- superassp:::new_asspdataobj(tr, sampleRate = 100.0, trackFormats = "REAL32",
                                     endRecord = 40L, origFreq = 44100,
                                     fileInfo = c(20L, 2L))
  expect_same_attrs(old, new)
})

test_that("origFreq without fileInfo matches inline", {
  tr <- make_tracks()
  old <- build_inline(tr, "REAL32", 100.0, 40L, origFreq = 16000)
  new <- superassp:::new_asspdataobj(tr, sampleRate = 100.0, trackFormats = "REAL32",
                                     endRecord = 40L, origFreq = 16000)
  expect_same_attrs(old, new)
})

test_that("filePath signature matches inline", {
  tr <- make_tracks()
  old <- build_inline(tr, "REAL32", 100.0, 40L, origFreq = 16000, filePath = "/tmp/x.wav")
  new <- superassp:::new_asspdataobj(tr, sampleRate = 100.0, trackFormats = "REAL32",
                                     endRecord = 40L, origFreq = 16000,
                                     filePath = "/tmp/x.wav")
  expect_same_attrs(old, new)
})

test_that("minimal signature (no origFreq/fileInfo/filePath) matches inline", {
  tr <- make_tracks()
  old <- build_inline(tr, "REAL32", 100.0, 40L)
  new <- superassp:::new_asspdataobj(tr, sampleRate = 100.0, trackFormats = "REAL32",
                                     endRecord = 40L)
  expect_same_attrs(old, new)
})

test_that("endRecord defaults to first-track row count", {
  tr <- make_tracks(n = 37L)
  new <- superassp:::new_asspdataobj(tr, sampleRate = 100.0, trackFormats = "REAL32")
  expect_equal(attr(new, "endRecord"), 37L)
})

test_that("written SSFF bytes are identical to inline construction", {
  skip_if(is.null(getNamespace("superassp")$write.AsspDataObj), "writer unavailable")
  tr <- make_tracks()
  old <- build_inline(tr, "REAL32", 100.0, 40L, origFreq = 44100, fileInfo = c(20L, 2L))
  new <- superassp:::new_asspdataobj(tr, sampleRate = 100.0, trackFormats = "REAL32",
                                     endRecord = 40L, origFreq = 44100,
                                     fileInfo = c(20L, 2L))
  f_old <- tempfile(fileext = ".f0"); f_new <- tempfile(fileext = ".f0")
  on.exit(unlink(c(f_old, f_new)))
  superassp:::write.AsspDataObj(old, f_old)
  superassp:::write.AsspDataObj(new, f_new)
  expect_identical(readBin(f_old, "raw", n = file.info(f_old)$size),
                   readBin(f_new, "raw", n = file.info(f_new)$size))
})
