# Live cross-package parity: superassp's SSFF reader/writer against the
# installed wrassp, on files wrassp generates from the bundled sample recording.
#
# The wrassp side runs in a subprocess (helper-wrassp-interop.R): loading wrassp
# into this session would register its print()/as_tibble() methods for
# AsspDataObj on top of superassp's and make dispatch order-dependent. Skips
# when wrassp is not installed; fails when it is installed but the interop run
# breaks.

interop <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) cached <<- wrassp_interop()
    cached
  }
})

skip_without_wrassp <- function() {
  data <- interop()
  if (is.null(data)) testthat::skip("wrassp is not installed")
  data
}

test_that("freshly generated wrassp analyses read identically through read_ssff()", {
  data <- skip_without_wrassp()
  # the interop run must not leave wrassp's namespace loaded: it would register
  # wrassp's as_tibble()/print() methods over superassp's for the rest of the
  # session (see helper-wrassp-interop.R)
  expect_false("wrassp" %in% loadedNamespaces())
  expect_length(data$expected, length(data$files))
  for (name in names(data$files)) {
    actual <- read_ssff(data$files[[name]])
    expected <- data$expected[[name]]
    expect_identical(names(actual), names(expected), info = name)
    for (track in names(expected)) {
      expect_identical(dim(actual[[track]]), dim(expected[[track]]), info = paste(name, track))
      expect_identical(typeof(actual[[track]]), typeof(expected[[track]]), info = paste(name, track))
      expect_equal(actual[[track]], expected[[track]], tolerance = 0, info = paste(name, track))
    }
    for (field in c("sampleRate", "startTime", "startRecord", "endRecord", "trackFormats")) {
      expect_equal(attr(actual, field), attr(expected, field), tolerance = 0,
                   info = paste(name, field))
    }
  }
})

test_that("windowed reads agree with wrassp record for record", {
  data <- skip_without_wrassp()
  path <- data$files[["f0_ksv"]]
  for (window in data$windows) {
    actual   <- read_ssff(path, begin = window$begin, end = window$end)
    expected <- window$obj
    label <- paste("begin", window$begin, "end", window$end)
    expect_identical(dim(actual$F0), dim(expected$F0), info = label)
    expect_equal(actual$F0, expected$F0, tolerance = 0, info = label)
    expect_equal(attr(actual, "startRecord"), attr(expected, "startRecord"), info = label)
    expect_equal(attr(actual, "endRecord"),   attr(expected, "endRecord"),   info = label)
    expect_equal(attr(actual, "startTime"),   attr(expected, "startTime"),   info = label)
    expect_equal(attr(actual, "sampleRate"),  attr(expected, "sampleRate"),  info = label)
  }
})

test_that("the legacy read.AsspDataObj() alias stays verbatim like wrassp", {
  data <- skip_without_wrassp()
  # read.AsspDataObj()/getAsspDataObj() are internal aliases (not exported, see
  # the export policy test); the test environment resolves them from the
  # package namespace.
  for (name in names(data$files)) {
    alias    <- read.AsspDataObj(data$files[[name]])
    expected <- data$expected[[name]]
    expect_identical(names(alias), names(expected), info = name)
    for (track in names(expected)) {
      expect_equal(alias[[track]], expected[[track]], tolerance = 0, info = paste(name, track))
    }
    expect_identical(alias[[1]], getAsspDataObj(data$files[[name]])[[1]], info = name)
  }
})

test_that("write_ssff() emits the same bytes as wrassp::write.AsspDataObj()", {
  data <- skip_without_wrassp()
  for (name in names(data$writers)) {
    ours <- tempfile(fileext = ".ssff")
    on.exit(unlink(ours), add = TRUE)
    write_ssff(data$writers[[name]], ours)
    theirs <- data$writer_files[[name]]
    expect_identical(readBin(ours, "raw", file.size(ours)),
                     readBin(theirs, "raw", file.size(theirs)), info = name)
  }
})

test_that("delTrack() agrees with wrassp's", {
  data <- skip_without_wrassp()
  ours     <- delTrack(read_ssff(data$files[["forest"]]), "bw")
  expected <- data$deltrack
  expect_identical(names(ours), names(expected))
  expect_equal(ours[[1]], expected[[1]], tolerance = 0)
  expect_identical(attr(ours, "trackFormats"), attr(expected, "trackFormats"))
})
