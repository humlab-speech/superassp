# Golden-file parity with wrassp.
#
# The fixtures in golden/ssff/ were produced by wrassp -- the reference
# implementation of the same SSFF reader/writer that superassp carries -- from
# the bundled sample recording; golden/ssff/PROVENANCE.md records how. Next to
# every fixture sits the object that wrassp::read.AsspDataObj() returned for it
# at generation time (.expected.rds, xz-compressed).
#
# These tests therefore pin superassp's reader and writer to an *external*
# implementation of SSFF, and they run everywhere: wrassp is not needed, only
# the committed fixtures are. The live cross-package checks (regenerating the
# files with the installed wrassp) live in test-ssff-wrassp-interop.R.

golden_dir <- testthat::test_path("golden", "ssff")

golden_cases <- function() {
  files <- list.files(golden_dir, pattern = "\\.ssff$", full.names = TRUE)
  stats::setNames(files, sub("\\.ssff$", "", basename(files)))
}

golden_expected <- function(name) {
  readRDS(file.path(golden_dir, paste0(name, ".expected.rds")))
}

skip_without_golden <- function() {
  if (!length(golden_cases())) testthat::skip("SSFF gold fixtures are unavailable")
}

# Compare an object read by superassp against the wrassp-produced expectation.
# filePath is deliberately not compared (it names the reader's own file).
expect_golden_match <- function(actual, expected, label = "") {
  expect_identical(names(actual), names(expected), info = label)
  for (track in names(expected)) {
    expect_identical(dim(actual[[track]]), dim(expected[[track]]),
                     info = paste(label, track, "dims"))
    expect_identical(typeof(actual[[track]]), typeof(expected[[track]]),
                     info = paste(label, track, "type"))
    expect_equal(actual[[track]], expected[[track]], tolerance = 0,
                 info = paste(label, track, "values"))
  }
  for (field in c("sampleRate", "startTime", "startRecord", "endRecord",
                  "trackFormats", "origFreq")) {
    expect_equal(attr(actual, field), attr(expected, field), tolerance = 0,
                 info = paste(label, field))
  }
  invisible(TRUE)
}

raw_bytes <- function(path) readBin(path, "raw", file.size(path))

test_that("the gold fixture set is complete and matches its provenance record", {
  skip_without_golden()
  provenance <- readLines(file.path(golden_dir, "PROVENANCE.md"), warn = FALSE)
  rows <- grep("^\\| `.*\\.ssff` \\|", provenance, value = TRUE)
  expect_length(rows, length(golden_cases()))

  for (row in rows) {
    fields <- trimws(strsplit(row, "|", fixed = TRUE)[[1]])
    name <- gsub("`", "", fields[[2]])
    size <- as.integer(fields[[3]])
    md5  <- gsub("`", "", fields[[4]])
    path <- file.path(golden_dir, name)
    expect_true(file.exists(path), info = name)
    expect_equal(file.size(path), size, info = name)
    expect_equal(unname(tools::md5sum(path)), md5, info = name)
  }
})

test_that("every gold fixture reads exactly as wrassp read it", {
  skip_without_golden()
  for (name in names(golden_cases())) {
    expect_golden_match(read_ssff(golden_cases()[[name]]), golden_expected(name), label = name)
  }
})

test_that("the gold fixtures cover the storage classes and shapes they claim", {
  skip_without_golden()
  shapes <- list(
    f0_ksv       = list(tracks = "F0",    formats = "REAL32", type = "double",  dim = c(805L, 1L)),
    f0_mhs       = list(tracks = "pitch", formats = "REAL32", type = "double",  dim = c(805L, 1L)),
    rms          = list(tracks = "rms",   formats = "REAL32", type = "double",  dim = c(805L, 1L)),
    spectrum_dft = list(tracks = "dft",   formats = "REAL32", type = "double",  dim = c(15L, 1025L)),
    acf48        = list(tracks = "acf",   formats = "REAL64", type = "double",  dim = c(60L, 48L)),
    forest       = list(tracks = c("fm", "bw"), formats = c("INT16", "INT16"),
                        type = "integer", dim = c(60L, 4L)),
    f0_ksv_be    = list(tracks = "F0",    formats = "REAL32", type = "double",  dim = c(805L, 1L))
  )
  for (name in names(shapes)) {
    sh <- shapes[[name]]
    obj <- read_ssff(golden_cases()[[name]])
    expect_identical(names(obj), sh$tracks, info = name)
    expect_identical(attr(obj, "trackFormats"), sh$formats, info = name)
    for (track in names(obj)) {
      expect_identical(dim(obj[[track]]), sh$dim, info = paste(name, track))
      expect_identical(typeof(obj[[track]]), sh$type, info = paste(name, track))
    }
  }
})

test_that("little-endian gold fixtures survive a read -> write round trip byte-identically", {
  skip_without_golden()
  # f0_ksv_be is re-encoded on write (we always write host byte order) and
  # na_written_by_wrassp holds NaN, which our writer normalises to 0; both are
  # covered by their own assertions below.
  cases <- setdiff(names(golden_cases()), c("f0_ksv_be", "na_written_by_wrassp"))
  skip_if(!length(cases), "no round-trippable fixtures")
  for (name in cases) {
    original <- golden_cases()[[name]]
    back <- tempfile(fileext = ".ssff")
    on.exit(unlink(back), add = TRUE)
    write_ssff(read_ssff(original), back)
    expect_identical(raw_bytes(back), raw_bytes(original), info = name)
  }
})

test_that("the big-endian fixture decodes exactly like its little-endian twin", {
  skip_without_golden()
  be <- read_ssff(golden_cases()[["f0_ksv_be"]])
  le <- read_ssff(golden_cases()[["f0_ksv"]])
  expect_equal(be$F0, le$F0, tolerance = 0)
  expect_golden_match(be, golden_expected("f0_ksv_be"), label = "f0_ksv_be")
  # re-writing the swapped file yields the little-endian encoding
  back <- tempfile(fileext = ".ssff")
  on.exit(unlink(back), add = TRUE)
  write_ssff(be, back)
  expect_identical(raw_bytes(back), raw_bytes(golden_cases()[["f0_ksv"]]))
})

test_that("read_track() masks exactly the stored zeros, read_ssff() keeps them", {
  skip_without_golden()
  path <- golden_cases()[["f0_ksv"]]
  expected <- golden_expected("f0_ksv")$F0

  expect_equal(sum(expected == 0), 272)          # fixture property (unvoiced frames)

  raw <- read_ssff(path)$F0
  expect_equal(raw, expected, tolerance = 0)
  expect_equal(sum(raw == 0), 272)

  masked <- read_track(path)$F0
  expect_equal(sum(is.na(masked)), 272)
  expect_identical(is.na(masked), expected == 0)
  expect_equal(masked[!is.na(masked)], expected[expected != 0], tolerance = 0)
  expect_equal(read_track(path, zero_to_na = FALSE)$F0, expected, tolerance = 0)

  # the windowed view masks the same way
  win <- read_track(path, begin = 0.5, end = 1.5)$F0
  expect_equal(sum(is.na(win)), sum(expected[101:300] == 0))
})

test_that("wrassp's NA/NaN encoding passes through our reader and is normalised by our writer", {
  skip_without_golden()
  path <- golden_cases()[["na_written_by_wrassp"]]
  expected <- golden_expected("na_written_by_wrassp")$f0
  expect_true(is.nan(expected[1]) && is.nan(expected[3]))

  raw <- read_ssff(path)$f0
  expect_true(is.nan(raw[1]) && is.nan(raw[3]))
  expect_equal(raw[c(2, 4)], c(1.5, -2), tolerance = 0)

  # read_track() does not turn NaN into NA: only a stored 0 means "missing"
  masked <- read_track(path)$f0
  expect_true(is.nan(masked[1]) && !is.na(masked[2]))

  # our writer stores NA/NaN as 0 (the documented difference from wrassp)
  back <- tempfile(fileext = ".ssff")
  on.exit(unlink(back), add = TRUE)
  write_ssff(read_ssff(path), back)
  expect_equal(as.numeric(read_ssff(back)$f0), c(0, 1.5, 0, -2), tolerance = 0)
})

test_that("tracks = and threads = behave on the gold files", {
  skip_without_golden()
  forest <- golden_cases()[["forest"]]
  full   <- read_ssff(forest)

  only <- read_ssff(forest, tracks = "bw")
  expect_identical(names(only), "bw")
  expect_equal(only$bw, full$bw, tolerance = 0)
  expect_identical(attr(only, "trackFormats"), "INT16")
  expect_error(read_ssff(forest, tracks = "nope"), "not found")
  expect_error(read_ssff(forest, tracks = "nope"), "fm, bw")

  wide <- golden_cases()[["spectrum_dft"]]
  expect_identical(read_ssff(wide, threads = 4L)$dft, read_ssff(wide, threads = 1L)$dft)
})

test_that("delTrack() and addTrack() keep the container consistent through write/read", {
  skip_without_golden()
  forest <- read_ssff(golden_cases()[["forest"]])
  out <- tempfile(fileext = ".ssff")
  on.exit(unlink(out), add = TRUE)

  reduced <- delTrack(forest, "bw")
  expect_identical(names(reduced), "fm")
  expect_identical(attr(reduced, "trackFormats"), "INT16")
  write_ssff(reduced, out)
  expect_equal(read_ssff(out)$fm, forest$fm, tolerance = 0)

  replaced <- addTrack(forest, "bw", forest$bw, format = "INT16", deleteExisting = TRUE)
  expect_identical(attr(replaced, "trackFormats"), c("INT16", "INT16"))
  write_ssff(replaced, out)
  expect_equal(read_ssff(out)$bw, forest$bw, tolerance = 0)

  # adding a *new* track registers its format (it used to be dropped, which made
  # the object unwritable)
  extra <- matrix(seq_len(nrow(forest$fm)), ncol = 1)
  extended <- addTrack(forest, "extra", extra, format = "REAL32")
  expect_identical(names(extended), c("fm", "bw", "extra"))
  expect_identical(attr(extended, "trackFormats"), c("INT16", "INT16", "REAL32"))
  expect_length(track_formats(extended), length(names(extended)))
  write_ssff(extended, out)
  back <- read_ssff(out)
  expect_identical(names(back), c("fm", "bw", "extra"))
  expect_identical(attr(back, "trackFormats"), c("INT16", "INT16", "REAL32"))
  expect_equal(back$extra[, 1], seq_len(nrow(forest$fm)), tolerance = 0)
})
