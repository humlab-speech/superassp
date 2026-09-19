# Tests for the SSFF read engine: single-pass conversion, mmap window reads,
# endianness swapping, track selection, threading and the 0-is-missing rule.
#
# The expected values are derived independently of the reader: the data section
# is located in R and decoded with readBin(), so these tests fail if the C
# kernels and the file format ever disagree.

# 0-based file offset of the first data byte (just past the header marker line).
ssff_data_offset <- function(path) {
  raw <- readBin(path, "raw", min(65536L, file.size(path)))
  needle <- charToRaw("-----------------\n")
  pos <- grepRaw(needle, raw, fixed = TRUE)
  if (length(pos) == 0) stop("no SSFF header marker in ", path)
  pos - 1L + length(needle)
}

# Decode 'n' values of 'size' bytes starting at 0-based file offset 'from'.
ssff_ref_values <- function(path, n, from = 0L, size = 4L, float = FALSE,
                            endian = "little", signed = TRUE) {
  con <- file(path, "rb")
  on.exit(close(con))
  seek(con, from, origin = "start")
  what <- if (float || size == 8L) "numeric" else "integer"
  as.numeric(readBin(con, what, n = n, size = size, endian = endian, signed = signed))
}

write_test_obj <- function(path, mats, formats, rate = 100) {
  obj <- mats
  attr(obj, "trackFormats") <- formats
  attr(obj, "sampleRate")   <- rate
  attr(obj, "origFreq")     <- 44100
  attr(obj, "startTime")    <- 0
  attr(obj, "startRecord")  <- 1L
  attr(obj, "endRecord")    <- nrow(mats[[1]])
  attr(obj, "fileInfo")     <- c(20L, 2L)
  class(obj) <- "AsspDataObj"
  write_ssff(obj, path)
  invisible(path)
}

test_that("reader matches an independent readBin() decode for every storage type", {
  # UINT16/UINT32 have no SSFF storage class, so the writer refuses them
  cases <- list(
    REAL32 = list(type = "REAL32", size = 4L, float = TRUE, tol = 1e-4,
                  vals = function(n) round(rnorm(n) * 100, 2)),
    REAL64 = list(type = "REAL64", size = 8L, float = TRUE, tol = 1e-12,
                  vals = function(n) rnorm(n)),
    INT16  = list(type = "INT16",  size = 2L, float = FALSE, tol = 0, signed = TRUE,
                  vals = function(n) sample(-30000:30000, n)),
    INT32  = list(type = "INT32",  size = 4L, float = FALSE, tol = 0, signed = TRUE,
                  vals = function(n) sample(-100000:100000, n)),
    UINT8  = list(type = "UINT8",  size = 1L, float = FALSE, tol = 0, signed = FALSE,
                  vals = function(n) sample(0:255, n))
  )
  set.seed(11)
  for (nm in names(cases)) {
    cs   <- cases[[nm]]
    nrec <- 40L
    v    <- cs$vals(nrec)
    path <- tempfile(fileext = ".ssff")
    write_test_obj(path, list(tr = matrix(as.numeric(v), ncol = 1)), cs$type)
    on.exit(unlink(path), add = TRUE)

    obj <- read_ssff(path)
    expect_named(obj, "tr")
    expect_equal(dim(obj$tr), c(nrec, 1L), info = nm)

    # must match the bytes exactly
    ref <- ssff_ref_values(path, nrec, from = ssff_data_offset(path),
                           size = cs$size, float = cs$float,
                           signed = is.null(cs$signed) || cs$signed)
    expect_equal(as.numeric(obj$tr), ref, tolerance = 0, info = paste(nm, "vs readBin"))

    # and the source values, within the precision of the storage type
    if (cs$tol == 0) {
      expect_equal(as.numeric(obj$tr), as.numeric(v), info = nm)
    } else {
      expect_equal(as.numeric(obj$tr), as.numeric(v), tolerance = cs$tol, info = nm)
    }
  }
})

test_that("multi-track and multi-field files convert with the right layout", {
  set.seed(12)
  nrec <- 25L
  a <- matrix(round(rnorm(nrec * 3) * 10, 3), nrec, 3)   # 3 fields
  b <- matrix(round(rnorm(nrec) * 1000, 3), nrec, 1)     # 1 field
  path <- tempfile(fileext = ".ssff")
  write_test_obj(path, list(alpha = a, beta = b), c("REAL32", "REAL32"))
  on.exit(unlink(path), add = TRUE)

  obj <- read_ssff(path)
  expect_named(obj, c("alpha", "beta"))
  expect_equal(unname(obj$alpha), a, tolerance = 1e-3)
  expect_equal(as.numeric(obj$beta), as.numeric(b), tolerance = 1e-3)

  # raw bytes: record = [a1 a2 a3 b], 16 bytes per record
  off     <- ssff_data_offset(path)
  recsize <- 4L * 4L                      # [a1 a2 a3 b] per record
  stride  <- recsize %/% 4L
  # read the whole interleaved run starting at a field and keep every stride-th value
  ref_field <- function(field) {
    vals <- ssff_ref_values(path, nrec * stride, from = off + field * 4L,
                            size = 4L, float = TRUE)
    vals[seq(1L, length.out = nrec, by = stride)]
  }
  for (f in 0:2) {
    expect_equal(obj$alpha[, f + 1L], ref_field(f), tolerance = 0,
                 info = paste("alpha field", f))
  }
  expect_equal(as.numeric(obj$beta), ref_field(3L), tolerance = 0)
})

test_that("big-endian files convert to the same values as little-endian files", {
  set.seed(18)
  v <- round(rnorm(9), 3)
  le <- tempfile(fileext = ".ssff")
  be <- tempfile(fileext = ".ssff")
  on.exit(unlink(c(le, be)))
  write_test_obj(le, list(tr = matrix(v, ncol = 1)), "REAL32")

  # swap every 4-byte group of the data section and declare MSB-first
  raw <- readBin(le, "raw", file.size(le))
  off <- ssff_data_offset(le)
  body <- raw[(off + 1L):length(raw)]
  swapped <- as.vector(matrix(as.integer(body), nrow = 4)[4:1, ])
  hdr_txt <- rawToChar(raw[1:off])
  hdr_txt <- sub("Machine IBM-PC", "Machine SPARC  ", hdr_txt, fixed = TRUE)
  writeBin(c(charToRaw(hdr_txt), as.raw(swapped)), be)

  expect_equal(as.numeric(read_ssff(be)$tr), as.numeric(read_ssff(le)$tr))
  expect_equal(as.numeric(read_ssff(be)$tr),
               ssff_ref_values(be, length(v), from = ssff_data_offset(be),
                               size = 4L, float = TRUE, endian = "big"),
               tolerance = 0)
})

test_that("windowed reads return the requested records", {
  set.seed(13)
  v <- rnorm(100)
  path <- tempfile(fileext = ".ssff")
  write_test_obj(path, list(tr = matrix(v, ncol = 1)), "REAL32")
  on.exit(unlink(path), add = TRUE)

  obj <- read_ssff(path, begin = 0.2, end = 0.4)   # 100 Hz -> records 20..40
  expect_equal(n_records(obj), 21L)
  expect_equal(attr(obj, "startRecord"), 21L)
  expect_equal(attr(obj, "endRecord"), 41L)
  expect_equal(as.numeric(obj$tr), v[21:41], tolerance = 1e-6)

  bySample <- read_ssff(path, begin = 10, end = 19, samples = TRUE)
  expect_equal(n_records(bySample), 10L)
  expect_equal(as.numeric(bySample$tr), v[11:20], tolerance = 1e-6)

  # a window large enough to take the mmap path must agree too
  bigv <- rnorm(40000)
  big  <- tempfile(fileext = ".ssff")
  write_test_obj(big, list(tr = matrix(bigv, ncol = 1)), "REAL32")
  on.exit(unlink(big), add = TRUE)
  win <- read_ssff(big, begin = 100, end = 300)    # records 10000..30000
  expect_equal(as.numeric(win$tr), bigv[10001:30001], tolerance = 1e-6)
  off <- ssff_data_offset(big)
  expect_equal(as.numeric(win$tr),
               ssff_ref_values(big, 20001L, from = off + 10000L * 4L, size = 4L, float = TRUE),
               tolerance = 0)
})

test_that("threads = 4 produces identical results to threads = 1", {
  set.seed(14)
  m <- matrix(rnorm(500 * 17), 500, 17)
  path <- tempfile(fileext = ".ssff")
  write_test_obj(path, list(tr = m), "REAL32")
  on.exit(unlink(path), add = TRUE)

  one  <- read_ssff(path, threads = 1L)
  four <- read_ssff(path, threads = 4L)
  expect_identical(unname(one$tr), unname(four$tr))
})

test_that("tracks = selects tracks and unknown names are reported", {
  set.seed(15)
  a <- matrix(rnorm(60), 30, 2)
  b <- matrix(rnorm(30), 30, 1)
  c <- matrix(rnorm(30), 30, 1)
  path <- tempfile(fileext = ".ssff")
  write_test_obj(path, list(a = a, b = b, c = c), c("REAL32", "REAL32", "REAL32"))
  on.exit(unlink(path), add = TRUE)

  sel <- read_ssff(path, tracks = c("c", "a"))
  expect_named(sel, c("a", "c"))                 # file order, not request order
  expect_equal(unname(sel$a), a, tolerance = 1e-3)     # REAL32 storage precision
  expect_equal(as.numeric(sel$c), as.numeric(c), tolerance = 1e-3)
  expect_equal(attr(sel, "trackFormats"), c("REAL32", "REAL32"))

  expect_error(read_ssff(path, tracks = "nope"), "not found")
  expect_error(read_ssff(path, tracks = "nope"), "a, b, c")
  expect_error(read_ssff(path, tracks = 42), "character vector")
})

test_that("read_track() maps stored zeros to NA for non-audio tracks only", {
  v <- c(0, 1.5, 0, -2.5, 0, 0, 3)
  path <- tempfile(fileext = ".f0")
  write_test_obj(path, list(f0 = matrix(v, ncol = 1)), "REAL32")
  on.exit(unlink(path), add = TRUE)

  trk <- read_track(path)
  expect_true(all(is.na(trk$f0[c(1, 3, 5, 6)])))
  expect_equal(trk$f0[c(2, 4, 7)], v[c(2, 4, 7)], tolerance = 1e-6)

  # the raw reader still reports the encoding
  raw <- read_ssff(path)
  expect_equal(as.numeric(raw$f0), v, tolerance = 1e-6)
  expect_false(anyNA(raw$f0))

  # opt out explicitly
  expect_equal(as.numeric(read_track(path, zero_to_na = FALSE)$f0), v, tolerance = 1e-6)

  # integer tracks follow the same convention
  ipath <- tempfile(fileext = ".ssff")
  write_test_obj(ipath, list(pm = matrix(c(0, 1, 0, 1), ncol = 1)), "INT16")
  on.exit(unlink(ipath), add = TRUE)
  expect_identical(as.vector(read_track(ipath)$pm), c(NA_integer_, 1L, NA_integer_, 1L))
  expect_identical(as.vector(read_ssff(ipath)$pm), c(0L, 1L, 0L, 1L))
})

test_that("audio tracks keep their zeros (0 is silence, not missing)", {
  wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(wav == "", "test wav not found")

  audio <- read_track(wav)$audio
  expect_false(anyNA(audio))
  expect_identical(audio, read_ssff(wav)$audio)
})

test_that("writing maps NA and NaN to 0 and the round trip is stable", {
  v <- c(NA, 1.25, NaN, -3.5, 0)
  path <- tempfile(fileext = ".f0")
  on.exit(unlink(path), add = TRUE)
  write_test_obj(path, list(f0 = matrix(v, ncol = 1)), "REAL32")

  # on disk: 0 for both NA and NaN
  expect_equal(as.numeric(read_ssff(path)$f0), c(0, 1.25, 0, -3.5, 0), tolerance = 1e-6)
  expect_equal(ssff_ref_values(path, 5L, from = ssff_data_offset(path), size = 4L, float = TRUE),
               c(0, 1.25, 0, -3.5, 0), tolerance = 1e-6)

  # read_track: values survive, non-values become NA
  trk <- read_track(path)
  expect_true(all(is.na(trk$f0[c(1, 3, 5)])))
  expect_equal(trk$f0[2], 1.25, tolerance = 1e-6)
  expect_equal(trk$f0[4], -3.5, tolerance = 1e-6)

  # a second round trip is a fixed point
  path2 <- tempfile(fileext = ".f0")
  on.exit(unlink(path2), add = TRUE)
  write_ssff(trk, path2)
  expect_equal(as.numeric(read_ssff(path2)$f0), c(0, 1.25, 0, -3.5, 0))
})

test_that("a header with many generic variables parses to the same data", {
  set.seed(17)
  v <- rnorm(50)
  path <- tempfile(fileext = ".ssff")
  write_test_obj(path, list(tr = matrix(v, ncol = 1)), "REAL32")
  on.exit(unlink(path), add = TRUE)

  # splice 300 generic variables into the header
  raw  <- readBin(path, "raw", file.size(path))
  off  <- ssff_data_offset(path)
  marker <- charToRaw("-----------------\n")
  filler <- charToRaw(paste0(paste(sprintf("META_%03d FLOAT 1 %d.0", 1:300, 1:300),
                                   collapse = "\n"), "\n"))
  head_txt <- raw[1:(off - length(marker))]     # header without the marker line
  body     <- raw[(off + 1L):length(raw)]
  writeBin(c(head_txt, filler, marker, body), path)

  obj <- read_ssff(path)
  expect_equal(as.numeric(obj$tr), v, tolerance = 1e-6)
  gv <- attr(obj, "genericVars")
  expect_false(is.null(gv))
  expect_true(all(sprintf("META_%03d", 1:300) %in% names(gv)))
})
