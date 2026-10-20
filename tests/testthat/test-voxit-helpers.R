# Pins the Voxit C++ helpers against pure-R transcriptions of the original
# (pre-optimisation) algorithms, so speedups cannot silently change results.

# Literal transcription of the original lz_complexity_cpp (string dictionary,
# substring rebuilt for every position and length).
lz_reference_r <- function(S, type = "exhaustive", normalize = TRUE) {
  S <- as.logical(S)
  n <- length(S)
  if (n < 1) return(0)
  ch <- ifelse(is.na(S) | S, "1", "0")
  sub <- function(i, len) paste0(ch[i + seq_len(len)], collapse = "")  # i 0-based
  dict <- character(0)
  c <- 0
  i <- 0L
  if (type == "exhaustive") {
    while (i < n) {
      maxlen <- 0L
      len <- 1L
      while (len <= n - i && len <= i) {
        if (sub(i, len) %in% dict) maxlen <- len
        len <- len + 1L
      }
      c <- c + 1
      if (maxlen > 0) {
        for (j in seq_len(maxlen)) dict <- union(dict, sub(i, j))
        i <- i + maxlen
      } else {
        dict <- union(dict, ch[i + 1L])
        i <- i + 1L
      }
    }
  } else if (type == "primitive") {
    while (i < n) {
      maxlen <- 0L
      for (len in seq_len(n - i)) {
        if (!(sub(i, len) %in% dict)) maxlen <- len
      }
      if (maxlen > 0) {
        c <- c + 1
        dict <- union(dict, sub(i, maxlen))
        i <- i + maxlen
      } else {
        i <- i + 1L
      }
    }
  } else {
    stop("Unknown type")
  }
  if (normalize) return(100 * c / (n / log2(n)))
  c
}

# Literal transcription of the original interp1_linear_cpp (linear scan).
interp1_reference_r <- function(xp, fp, xi) {
  m <- length(xp)
  vapply(xi, function(x) {
    idx <- 0L
    while (idx < m - 1L && isTRUE(xp[idx + 2L] <= x)) idx <- idx + 1L
    if (idx == m - 1L) {
      fp[m]
    } else if (idx == 0L && isTRUE(x < xp[1L])) {
      fp[1L]
    } else {
      x0 <- xp[idx + 1L]; x1 <- xp[idx + 2L]
      y0 <- fp[idx + 1L]; y1 <- fp[idx + 2L]
      y0 + ((x - x0) / (x1 - x0)) * (y1 - y0)
    }
  }, numeric(1))
}

test_that("lz_complexity_cpp matches the reference algorithm", {
  lz <- superassp:::lz_complexity_cpp
  set.seed(20261020)
  seqs <- list(logical(0), TRUE, FALSE, c(TRUE, FALSE), rep(TRUE, 40),
               rep(c(TRUE, FALSE), 30), c(TRUE, NA, FALSE, TRUE),
               sample(0:3, 60, TRUE))
  for (n in c(2, 5, 17, 64, 150)) for (p in c(0.2, 0.5, 0.8)) {
    seqs[[length(seqs) + 1L]] <- runif(n) < p
  }
  for (s in seqs) for (ty in c("exhaustive", "primitive")) for (nz in c(TRUE, FALSE)) {
    expect_identical(lz(s, type = ty, normalize = nz),
                     lz_reference_r(s, type = ty, normalize = nz))
  }
  expect_error(lz(c(TRUE, FALSE), type = "bogus"), "Unknown type")
})

test_that("interp1_linear_cpp matches the reference linear scan", {
  ip <- superassp:::interp1_linear_cpp
  set.seed(7)
  xp <- sort(runif(40) * 10); fp <- rnorm(40)
  cases <- list(
    list(xp, fp, c(-1, 0, xp, 10, 11, NaN, Inf, -Inf, runif(50) * 10)),
    list(c(1, 2, 2, 3), c(10, 20, 30, 40), c(0, 1, 1.5, 2, 2.5, 3, 4)),
    list(c(3, 1, 2, 5), c(1, 2, 3, 4), c(0, 1, 2, 2.5, 3, 4, 5, 6)),
    list(c(1, NA, 3), c(1, 2, 3), c(0, 1, 2, 3, 4)),
    list(5, 7, c(4, 5, 6)),
    list(runif(60), rnorm(60), runif(100)))
  for (a in cases) {
    # expect_equal, not identical: the compiler may fuse y0 + f * (y1 - y0)
    # into an FMA, which R arithmetic cannot reproduce bit-for-bit.
    expect_equal(ip(a[[1]], a[[2]], a[[3]]),
                 interp1_reference_r(a[[1]], a[[2]], a[[3]]))
  }
})

test_that("contiguous_runs_cpp returns 1-indexed runs", {
  cr <- superassp:::contiguous_runs_cpp
  expect_identical(cr(integer(0)),
                   list(starts = integer(0), stops = integer(0), lengths = integer(0)))
  expect_identical(cr(c(1L, 1L, 0L, 1L, 2L, 1L)),
                   list(starts = c(1L, 4L, 6L), stops = c(2L, 4L, 6L),
                        lengths = c(2L, 1L, 1L)))
  set.seed(3)
  x <- sample(0:2, 300, TRUE)
  r <- rle(x)
  ends <- cumsum(r$lengths)
  sel <- r$values == 1L
  expect_identical(cr(x),
                   list(starts = as.integer(ends[sel] - r$lengths[sel] + 1L),
                        stops = as.integer(ends[sel]),
                        lengths = as.integer(r$lengths[sel])))
})
