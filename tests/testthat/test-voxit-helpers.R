# Pins the Voxit C++ helpers against pure-R transcriptions of the original
# (pre-optimisation) algorithms, so speedups cannot silently change results.

# Verbatim R port of calc_lz_complexity.m (Quang Thai), exhaustive branch,
# as embedded in MillerLab-UCDavis/Voxit voxitAnalysis.m. 1-based indexing
# kept identical to the MATLAB source (gs[n] holds gs(n-1) of the paper).
lz_matlab_ref <- function(S, normalize = TRUE) {
  S <- as.logical(S)
  S[is.na(S)] <- TRUE
  N <- length(S)
  gs <- numeric(N + 1)
  gs[1] <- 0
  s <- paste(ifelse(S, "1", "0"), collapse = "")
  sub <- function(a, b) substr(s, a, b)
  found <- function(m, n) grepl(sub(m, n), sub(1, n - 1), fixed = TRUE)
  gs[2] <- 1
  if (N >= 2) for (n in 2:N) {
    ok <- FALSE
    idx <- (gs[n] + 1):n
    for (k in seq_len(ceiling(length(idx) / 2))) {
      mu <- idx[length(idx) - k + 1]
      if (!found(mu, n)) { gs[n + 1] <- mu; ok <- TRUE; break }
      ml <- idx[k]
      if (found(ml, n)) { gs[n + 1] <- ml - 1; ok <- TRUE; break }
      else if (mu == ml + 1) { gs[n + 1] <- ml; ok <- TRUE; break }
    }
    if (!ok) stop("Internal error: could not find eigenvalue")
  }
  h <- numeric(length(gs)); hl <- 1; hp <- 0
  repeat {
    cand <- which(gs[(hp + 2):length(gs)] > hp)
    if (hp + 2 > length(gs) || !length(cand)) break
    hl <- hl + 1; hp <- hp + cand[1]; h[hl] <- hp
  }
  if (h[hl] != N) { hl <- hl + 1; h[hl] <- N }
  if (gs[h[hl] - 1 + 1] > h[hl - 1]) stop("Check failed for exhaustive sequence")
  C <- hl - 1
  if (normalize) C / (N / log2(N)) else C
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

test_that("lz_complexity_cpp matches Voxit's calc_lz_complexity (LZ76)", {
  lz <- superassp:::lz_complexity_cpp
  # Lempel & Ziv (1976) / Kaspar & Schuster (1987) worked example:
  # 0 | 001 | 10 | 100 | 1000 | 101 -> 6 components
  ex <- as.integer(strsplit("0001101001000101", "")[[1]])
  expect_identical(lz(ex, normalize = FALSE), 6)

  set.seed(20261020)
  seqs <- list(logical(0), TRUE, FALSE, c(TRUE, FALSE), rep(TRUE, 40),
               rep(c(TRUE, FALSE), 30), c(TRUE, NA, FALSE, TRUE),
               sample(0:3, 60, TRUE), rep(runif(12) < 0.5, each = 9))
  for (n in c(2, 3, 5, 17, 64, 150)) for (p in c(0.1, 0.5, 0.8)) {
    seqs[[length(seqs) + 1L]] <- runif(n) < p
  }
  for (s in seqs[-1]) {
    expect_identical(lz(s, normalize = FALSE), lz_matlab_ref(s, FALSE))
    expect_equal(lz(s, normalize = TRUE), 100 * lz_matlab_ref(s, TRUE))
  }
  expect_identical(lz(logical(0)), 0)
  expect_error(lz(c(TRUE, FALSE), type = "primitive"), "Unknown type")
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

test_that(".voxit_rhythmic_complexity matches Voxit's GentleComplexityAllPauses", {
  # Literal transcription of voxitAnalysis.m lines 281-286, 313
  gsps_ref <- function(starts, ends) {
    fsG <- 100
    g <- rep(1, round(ends[length(ends)] * fsG))
    d <- starts[-1] - ends[-length(ends)]
    for (i in which(d >= .1 & d <= 3)) {
      a <- round(ends[i] * fsG) + 1
      b <- round(starts[i + 1] * fsG) - 1
      if (b >= a) g[a:b] <- 0
    }
    100 * lz_matlab_ref(g[(starts[1] * fsG + 1):length(g)])
  }
  st <- c(0.53, 1.20, 1.71, 3.90, 4.12); en <- c(1.02, 1.55, 2.40, 4.05, 4.80)
  expect_equal(superassp:::.voxit_rhythmic_complexity(st, en, 4.80), gsps_ref(st, en))
  st <- c(0, 0.4, 1.1); en <- c(0.3, 0.9, 1.5)
  expect_equal(superassp:::.voxit_rhythmic_complexity(st, en, 1.5), gsps_ref(st, en))
})
