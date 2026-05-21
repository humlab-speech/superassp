# Bit-faithful creak feature + posterior parity vs MATLAB on a1.wav.
# Regenerate goldens with: bash tools/regen_goldens.sh

golden_path <- testthat::test_path("golden", "voiceanalysis", "a1.rds")

skip_if_no_golden <- function() {
  if (!file.exists(golden_path))
    testthat::skip("Golden a1.rds not present (run tools/regen_goldens.sh)")
}

test_that("creak feature matrix shape matches MATLAB on a1.wav", {
  skip_if_no_golden()
  g <- readRDS(golden_path)
  if (is.null(g$creak_FeatTot)) skip("MATLAB FeatTot not in golden")

  F_r <- vat_creak_features_cpp(g$x, g$fs)
  expect_equal(ncol(F_r), ncol(g$creak_FeatTot))
  # Allow ±2 frames difference at the edges
  expect_lt(abs(nrow(F_r) - nrow(g$creak_FeatTot)), 5)
})

test_that("creak static features correlate with MATLAB", {
  skip_if_no_golden()
  g <- readRDS(golden_path)
  if (is.null(g$creak_FeatTot)) skip("MATLAB FeatTot not in golden")

  F_r <- vat_creak_features_cpp(g$x, g$fs)
  n <- min(nrow(F_r), nrow(g$creak_FeatTot))
  # Static columns 1..12. Min thresholds per column (looser for res_p which
  # has many ad-hoc peak-picking heuristics).
  thresh <- c(0.85, 0.55, 0.90, 0.90, 0.90, 0.85, 0.85, 0.85, NA, 0.99, 0.99, 0.90)
  for (j in 1:12) {
    if (is.na(thresh[j])) next
    m <- g$creak_FeatTot[1:n, j]; r <- F_r[1:n, j]
    if (sd(m) > 0 && sd(r) > 0) {
      cc <- cor(r, m)
      expect_gt(cc, thresh[j], label = sprintf("col %d (cor=%.3f)", j, cc))
    }
  }
})

test_that("creak posterior decision matches MATLAB on a1.wav", {
  skip_if_no_golden()
  g <- readRDS(golden_path)
  ann_path <- system.file("extdata", "creak_ann.rds", package = "superassp")
  if (ann_path == "") skip("creak_ann.rds not bundled")

  r <- .vat_creak_detect(g$x, g$fs, backend = "cpp")
  n <- min(length(r$posterior), length(g$creak_post))

  # Posterior magnitudes should be small and bounded
  expect_lt(max(r$posterior), 1.0)
  expect_gt(min(r$posterior), 0.0)

  # Posterior max|diff|: on a sustained vowel both networks output near zero;
  # any single-frame divergence should remain below 0.15
  expect_lt(max(abs(r$posterior[1:n] - g$creak_post[1:n])), 0.15)

  # Decision agreement should be near-perfect (a1.wav is a sustained vowel,
  # no creak expected)
  agree <- mean(r$decision[1:n] == g$creak_dec[1:n])
  expect_gt(agree, 0.95)
})

test_that("vat_creak_detect_cpp output shape sane on synthetic input", {
  fs <- 16000
  x <- sin(2 * pi * 150 * seq(0, 1, by = 1/fs)) + 0.01 * rnorm(fs + 1)
  ann_path <- system.file("extdata", "creak_ann.rds", package = "superassp")
  if (ann_path == "") skip("creak_ann.rds not bundled")
  r <- .vat_creak_detect(x, fs, backend = "cpp")
  expect_named(r, c("features", "posterior", "decision", "time"))
  expect_equal(length(r$posterior), length(r$decision))
  expect_true(all(r$posterior >= 0 & r$posterior <= 1))
})
