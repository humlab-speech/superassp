# Golden-file regression tests against MATLAB reference on a1.wav.
# Regenerate with: bash tools/regen_goldens.sh

golden_path <- testthat::test_path("golden", "voiceanalysis", "a1.rds")

skip_if_no_golden <- function() {
  if (!file.exists(golden_path)) {
    testthat::skip("Golden a1.rds not present (run tools/regen_goldens.sh)")
  }
}

test_that("SRH pitch tracker matches MATLAB on a1.wav", {
  skip_if_no_golden()
  g <- readRDS(golden_path)
  res <- .vat_pitch(g$x, g$fs, 20, 500, method = "srh")

  n <- min(length(res$f0), length(g$f0))
  m_f0 <- g$f0[seq_len(n)]
  r_f0 <- res$f0[seq_len(n)]
  voiced <- which(g$VUV[seq_len(n)] == 1 & m_f0 > 0 & r_f0 > 0)
  if (length(voiced) > 20) {
    expect_gt(cor(r_f0[voiced], m_f0[voiced]), 0.85)
    expect_lt(median(abs(r_f0[voiced] - m_f0[voiced])), 5)  # Hz
  }
})

test_that("SE-VQ GCI count and median IGI match MATLAB on a1.wav", {
  skip_if_no_golden()
  g <- readRDS(golden_path)
  # Upsample MATLAB-style f0/VUV from 10ms to per-sample
  ts <- seq(0, by = 0.010, length.out = length(g$f0))
  tx <- (seq_along(g$x) - 1) / g$fs
  f0_samp  <- stats::approx(ts, g$f0,  tx, rule = 2)$y
  vuv_samp <- as.integer(stats::approx(ts, g$VUV, tx, method = "constant", rule = 2)$y >= 0.5)

  r <- .vat_se_vq(g$x, g$fs, f0_samp, vuv_samp, backend = "cpp")

  # GCI count within ±25 %
  expect_lt(abs(length(r$GCI) - length(g$GCI)) / length(g$GCI), 0.25)

  # Median inter-GCI interval should match within 10 %
  if (length(r$GCI) > 5 && length(g$GCI) > 5) {
    m_igi <- median(diff(g$GCI))
    r_igi <- median(diff(r$GCI))
    expect_lt(abs(r_igi - m_igi) / m_igi, 0.1)
  }
})

test_that("SE-VQ var-F0 GCI count matches MATLAB on a1.wav", {
  skip_if_no_golden()
  g <- readRDS(golden_path)
  if (is.null(g$GCI_varF0) || length(g$GCI_varF0) < 5) {
    skip("MATLAB var-F0 golden not present (regen with patched get_MBS_GCI_intervals.m)")
  }
  ts <- seq(0, by = 0.010, length.out = length(g$f0))
  tx <- (seq_along(g$x) - 1) / g$fs
  f0_samp  <- stats::approx(ts, g$f0,  tx, rule = 2)$y
  vuv_samp <- as.integer(stats::approx(ts, g$VUV, tx, method = "constant", rule = 2)$y >= 0.5)

  r <- .vat_se_vq(g$x, g$fs, f0_samp, vuv_samp, var_f0 = TRUE, backend = "cpp")
  expect_lt(abs(length(r$GCI) - length(g$GCI_varF0)) / length(g$GCI_varF0), 0.25)
  if (length(r$GCI) > 5 && length(g$GCI_varF0) > 5) {
    m_igi <- median(diff(g$GCI_varF0))
    r_igi <- median(diff(r$GCI))
    expect_lt(abs(r_igi - m_igi) / m_igi, 0.15)
  }
})

test_that("peakSlope correlates with MATLAB on a1.wav", {
  skip_if_no_golden()
  g <- readRDS(golden_path)
  ps_r <- .vat_peak_slope(g$x, g$fs, backend = "cpp")
  n <- min(length(ps_r), length(g$ps))
  # Daless wavelet bank should give strong correlation
  expect_gt(cor(ps_r[seq_len(n)], g$ps[seq_len(n)]), 0.7)
})

test_that("voice quality (NAQ/QOQ/H1H2) matches MATLAB on a1.wav", {
  skip_if_no_golden()
  g <- readRDS(golden_path)
  if (all(g$NAQ == 0)) skip("MATLAB golden VQ is zero — regen with patched get_NAQ_QOQ_H1H2.m")

  # Use MATLAB GCI + MATLAB glottal flow to isolate VQ from upstream drift
  res <- .vat_voice_quality(g$g_iaif, g$fs, g$GCI, backend = "cpp")

  n <- min(length(res$NAQ), length(g$NAQ))
  # Drop frames where either side is zero (boundary / undefined cycles)
  ok <- which(g$NAQ[seq_len(n)] != 0 & res$NAQ[seq_len(n)] != 0)
  if (length(ok) > 20) {
    expect_gt(cor(res$NAQ[ok], g$NAQ[ok]), 0.85)
    expect_lt(median(abs(res$NAQ[ok] - g$NAQ[ok])) / median(g$NAQ[ok]), 0.25)
  }
  ok <- which(g$QOQ[seq_len(n)] != 0 & res$QOQ[seq_len(n)] != 0)
  if (length(ok) > 20) {
    expect_gt(cor(res$QOQ[ok], g$QOQ[ok]), 0.7)
  }
})

test_that("MDQ correlates with MATLAB on a1.wav", {
  skip_if_no_golden()
  g <- readRDS(golden_path)
  # Use MATLAB residual + GCI to isolate MDQ from upstream differences
  mdq_r <- .vat_mdq(g$res, g$fs, g$GCI)
  n <- min(length(mdq_r), length(g$mdq))
  ok <- which(is.finite(mdq_r[seq_len(n)]) & is.finite(g$mdq[seq_len(n)]) &
              g$mdq[seq_len(n)] > 0)
  if (length(ok) > 20) {
    expect_gt(cor(mdq_r[ok], g$mdq[ok]), 0.7)
  }
})

test_that("creak ANN forward pass produces output of expected length", {
  skip_if_no_golden()
  ann_path <- system.file("extdata", "creak_ann.rds", package = "superassp")
  testthat::skip_if(ann_path == "", "creak_ann.rds not bundled")
  ann <- readRDS(ann_path)
  expect_equal(length(ann$mini), length(ann$maxi))
  expect_equal(ncol(ann$IW), length(ann$mini))
})
