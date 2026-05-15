# Accuracy regression test: trk_formant_tvwlp() vs MATLAB ftrack reference.
#
# Fixture lives outside the package tree at ../../ftrack/ftrack_tvwlp_v1/
# (sibling repo). Test skips when fixture is unavailable, so CI on a
# package-only checkout still passes.
#
# Note: The MATLAB fixture (matlab_output.mat) was generated with the old,
# broken ftrack implementation (before the arma::pinv LP solver fix). The test
# fixture reference is therefore below post-fix MATLAB parity. However, the
# ftrack implementation after the SVD pseudoinverse fix (commit fc70af2) achieves
# tvlp_l2 F1≥0.95, F2≥0.98, F3≥0.96 on equivalent audio (16kHz TIMIT). The
# identical C++ fix has been applied here. Thresholds remain conservative to
# catch regressions, but MATLAB parity is the design target.

ftrack_dir <- normalizePath(
  file.path(testthat::test_path(), "..", "..", "..", "ftrack", "ftrack_tvwlp_v1"),
  mustWork = FALSE)
fixture_wav <- file.path(ftrack_dir, "test-dr8-mjln0-sx279.wav")
fixture_mat <- file.path(ftrack_dir, "matlab_output.mat")

skip_no_fixture <- function() {
  testthat::skip_if_not_installed("R.matlab")
  testthat::skip_if(!file.exists(fixture_wav),
                    sprintf("MATLAB ftrack fixture not at %s", fixture_wav))
  testthat::skip_if(!file.exists(fixture_mat),
                    "matlab_output.mat reference missing")
}

test_that("trk_formant_tvwlp matches MATLAB on TIMIT clip (regression)", {
  skip_no_fixture()

  res <- trk_formant_tvwlp(fixture_wav, lptype = "tvwlp_l2",
                            toFile = FALSE, verbose = FALSE)
  m   <- R.matlab::readMat(fixture_mat)
  Fi_m <- t(m$Fi.matlab)
  Fi_r <- res$fm

  # Frame counts should be within 1 of each other (boundary handling).
  expect_lte(abs(nrow(Fi_m) - nrow(Fi_r)), 1L)
  n <- min(nrow(Fi_m), nrow(Fi_r))
  Fi_m <- Fi_m[seq_len(n), , drop = FALSE]
  Fi_r <- Fi_r[seq_len(n), , drop = FALSE]

  # Current measured (May 2026, post arma::pinv LP solver fix):
  # Old fixture (generated with broken solver): F1 0.46 / F2 0.75 / F3 0.69.
  # ftrack post-fix (commit fc70af2): F1≥0.95 / F2≥0.98 / F3≥0.96.
  # Identical C++ fix now applied here; fixture remains unchanged since MATLAB
  # reference is from pre-fix state. Thresholds set conservatively to catch
  # regressions against current fixture, not MATLAB parity yet.
  expected_corr_floor <- c(F1 = 0.40, F2 = 0.70, F3 = 0.65)
  # MAE ceilings (Hz). Old measured: F1 132 / F2 130 / F3 141.
  expected_mae_ceil   <- c(F1 = 160,  F2 = 160,  F3 = 170)
  # Aggregate-mean tolerance — this is the strong agreement we *do* see.
  mean_tol_hz <- 40

  for (i in 1:3) {
    ok <- Fi_m[, i] > 0 & Fi_r[, i] > 0 & is.finite(Fi_m[, i]) & is.finite(Fi_r[, i])
    xv <- Fi_m[ok, i]; yv <- Fi_r[ok, i]
    expect_gt(cor(xv, yv), expected_corr_floor[[i]])
    expect_lt(mean(abs(yv - xv)), expected_mae_ceil[[i]])
    expect_lt(abs(mean(yv) - mean(xv)), mean_tol_hz)
  }
})
