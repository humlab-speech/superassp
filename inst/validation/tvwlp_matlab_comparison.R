# TVWLP Formant Tracker — MATLAB vs Rcpp Validation Harness
#
# Compares trk_formant_tvwlp() against the MATLAB reference implementation
# in ../ftrack/ftrack_tvwlp_v1.
#
# Assumptions:
#   - sibling repo at <superassp>/../ftrack/ftrack_tvwlp_v1 with run_matlab_test.m
#   - matlab_output.mat present (or MATLAB available to regenerate it)
#
# Usage:
#   Rscript inst/validation/tvwlp_matlab_comparison.R
#
# Outputs:
#   inst/validation/tvwlp_validation_report.md
#   inst/validation/tvwlp_validation_data.rds

suppressPackageStartupMessages({
  devtools::load_all(here::here())
  library(R.matlab)
})

ftrack_dir <- normalizePath(file.path(here::here(), "..", "ftrack", "ftrack_tvwlp_v1"),
                            mustWork = FALSE)
out_dir    <- file.path(here::here(), "inst", "validation")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stopifnot(dir.exists(ftrack_dir))

wav_timit  <- file.path(ftrack_dir, "test-dr8-mjln0-sx279.wav")
mat_tvwlp  <- file.path(ftrack_dir, "matlab_output.mat")
mat_tvlp   <- file.path(ftrack_dir, "matlab_tvlp_l2.mat")

# ---------------- Comparison helpers ----------------

compare_one <- function(Fi_m, Fi_r, label) {
  n <- min(nrow(Fi_m), nrow(Fi_r))
  Fi_m <- Fi_m[seq_len(n), , drop = FALSE]
  Fi_r <- Fi_r[seq_len(n), , drop = FALSE]
  out <- vector("list", 3)
  for (i in 1:3) {
    ok <- Fi_m[, i] > 0 & Fi_r[, i] > 0 & is.finite(Fi_m[, i]) & is.finite(Fi_r[, i])
    xv <- Fi_m[ok, i]; yv <- Fi_r[ok, i]; d <- yv - xv
    out[[i]] <- data.frame(
      label      = label, formant = sprintf("F%d", i),
      n_valid    = sum(ok), n_total = length(ok),
      mean_matlab = mean(xv), mean_r = mean(yv),
      corr       = cor(xv, yv),
      rmse       = sqrt(mean(d^2)),
      mae        = mean(abs(d)),
      bias       = mean(d),
      within_50  = 100 * mean(abs(d) < 50),
      within_100 = 100 * mean(abs(d) < 100)
    )
  }
  do.call(rbind, out)
}

# ---------------- Run Rcpp port ----------------

cat("Running Rcpp tvwlp_l2 (weighted) ...\n")
res_tvwlp <- trk_formant_tvwlp(wav_timit, lptype = "tvwlp_l2",
                                toFile = FALSE, verbose = FALSE)

cat("Running Rcpp tvlp_l2 (unweighted) ...\n")
res_tvlp  <- trk_formant_tvwlp(wav_timit, lptype = "tvlp_l2",
                                toFile = FALSE, verbose = FALSE)

# ---------------- Load MATLAB reference ----------------

stopifnot(file.exists(mat_tvwlp))
m1 <- readMat(mat_tvwlp)
Fi_m_tvwlp <- t(m1$Fi.matlab)

stopifnot(file.exists(mat_tvlp))
m2 <- readMat(mat_tvlp)
Fi_m_tvlp  <- t(m2$Fi)

# ---------------- Compare ----------------

df_tvwlp <- compare_one(Fi_m_tvwlp, res_tvwlp$fm, "tvwlp_l2 (weighted)")
df_tvlp  <- compare_one(Fi_m_tvlp,  res_tvlp$fm,  "tvlp_l2 (unweighted)")
results  <- rbind(df_tvwlp, df_tvlp)

print(results, row.names = FALSE)

saveRDS(list(results = results,
             rcpp_tvwlp = res_tvwlp, rcpp_tvlp = res_tvlp,
             matlab_tvwlp = Fi_m_tvwlp, matlab_tvlp = Fi_m_tvlp),
        file.path(out_dir, "tvwlp_validation_data.rds"))

# ---------------- Markdown report ----------------

report_path <- file.path(out_dir, "tvwlp_validation_report.md")
fmt_row <- function(r) sprintf(
  "| %s | %s | %d/%d | %.1f | %.1f | %.4f | %.1f | %.1f | %+.1f | %.1f%% |",
  r$label, r$formant, r$n_valid, r$n_total,
  r$mean_matlab, r$mean_r, r$corr, r$rmse, r$mae, r$bias, r$within_50)

lines <- c(
  "# TVWLP Rcpp vs MATLAB Validation Report",
  "",
  sprintf("Generated: %s", format(Sys.time())),
  sprintf("Test audio: `%s` (16 kHz, 3.59 s TIMIT sentence)", basename(wav_timit)),
  "",
  "| LP method | Formant | n_valid | MATLAB mean Hz | R mean Hz | corr | RMSE Hz | MAE Hz | bias Hz | within 50 Hz |",
  "|---|---|---|---|---|---|---|---|---|---|",
  vapply(seq_len(nrow(results)), function(i) fmt_row(results[i, ]), character(1)),
  "",
  "## Interpretation",
  "",
  "Post-fix state (May 2026). Three bugs fixed against MATLAB reference:",
  "1. **Levinson-Durbin sign error** in `levinson_durbin()` (tvwlp_core.cpp):",
  "   reflection coefficient + predictor update had wrong signs. LPC residual",
  "   correlation jumped from 0.16 to 0.998 against MATLAB `lpc()`.",
  "2. **Wrong Blackman window** in SRH and SEDREAMS: `av::blackman()` is the",
  "   periodic DFT-even variant; MATLAB uses the symmetric form. GCI count",
  "   dropped from 581 (231 spurious) to 349, matching MATLAB's 350.",
  "3. **Wrong Hann window** convention in get_lpc_residual_cpp (minor).",
  "",
  "Remaining gap vs Python port (0.85-0.95 corr): tvwlp_l2 F1=0.46/F2=0.75/F3=0.69.",
  "Probable cause is Armadillo `arma::solve` behaviour on near-singular per-frame",
  "LP systems (the runtime emits `solve(): system is singular` warnings). F1 in",
  "tvlp_l2 (unweighted) stays at 0.34 even though signal correlation is 0.9999,",
  "suggesting the residual gap is in the LP solver, not the preprocessing.",
  "",
  "## Caveats",
  "",
  "- MATLAB does not output bandwidths; superassp's `bw` track cannot be ground-truthed here.",
  "- R produces 358 frames vs MATLAB 359 (one-frame boundary difference).",
  "- FFmpeg vs MATLAB resampling: signal correlation 0.9999, negligible contribution."
)
writeLines(lines, report_path)
cat("Report written to:", report_path, "\n")
