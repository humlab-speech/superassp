#!/usr/bin/env Rscript
# Convert MATLAB-generated golden .mat files to .rds for testthat consumption.

suppressPackageStartupMessages({ library(R.matlab) })

golden_dir <- "tests/testthat/golden"

# Main goldens
mat <- readMat(file.path(golden_dir, "a1.mat"))
g <- list(
  x        = as.numeric(mat$x),
  fs       = as.numeric(mat$fs),
  f0       = as.numeric(mat$f0),
  VUV      = as.integer(mat$VUV),
  SRHVal   = as.numeric(mat$SRHVal),
  t_f0     = as.numeric(mat$t.f0),
  creak_pp = as.integer(mat$creak.pp),
  creak_post = as.numeric(mat$creak.post),
  creak_dec  = as.integer(mat$creak.dec),
  creak_FeatTot = mat$creak.FeatTot,
  H2H1     = as.numeric(mat$H2H1),
  res_p    = as.numeric(mat$res.p),
  GCI      = as.integer(mat$GCI),
  GCI_varF0 = as.integer(mat$GCI.varF0),
  f0_samp  = as.numeric(mat$f0.samp),
  rep      = as.numeric(mat$rep),
  res      = as.numeric(mat$res),
  MBS      = as.numeric(mat$MBS),
  g_iaif   = as.numeric(mat$g.iaif),
  ar_lpc   = mat$ar.lpc,
  e_lpc    = as.numeric(mat$e.lpc),
  NAQ      = as.numeric(mat$NAQ),
  QOQ      = as.numeric(mat$QOQ),
  H1H2_param = as.numeric(mat$H1H2.param),
  HRF      = as.numeric(mat$HRF),
  ps       = as.numeric(mat$ps),
  mdq      = as.numeric(mat$mdq)
)
saveRDS(g, file.path(golden_dir, "a1.rds"))
cat("Wrote", file.path(golden_dir, "a1.rds"), "\n")

# ANN weights
ann_path <- file.path(golden_dir, "creak_ann_flat.mat")
if (file.exists(ann_path)) {
  ann_mat <- readMat(ann_path)
  ann <- list(
    IW   = ann_mat$IW,
    LW   = ann_mat$LW,
    b_h  = as.numeric(ann_mat$b.h),
    b_o  = as.numeric(ann_mat$b.o),
    mini = as.numeric(ann_mat$Minis),
    maxi = as.numeric(ann_mat$Maxis),
    out_xmin = as.numeric(ann_mat$out.xmin),
    out_xmax = as.numeric(ann_mat$out.xmax),
    out_ymin = as.numeric(ann_mat$out.ymin),
    out_ymax = as.numeric(ann_mat$out.ymax)
  )
  dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
  saveRDS(ann, "inst/extdata/creak_ann.rds")
  cat("Wrote inst/extdata/creak_ann.rds (", nrow(ann$IW), "hidden,",
      length(ann$mini), "features)\n")
}
