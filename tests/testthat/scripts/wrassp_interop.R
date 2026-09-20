#!/usr/bin/env Rscript
#
# wrassp side of the cross-package SSFF parity checks (see
# helper-wrassp-interop.R for why this runs in a subprocess).
#
# Usage: Rscript wrassp_interop.R <input wav> <output dir>
#
# Writes every file it produces plus an RDS bundle containing
#   files          - named character vector of the generated SSFF files
#   expected       - named list: wrassp::read.AsspDataObj() of each file
#   windows        - list of list(begin, end, obj): windowed wrassp reads
#   deltrack       - wrassp::delTrack() applied to the two-track fixture
#   writers        - named list of synthetic AsspDataObj to write on both sides
#   writer_files   - named character vector: wrassp-written bytes for those
#
# The caller (test-ssff-wrassp-interop.R) compares superassp's reader and
# writer against this bundle.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("usage: wrassp_interop.R <input wav> <output dir>")
wav <- args[[1]]
out <- args[[2]]
stopifnot(file.exists(wav))
dir.create(out, recursive = TRUE, showWarnings = FALSE)

suppressMessages(library(wrassp))

## ------------------------------------------------------------ analyses

producers <- list(
  f0_ksv       = function(o) ksvF0(wav, toFile = TRUE, outputDirectory = o),
  f0_mhs       = function(o) mhsF0(wav, toFile = TRUE, outputDirectory = o),
  rms          = function(o) rmsana(wav, toFile = TRUE, outputDirectory = o),
  spectrum_dft = function(o) dftSpectrum(wav, toFile = TRUE, outputDirectory = o,
                                         beginTime = 0, endTime = 0.3, windowShift = 20),
  acf48        = function(o) acfana(wav, toFile = TRUE, outputDirectory = o,
                                    beginTime = 0, endTime = 0.3),
  forest       = function(o) forest(wav, toFile = TRUE, outputDirectory = o,
                                    beginTime = 0, endTime = 0.3)
)

files <- character()
for (name in names(producers)) {
  dir <- file.path(out, "runs", name)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  suppressMessages(producers[[name]](dir))
  produced <- list.files(dir, full.names = TRUE)
  if (length(produced) != 1L) stop("expected one output for '", name, "', got ", length(produced))
  files[[name]] <- produced[[1]]
}

expected <- lapply(files, function(f) suppressMessages(read.AsspDataObj(f)))

## ------------------------------------------------------------- windows

window_ranges <- list(c(0.5, 1.5), c(2.0, 3.0))
windows <- lapply(window_ranges, function(w) {
  list(begin = w[[1]], end = w[[2]],
       obj = suppressMessages(read.AsspDataObj(files[["f0_ksv"]], begin = w[[1]], end = w[[2]])))
})

## ------------------------------------------------------------- delTrack

deltrack <- delTrack(expected[["forest"]], "bw")

# adding a *new* track: compared against superassp::addTrack() by the caller
addtrack <- addTrack(expected[["forest"]], "extra",
                     matrix(seq_len(nrow(expected[["forest"]]$fm)), ncol = 1), "REAL32")

## ------------------------------------------------- writer byte parity

mk <- function(mats, formats) {
  obj <- mats
  attr(obj, "trackFormats") <- formats
  attr(obj, "sampleRate")   <- 100
  attr(obj, "origFreq")     <- 44100
  attr(obj, "startTime")    <- 0
  attr(obj, "startRecord")  <- 1L
  attr(obj, "endRecord")    <- nrow(mats[[1]])
  attr(obj, "fileInfo")     <- c(20L, 2L)
  class(obj) <- "AsspDataObj"
  obj
}

set.seed(99)
writers <- list(
  f32   = mk(list(f0 = matrix(round(rnorm(24), 5), ncol = 1)), "REAL32"),
  f64   = mk(list(fm = matrix(rnorm(32), ncol = 4)), "REAL64"),
  i16   = mk(list(pm = matrix(sample(0:1, 24, TRUE), ncol = 1)), "INT16"),
  multi = mk(list(a = matrix(round(rnorm(24), 3), ncol = 1),
                  b = matrix(rnorm(48), ncol = 2)), c("REAL32", "REAL64"))
)

writer_files <- character()
for (name in names(writers)) {
  path <- file.path(out, paste0("writer_", name, ".ssff"))
  invisible(suppressMessages(write.AsspDataObj(writers[[name]], path)))
  writer_files[[name]] <- path
}

saveRDS(list(files = files, expected = expected, windows = windows,
             deltrack = deltrack, addtrack = addtrack,
             writers = writers, writer_files = writer_files),
        file.path(out, "interop.rds"))
cat("interop bundle written\n")
