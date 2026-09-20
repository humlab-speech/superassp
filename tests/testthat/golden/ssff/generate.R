#!/usr/bin/env Rscript
#
# Generate the SSFF gold fixtures used by test-ssff-wrassp-golden.R.
#
# The fixtures are produced by wrassp -- the reference implementation of the
# very SSFF reader/writer that superassp carries -- from the bundled sample
# recording. Next to every fixture this script stores the object that
# wrassp::read.AsspDataObj() returns for it, so the tests can compare
# superassp's reader against an external expectation without loading wrassp
# into the test session.
#
# Usage:
#   Rscript tests/testthat/golden/ssff/generate.R           # write (refuses to overwrite)
#   Rscript tests/testthat/golden/ssff/generate.R --force   # overwrite
#   Rscript tests/testthat/golden/ssff/generate.R --check   # report drift, write nothing
#
# Requires the wrassp package (Suggests). Does NOT require superassp to be
# installed: the input recording is located relative to this script.

args    <- commandArgs(trailingOnly = TRUE)
force   <- "--force" %in% args
check   <- "--check" %in% args

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
here        <- normalizePath(dirname(script_file[1]))
pkg_root    <- normalizePath(file.path(here, "..", "..", "..", ".."))
wav         <- file.path(pkg_root, "inst", "samples", "sustained", "a1.wav")

if (!requireNamespace("wrassp", quietly = TRUE)) {
  stop("wrassp is required to generate the fixtures (install.packages('wrassp'))")
}
if (!file.exists(wav)) {
  stop("sample recording not found: ", wav)
}

md5 <- function(path) unname(tools::md5sum(path))
msg <- function(...) cat(..., "\n", sep = "")
suppressMessages(library(wrassp))

## ---------------------------------------------------------------- generation

work <- file.path(tempdir(), "ssff-gold")
unlink(work, recursive = TRUE)
dir.create(work, recursive = TRUE)

# Every producer gets its own output directory; the single file it writes is
# renamed to <fixture>.ssff.
produce <- function(name, fun) {
  out <- file.path(work, "run", name)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  suppressMessages(fun(out))
  files <- list.files(out, full.names = TRUE)
  if (length(files) != 1L) {
    stop("expected exactly one file for '", name, "', got ", length(files))
  }
  file.rename(files, file.path(work, paste0(name, ".ssff")))
}

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

for (name in names(producers)) {
  produce(name, producers[[name]])
  msg("generated ", name, ".ssff (", file.size(file.path(work, paste0(name, ".ssff"))), " bytes)")
}

# Big-endian derivative of the F0 fixture: the data section of a SSFF file is a
# stream of values in the byte order named by the header's Machine line, so
# reversing every 4-byte group and saying "SPARC" is a valid re-encoding.
be_body <- function(from, to) {
  raw  <- readBin(from, "raw", file.size(from))
  end  <- grepRaw(charToRaw("-----------------\n"), raw, fixed = TRUE) + 17L
  body <- raw[(end + 1L):length(raw)]
  stopifnot(length(body) %% 4L == 0L)
  header <- sub("Machine IBM-PC", "Machine SPARC  ", rawToChar(raw[1:end]), fixed = TRUE)
  writeBin(c(charToRaw(header), as.raw(as.vector(matrix(as.integer(body), nrow = 4)[4:1, ]))), to)
}
be_body(file.path(work, "f0_ksv.ssff"), file.path(work, "f0_ksv_be.ssff"))
msg("generated f0_ksv_be.ssff (big-endian derivative of f0_ksv.ssff)")

# A wrassp-written file whose track contains NA/NaN. wrassp stores those values
# as NaN; superassp's writer stores them as 0. Both encodings are pinned by the
# tests so the documented difference cannot drift unnoticed.
na_obj <- list(f0 = matrix(c(NA, 1.5, NaN, -2), ncol = 1))
attr(na_obj, "trackFormats") <- "REAL32"
attr(na_obj, "sampleRate")   <- 100
attr(na_obj, "origFreq")     <- 44100
attr(na_obj, "startTime")    <- 0
attr(na_obj, "startRecord")  <- 1L
attr(na_obj, "endRecord")    <- 4L
attr(na_obj, "fileInfo")     <- c(20L, 2L)
class(na_obj) <- "AsspDataObj"
invisible(suppressMessages(write.AsspDataObj(na_obj, file.path(work, "na_written_by_wrassp.ssff"))))
msg("generated na_written_by_wrassp.ssff")

## ------------------------------------------------------------- expectations

fixtures <- list.files(work, pattern = "\\.ssff$", full.names = TRUE)
for (fixture in fixtures) {
  expected <- suppressMessages(read.AsspDataObj(fixture))
  # filePath is generation-directory dependent; store the file name instead so
  # the expectations are reproducible and --check can compare them.
  attr(expected, "filePath") <- basename(fixture)
  saveRDS(expected, paste0(tools::file_path_sans_ext(fixture), ".expected.rds"),
          compress = "xz")
}
msg("captured ", length(fixtures), " wrassp::read.AsspDataObj() expectations")

## --------------------------------------------------------------- provenance

provenance <- c(
  "# SSFF gold fixtures -- provenance",
  "",
  "Generated by `tests/testthat/golden/ssff/generate.R`; do not edit by hand.",
  "",
  "| item | value |",
  "|---|---|",
  paste0("| generated | ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), " |"),
  paste0("| R | ", R.version.string, " |"),
  paste0("| platform | ", R.version$platform, " |"),
  paste0("| wrassp | ", as.character(utils::packageVersion("wrassp")), " |"),
  paste0("| input recording | `inst/samples/sustained/a1.wav` (md5 `", md5(wav), "`) |"),
  "",
  "The fixtures are byte-stable: re-running the generator reproduces them",
  "(`generate.R --check` verifies this). The `.expected.rds` files hold the",
  "object that `wrassp::read.AsspDataObj()` returned for the fixture at",
  "generation time -- they are the external reference the tests compare against.",
  "",
  "| fixture | size (bytes) | md5 | producer |",
  "|---|---|---|---|"
)

for (fixture in fixtures) {
  name <- basename(fixture)
  producer <- switch(sub("\\.ssff$", "", name),
    f0_ksv       = "`wrassp::ksvF0()` on the full recording",
    f0_mhs       = "`wrassp::mhsF0()` on the full recording",
    rms          = "`wrassp::rmsana()` on the full recording",
    spectrum_dft = "`wrassp::dftSpectrum(beginTime = 0, endTime = 0.3, windowShift = 20)`",
    acf48        = "`wrassp::acfana(beginTime = 0, endTime = 0.3)`",
    forest       = "`wrassp::forest(beginTime = 0, endTime = 0.3)`",
    f0_ksv_be    = "derived from `f0_ksv.ssff` (4-byte groups reversed, `Machine SPARC`)",
    na_written_by_wrassp = "`wrassp::write.AsspDataObj()` of a track holding `NA`/`NaN`",
    name
  )
  provenance <- c(provenance, sprintf("| `%s` | %d | `%s` | %s |",
                                      name, file.size(fixture), md5(fixture), producer))
}

provenance <- c(provenance, "",
  "Regenerate with `Rscript tests/testthat/golden/ssff/generate.R --force`",
  "after a deliberate change; `--check` reports drift without writing.")

## ------------------------------------------------------------------ install

shipped <- list.files(work, pattern = "\\.(ssff|rds)$", full.names = TRUE)
names(shipped) <- basename(shipped)

if (check) {
  drift <- character()
  for (name in names(shipped)) {
    target <- file.path(here, name)
    if (!file.exists(target)) { drift <- c(drift, paste0(name, " (missing)")); next }
    if (grepl("\\.rds$", name)) {
      if (!identical(readRDS(target), readRDS(shipped[[name]]))) drift <- c(drift, name)
    } else if (!identical(md5(target), md5(shipped[[name]]))) {
      drift <- c(drift, name)
    }
  }
  if (length(drift)) {
    stop("fixture drift detected: ", paste(drift, collapse = ", "))
  }
  msg("--check: ", length(shipped), " fixtures match")
} else {
  present <- list.files(here, pattern = "\\.(ssff|rds)$")
  if (length(present) && !force) {
    stop("fixtures already present in ", here,
         " -- re-run with --force to overwrite (or --check to verify)")
  }
  file.copy(shipped, here, overwrite = TRUE)
  writeLines(provenance, file.path(here, "PROVENANCE.md"))
  msg("installed ", length(shipped), " files into ", here)
}
