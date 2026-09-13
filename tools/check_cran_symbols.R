#!/usr/bin/env Rscript
#
# check_cran_symbols.R -- gate for CRAN's "checking compiled code" scan.
#
# R's check scans <pkg>/libs/ for undefined symbols that are known to terminate R
# or to write outside R's console (the table is tools:::so_symbol_names_table).
# This script runs the same scan and exits non-zero on any hit, so the fix set
# cannot silently regress.
#
# Usage:  Rscript tools/check_cran_symbols.R [installed-package-dir]
#         (default: the installed superassp)

args <- commandArgs(trailingOnly = TRUE)
pkgdir <- if (length(args)) args[[1L]] else find.package("superassp")

if (!dir.exists(file.path(pkgdir, "libs"))) {
  stop("no libs/ directory in ", pkgdir,
       " -- check an installed package, not a source tree", call. = FALSE)
}

res <- tools:::check_compiled_code(pkgdir)
txt <- format(res)
cat(txt, sep = "\n")

# check_compiled_code() returns an empty list when clean, otherwise one entry per
# shared object with the offending symbols.
bad <- length(res) > 0L || any(nzchar(txt))
if (bad) {
  cat("\nFAIL: banned entry points in compiled code (see above)\n")
  quit(status = 1L)
}
cat("\nOK: no banned entry points in compiled code\n")
quit(status = 0L)
