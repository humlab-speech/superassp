# Contract-attribute audit for trk_* functions.
#
# The CLAUDE.md project spec requires every user-facing trk_* function to:
#   1. Carry attributes ext, tracks, outputType, nativeFiletypes.
#   2. Default toFile = FALSE in its signature.
#
# This file runs both checks against the live namespace.

required <- c("ext", "tracks", "outputType", "nativeFiletypes")

collect_attr_problems <- function() {
  ns      <- asNamespace("superassp")
  exports <- getNamespaceExports("superassp")
  trk     <- exports[grepl("^trk_", exports)]
  problems <- character(0)
  for (nm in trk) {
    fn   <- get(nm, envir = ns)
    have <- names(attributes(fn))
    miss <- setdiff(required, have)
    if (length(miss) > 0) {
      problems <- c(problems, sprintf("%s missing: %s", nm,
                                      paste(miss, collapse = ", ")))
    }
  }
  problems
}

collect_toFile_problems <- function() {
  ns      <- asNamespace("superassp")
  exports <- getNamespaceExports("superassp")
  trk     <- exports[grepl("^trk_", exports)]
  problems <- character(0)
  for (nm in trk) {
    fn   <- get(nm, envir = ns)
    fmls <- formals(fn)
    if (!"toFile" %in% names(fmls)) next
    if (!identical(fmls$toFile, FALSE)) {
      problems <- c(problems, sprintf("%s: toFile default = %s",
                                      nm, deparse(fmls$toFile)))
    }
  }
  problems
}

test_that("every exported trk_* carries the required contract attributes", {
  problems <- collect_attr_problems()
  expect_equal(problems, character(0),
               info = paste(problems, collapse = "\n"))
})

test_that("every exported trk_* defaults toFile=FALSE", {
  problems <- collect_toFile_problems()
  expect_equal(problems, character(0),
               info = paste(problems, collapse = "\n"))
})
