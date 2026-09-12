# Shared guard for tests that exercise pladdrr-backed functions.
#
# pladdrr is an optional, GitHub-only dependency (see DESCRIPTION): it is not
# installed on CRAN or win-builder check machines, so every test that calls an
# lst_*/trk_* function backed by it must skip rather than fail. Use this
# helper instead of writing the skip condition inline, so that a missing
# pladdrr always produces the same SKIP reason.
skip_without_pladdrr <- function() {
  testthat::skip_if_not(
    superassp:::pladdrr_available(),
    "pladdrr not available"
  )
}
