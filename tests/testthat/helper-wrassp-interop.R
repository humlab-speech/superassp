# Runs the wrassp side of the cross-package SSFF checks in a subprocess.
#
# Why a subprocess: loading wrassp into the test session registers its
# as_tibble()/print() methods for AsspDataObj on top of superassp's (R prints
# "Registered S3 methods overwritten by 'wrassp'"), which would make the rest of
# the session's dispatch depend on test order. The subprocess also keeps the
# checked-against package at arm's length from the code under test.
#
# Returns NULL when wrassp is not installed (the caller skips) and errors when
# the subprocess fails while wrassp is available (a real problem, not a skip).
#
# NOTE: availability is checked with system.file(), never with
# requireNamespace()/is_installed(). Loading wrassp's namespace -- even without
# attaching it -- registers wrassp's as_tibble.AsspDataObj()/print.AsspDataObj()
# over superassp's for the rest of the session, which silently breaks
# test-track-naming-phase2.R and any other test that consumes the tibble view.

wrassp_installed <- function() nzchar(system.file(package = "wrassp"))

wrassp_interop_script <- function() {
  testthat::test_path("scripts", "wrassp_interop.R")
}

wrassp_interop <- function() {
  if (!wrassp_installed()) return(NULL)
  script <- wrassp_interop_script()
  if (!file.exists(script)) return(NULL)

  wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  if (!nzchar(wav) || !file.exists(wav)) return(NULL)

  out <- file.path(tempdir(), "wrassp-interop")
  unlink(out, recursive = TRUE)
  dir.create(out, recursive = TRUE)

  # The library path goes through the environment, not system2(env = ): on
  # Windows system2() pastes `env` after the command as plain arguments, so
  # Rscript received "R_LIBS=..." as its script and died without output.
  rscript <- file.path(R.home("bin"), "Rscript")
  log <- withr::with_envvar(
    c(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep)),
    suppressWarnings(system2(
      rscript,
      args = c(shQuote(script), shQuote(wav), shQuote(out)),
      stdout = TRUE, stderr = TRUE
    ))
  )

  rds <- file.path(out, "interop.rds")
  if (!file.exists(rds)) {
    status <- attr(log, "status")
    stop("wrassp interop subprocess failed (status ",
         if (is.null(status)) 0L else status, "):\n",
         paste(log, collapse = "\n"))
  }
  readRDS(rds)
}
