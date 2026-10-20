# PSOCK cluster reused across calls (Windows has no fork). Starting R workers
# and loading superassp in each costs ~1-3 s, which used to be paid on every
# batch call. The cluster is kept until a different size is requested, a
# worker stops answering, or the package unloads (.onUnload).
.cluster_cache <- new.env(parent = emptyenv())

.superassp_cluster <- function(n_cores) {
  cl <- .cluster_cache$cl
  if (!is.null(cl) && length(cl) == n_cores &&
      isTRUE(tryCatch(all(unlist(parallel::clusterEvalQ(cl, TRUE))),
                      error = function(e) FALSE))) {
    return(cl)
  }
  .stop_superassp_cluster()
  cl <- parallel::makeCluster(n_cores)
  parallel::clusterEvalQ(cl, library(superassp))
  .cluster_cache$cl <- cl
  cl
}

.stop_superassp_cluster <- function() {
  if (!is.null(.cluster_cache$cl)) {
    try(parallel::stopCluster(.cluster_cache$cl), silent = TRUE)
    .cluster_cache$cl <- NULL
  }
}

#' Run a per-file closure sequentially or in parallel
#'
#' Shared dispatch logic used by batch DSP wrappers: auto-enables parallel
#' processing for 2+ files, picks a fork cluster (Unix) or PSOCK cluster
#' (Windows), and falls back to a plain sequential loop (with an optional
#' progress bar) otherwise. Extracted from the working implementation in
#' `trk_pitchmark_estk()` so every batch wrapper shares one dispatch path.
#'
#' @param n_files Integer number of files to process.
#' @param process_single_file Function of one argument (integer index `i`)
#'   that processes file `i` and returns its result. Must be self-contained
#'   (read all inputs from its enclosing environment; write nothing to it).
#' @param parallel Logical or `NULL`. `NULL` (default) auto-enables for
#'   `n_files > 1`.
#' @param n_cores Integer or `NULL`. `NULL` (default) uses
#'   `parallel::detectCores() - 1` (minimum 1).
#' @param verbose Logical. Show a progress bar (sequential path) or a
#'   progress-aware parallel apply (`pbapply`/`pbmcapply`, if installed).
#' @param export_vars Character vector of variable names `process_single_file`
#'   references from its enclosing environment. Required for the Windows
#'   PSOCK path (`parallel::clusterExport()`); ignored on the fork path.
#' @param export_env Environment to export `export_vars` from. Default the
#'   caller's environment.
#' @param progress_label Character. Label for the sequential-path progress bar.
#'
#' @return A `list` of length `n_files`, one element per file, in input order.
#' @keywords internal
#' @noRd
run_parallel_files <- function(n_files, process_single_file,
                                parallel = NULL, n_cores = NULL,
                                verbose = TRUE, export_vars = character(0),
                                export_env = parent.frame(),
                                progress_label = "Processing files") {

  if (is.null(parallel)) parallel <- n_files > 1

  if (is.null(n_cores)) {
    n_cores <- parallel::detectCores() - 1
    if (is.na(n_cores) || n_cores < 1) n_cores <- 1
  }
  # CRAN/R CMD check caps mclapply() to 2 cores (parallel:::.check_ncores())
  if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))) {
    n_cores <- min(n_cores, 2L)
  }

  use_parallel <- parallel && n_files > 1 && n_cores > 1

  if (use_parallel) {
    if (verbose) {
      cli::cli_inform("Using parallel processing on {n_cores} core{?s}")
    }
    if (.Platform$OS.type == "windows") {
      cl <- .superassp_cluster(n_cores)
      parallel::clusterExport(cl, export_vars, envir = export_env)

      if (verbose && requireNamespace("pbapply", quietly = TRUE)) {
        results <- pbapply::pblapply(seq_len(n_files), process_single_file, cl = cl)
      } else {
        results <- parallel::parLapply(cl, seq_len(n_files), process_single_file)
      }
    } else {
      if (verbose && requireNamespace("pbmcapply", quietly = TRUE)) {
        results <- pbmcapply::pbmclapply(
          seq_len(n_files),
          process_single_file,
          mc.cores = n_cores,
          mc.preschedule = TRUE
        )
      } else {
        results <- parallel::mclapply(
          seq_len(n_files),
          process_single_file,
          mc.cores = n_cores,
          mc.preschedule = TRUE
        )
      }
    }

    failed <- vapply(results, inherits, logical(1), "try-error")
    if (any(failed)) {
      cli::cli_warn("{sum(failed)} worker{?s} failed unexpectedly (e.g. crashed or was killed); treating as missing result{?s}")
      results[failed] <- list(NULL)
    }
  } else {
    results <- vector("list", n_files)
    if (verbose && n_files > 1) {
      cli::cli_progress_bar(
        progress_label,
        total = n_files,
        format = "{cli::pb_spin} {cli::pb_current}/{cli::pb_total} | ETA: {cli::pb_eta}"
      )
      for (i in seq_len(n_files)) {
        results[[i]] <- process_single_file(i)
        cli::cli_progress_update()
      }
      cli::cli_progress_done()
    } else {
      for (i in seq_len(n_files)) {
        results[[i]] <- process_single_file(i)
      }
    }
  }

  results
}
