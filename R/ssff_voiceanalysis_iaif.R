##' Estimate glottal flow via IAIF using voiceanalysis
##'
##' Iterative Adaptive Inverse Filtering \insertCite{Alku1992}{superassp},
##' driven by the bit-faithful Rcpp port in the \pkg{voiceanalysis} package.
##' GCIs are detected internally using SE-VQ.
##'
##' Alternative to \code{\link{trk_covarep_iaif}}, which uses superassp's
##' own native C++ IAIF implementation. The two implementations differ in
##' parameterisation (\code{trk_iaif_vat} uses the MATLAB-VAT defaults and
##' the SE-VQ GCI grid; \code{trk_covarep_iaif} uses COVAREP defaults and a
##' contiguous frame grid).
##'
##' @inheritParams trk_covarep_iaif
##' @param p LPC prediction order. \code{NULL} (default) sets \code{round(fs/1000)+2}.
##'
##' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with tracks:
##'   \describe{
##'     \item{\code{glottal_flow}}{REAL64, integrated glottal flow,
##'       n_samples × 1.}
##'     \item{\code{glottal_derivative}}{REAL64, glottal flow derivative,
##'       n_samples × 1.}
##'   }
##'   Frame rate equals the audio sample rate. Schema matches
##'   \code{trk_covarep_iaif()} for swap-in convenience.
##'   If \code{toFile = TRUE}: invisibly returns the count of files written.
##'
##' @references
##' \insertCite{Alku1992}{superassp}
##' \insertCite{KaneGobl2013}{superassp}
##' @seealso \code{\link{trk_covarep_iaif}}, \code{\link{trk_gci_vat}}
##' @examples
##' \dontrun{
##' trk_iaif_vat(
##'   system.file("samples", "sustained", "a1.wav", package = "superassp"),
##'   toFile = FALSE
##' )
##' }
##' @export
trk_iaif_vat <- function(listOfFiles,
                          beginTime = 0.0,
                          endTime = 0.0,
                          p = NULL,
                          toFile = TRUE,
                          explicitExt = "glv",
                          outputDirectory = NULL,
                          verbose = TRUE,
                          ...) {

  if (FALSE)
    cli::cli_abort(c("Package {.pkg voiceanalysis} is required.",
                     "i" = "Install via {.code pak::pkg_install('jckane/Voice_Analysis_Toolkit/voiceanalysis')}"))

  if (is.null(listOfFiles) || length(listOfFiles) == 0)
    cli::cli_abort("No input files specified in {.arg listOfFiles}")
  listOfFiles <- normalizePath(path.expand(fast_strip_file_protocol(listOfFiles)), mustWork = FALSE)
  if (!all(file.exists(listOfFiles)))
    cli::cli_abort("Some files do not exist.")

  n_files <- length(listOfFiles)
  beginTime <- if (length(beginTime) == 1) rep(beginTime, n_files) else beginTime
  endTime   <- if (length(endTime)   == 1) rep(endTime,   n_files) else endTime
  makeOutputDirectory(outputDirectory, FALSE, "trk_iaif_vat")
  if (verbose) format_apply_msg("trk_iaif_vat", n_files, beginTime, endTime)

  results <- vector("list", n_files)
  for (i in seq_len(n_files)) {
    file_path <- listOfFiles[i]
    bt <- beginTime[i]; et <- endTime[i]
    tryCatch({
      invisible(utils::capture.output(
        audio_data <- av::read_audio_bin(
          audio = file_path,
          start_time = if (bt > 0) bt else NULL,
          end_time   = if (et > 0) et else NULL,
          channels = 1
        ), type = "message"))
      fs <- attr(audio_data, "sample_rate")
      wave <- as.numeric(audio_data)
      mx <- max(abs(wave)); if (mx > 1) wave <- wave / mx

      gci_res <- .vat_se_vq(wave, fs)
      gci <- as.integer(gci_res$GCI)
      if (length(gci) < 3) {
        cli::cli_warn("Too few GCIs ({length(gci)}) for IAIF on {.file {basename(file_path)}}")
        results[[i]] <- if (toFile) FALSE else NULL; next
      }
      res <- .vat_iaif(wave, fs, GCI = gci, p = p)
      n_samples <- length(res$g)
      if (n_samples == 0) {
        cli::cli_warn("IAIF failed for {.file {basename(file_path)}}")
        results[[i]] <- if (toFile) FALSE else NULL; next
      }

      out_obj <- list(
        glottal_flow       = matrix(as.numeric(res$g),  ncol = 1),
        glottal_derivative = matrix(as.numeric(res$dg), ncol = 1)
      )
      attr(out_obj, "trackFormats") <- c("REAL64", "REAL64")
      attr(out_obj, "sampleRate")   <- as.numeric(fs)
      attr(out_obj, "origFreq")     <- as.numeric(fs)
      attr(out_obj, "startTime")    <- as.numeric(bt)
      attr(out_obj, "startRecord")  <- 1L
      attr(out_obj, "endRecord")    <- as.integer(n_samples)
      attr(out_obj, "fileInfo")     <- c(20L, 2L)
      class(out_obj) <- "AsspDataObj"

      if (toFile) {
        out_file <- generate_output_path(file_path, explicitExt, outputDirectory)
        write.AsspDataObj(out_obj, out_file)
        results[[i]] <- TRUE
      } else {
        results[[i]] <- out_obj
      }
    }, error = function(e) {
      warning("Error processing ", basename(file_path), ": ", e$message, call. = FALSE)
      results[i] <- list(NULL)
    })
  }

  if (toFile) invisible(sum(unlist(results), na.rm = TRUE))
  else if (n_files == 1) results[[1]] else results
}

attr(trk_iaif_vat, "ext")              <- "glv"
attr(trk_iaif_vat, "tracks")           <- c("glottal_flow", "glottal_derivative")
attr(trk_iaif_vat, "outputType")       <- "SSFF"
attr(trk_iaif_vat, "nativeFiletypes")  <- c("wav", "flac", "mp3", "mp4")
attr(trk_iaif_vat, "suggestCaching")   <- FALSE
