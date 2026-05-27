##' Detect creaky voice via voiceanalysis (Kane-Drugman ANN)
##'
##' Bit-faithful Rcpp port of the Kane, Drugman & Gobl (2013) creak
##' detector \insertCite{KaneDrugmanGobl2013}{superassp}: 36-feature
##' Kane-Drugman + Ishi pipeline through a shallow ANN trained on the
##' original MATLAB VAT data. Returns a posterior probability per 10 ms
##' frame plus a binary decision (threshold 0.3 by default).
##'
##' Alternative to \code{\link{trk_covarep_creak}}, which uses a partial
##' R-side reimplementation. \code{trk_creak_vat} restores the full
##' 12-base × (static + delta + delta-delta) feature stack and the
##' logistic-output ANN that MATLAB's \code{patternnet} uses.
##'
##' @inheritParams trk_covarep_creak
##' @param threshold Numeric decision threshold for binarising the posterior
##'   (default 0.3, matches MATLAB).
##'
##' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with tracks:
##'   \describe{
##'     \item{\code{creak_pp}}{REAL32, posterior probability of creaky
##'       phonation, n_frames × 1.}
##'     \item{\code{creak_bin}}{REAL32, binary creak decision (0 / 1),
##'       n_frames × 1.}
##'   }
##'   Frame rate: 100 Hz (10 ms hop). Schema mirrors
##'   \code{trk_covarep_creak()} for swap-in convenience.
##'
##' @references
##' \insertCite{KaneDrugmanGobl2013}{superassp}
##' @seealso \code{\link{trk_covarep_creak}}
##' @examples
##' \dontrun{
##' trk_creak_vat(
##'   system.file("samples", "sustained", "a1.wav", package = "superassp"),
##'   toFile = FALSE
##' )
##' }
##' @export
trk_creak_vat <- function(listOfFiles,
                          beginTime = 0.0,
                          endTime = 0.0,
                          threshold = 0.3,
                          toFile = FALSE,
                          explicitExt = "crv",
                          outputDirectory = NULL,
                          verbose = TRUE) {

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
  makeOutputDirectory(outputDirectory, FALSE, "trk_creak_vat")
  if (verbose) format_apply_msg("trk_creak_vat", n_files, beginTime, endTime)

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

      cd <- .vat_creak_detect(wave, fs, threshold = threshold)
      n_frames <- length(cd$posterior)
      if (n_frames == 0) {
        cli::cli_warn("Creak detection returned no frames for {.file {basename(file_path)}}")
        results[[i]] <- if (toFile) FALSE else NULL; next
      }

      out_obj <- list(
        creak_pp  = matrix(as.numeric(cd$posterior), ncol = 1),
        creak_bin = matrix(as.numeric(cd$decision),  ncol = 1)
      )
      attr(out_obj, "trackFormats") <- c("REAL32", "REAL32")
      attr(out_obj, "sampleRate")   <- 100  # 10 ms hop
      attr(out_obj, "origFreq")     <- as.numeric(fs)
      attr(out_obj, "startTime")    <- as.numeric(bt)
      attr(out_obj, "startRecord")  <- 1L
      attr(out_obj, "endRecord")    <- as.integer(n_frames)
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

attr(trk_creak_vat, "ext")              <- "crv"
attr(trk_creak_vat, "tracks")           <- c("creak_pp", "creak_bin")
attr(trk_creak_vat, "outputType")       <- "SSFF"
attr(trk_creak_vat, "nativeFiletypes")  <- c("wav", "flac", "mp3", "mp4")
attr(trk_creak_vat, "suggestCaching")   <- FALSE
