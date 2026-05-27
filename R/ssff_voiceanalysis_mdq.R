##' Maxima Dispersion Quotient (MDQ) via voiceanalysis
##'
##' Computes the Maxima Dispersion Quotient
##' \insertCite{KaneGobl2013MDQ}{superassp} for breathy / tense voice
##' discrimination. Uses the bit-faithful Daless wavelet bank + per-GCI
##' dispersion measurement from \pkg{voiceanalysis}.
##'
##' MDQ is computed per-GCI, then resampled to a fixed 100 Hz frame grid
##' (10 ms hop) so it can sit alongside other \code{trk_*} tracks.
##'
##' @inheritParams trk_covarep_creak
##'
##' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with one track:
##'   \describe{
##'     \item{\code{mdq}}{REAL32, MDQ value, n_frames × 1. Higher values
##'       indicate more dispersed wavelet maxima (breathier voice).}
##'   }
##'   Frame rate: 100 Hz (10 ms hop).
##'
##' @references
##' \insertCite{KaneGobl2013MDQ}{superassp}
##' @examples
##' \dontrun{
##' trk_mdq_vat(
##'   system.file("samples", "sustained", "a1.wav", package = "superassp"),
##'   toFile = FALSE
##' )
##' }
##' @export
trk_mdq_vat <- function(listOfFiles,
                        beginTime = 0.0,
                        endTime = 0.0,
                        toFile = FALSE,
                        explicitExt = "mdq",
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
  makeOutputDirectory(outputDirectory, FALSE, "trk_mdq_vat")
  if (verbose) format_apply_msg("trk_mdq_vat", n_files, beginTime, endTime)

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

      se <- .vat_se_vq(wave, fs)
      if (length(se$GCI) < 3) {
        cli::cli_warn("Too few GCIs for MDQ on {.file {basename(file_path)}}")
        results[[i]] <- if (toFile) FALSE else NULL; next
      }
      mdq_per_gci <- .vat_mdq(se$res, fs, se$GCI)

      # Resample per-GCI -> 10 ms grid via stats::approx using GCI times
      n_frames <- max(1L, floor(length(wave) / round(0.010 * fs)))
      gci_t <- (se$GCI - 1) / fs  # seconds
      target_t <- seq(0, by = 0.010, length.out = n_frames)
      mdq_t <- stats::approx(gci_t, as.numeric(mdq_per_gci),
                              xout = target_t, rule = 2)$y
      mdq_t[is.na(mdq_t)] <- 0

      out_obj <- list(mdq = matrix(mdq_t, ncol = 1))
      attr(out_obj, "trackFormats") <- c("REAL32")
      attr(out_obj, "sampleRate")   <- 100
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

attr(trk_mdq_vat, "ext")              <- "mdq"
attr(trk_mdq_vat, "tracks")           <- c("mdq")
attr(trk_mdq_vat, "outputType")       <- "SSFF"
attr(trk_mdq_vat, "nativeFiletypes")  <- c("wav", "flac", "mp3", "mp4")
attr(trk_mdq_vat, "suggestCaching")   <- FALSE
