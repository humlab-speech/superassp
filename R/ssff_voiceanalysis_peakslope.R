##' Peak slope via voiceanalysis Daless wavelet bank
##'
##' Computes the peak slope acoustic parameter \insertCite{KaneGobl2011}{superassp}
##' using the bit-faithful Daless wavelet bank in \pkg{voiceanalysis}. Faster
##' and more numerically stable than the pure-R \code{\link{trk_peakslope}}
##' (which approximates Daless via db4 wavelets) for the same algorithm.
##'
##' @inheritParams trk_covarep_creak
##'
##' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with one track:
##'   \describe{
##'     \item{\code{peak_slope}}{REAL32, peak-slope coefficient per frame,
##'       n_frames × 1.}
##'   }
##'   Frame rate: 100 Hz (10 ms hop).
##'
##' @references
##' \insertCite{KaneGobl2011}{superassp}
##' @seealso \code{\link{trk_peakslope}}
##' @export
trk_peakslope_vat <- function(listOfFiles,
                              beginTime = 0.0,
                              endTime = 0.0,
                              toFile = FALSE,
                              explicitExt = "psv",
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
  makeOutputDirectory(outputDirectory, FALSE, "trk_peakslope_vat")
  if (verbose) format_apply_msg("trk_peakslope_vat", n_files, beginTime, endTime)

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

      ps <- .vat_peak_slope(wave, fs)
      n_frames <- length(ps)
      if (n_frames == 0) {
        cli::cli_warn("peak slope returned 0 frames for {.file {basename(file_path)}}")
        results[[i]] <- if (toFile) FALSE else NULL; next
      }
      out_obj <- list(peak_slope = matrix(as.numeric(ps), ncol = 1))
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

attr(trk_peakslope_vat, "ext")              <- "psv"
attr(trk_peakslope_vat, "tracks")           <- c("peak_slope")
attr(trk_peakslope_vat, "outputType")       <- "SSFF"
attr(trk_peakslope_vat, "nativeFiletypes")  <- c("wav", "flac", "mp3", "mp4")
attr(trk_peakslope_vat, "suggestCaching")   <- FALSE
