##' Track fundamental frequency using SRH via the voiceanalysis package
##'
##' Returns F0, voiced/unvoiced decisions, and SRH amplitude per frame using
##' the SRH algorithm \insertCite{Drugman2011SRH}{superassp}. Prefer this over
##' \code{\link{trk_pitch_srh}} when MATLAB-VAT parity matters.
##'
##' @details
##' Bit-faithful Rcpp port of the Summation of Residual Harmonics pitch tracker
##' from the original Kane MATLAB Voice Analysis Toolkit. The two implementations
##' differ in framing and smoothing details; use \code{trk_pitch_srh} for the
##' native superassp version.
##'
##' @inheritParams trk_pitch_rapt
##' @param minF Minimum F0 in Hz (default 50).
##' @param maxF Maximum F0 in Hz (default 500).
##'
##' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with tracks:
##'   \describe{
##'     \item{\code{f0}}{REAL32, fundamental frequency in Hz, n_frames × 1.
##'       Zero indicates unvoiced frames.}
##'     \item{\code{vad}}{REAL32, voiced/unvoiced decision (0 = unvoiced,
##'       1 = voiced), n_frames × 1.}
##'     \item{\code{srh_val}}{REAL32, SRH amplitude per frame, n_frames × 1.}
##'   }
##'   Frame rate: 100 Hz (fixed 10 ms hop). Audio is resampled to 16 kHz
##'   internally to match the MATLAB pipeline.
##'   If \code{toFile = TRUE}: invisibly returns the count of files written.
##'
##' @references
##' \insertCite{Drugman2011SRH}{superassp}
##' \insertCite{KaneGobl2013}{superassp}
##'
##' @seealso \code{\link{trk_pitch_srh}}
##' @examples
##' \dontrun{
##' trk_pitch_vat(
##'   system.file("samples", "sustained", "a1.wav", package = "superassp"),
##'   toFile = FALSE
##' )
##' }
##' @export
trk_pitch_vat <- function(listOfFiles,
                          beginTime = 0.0,
                          endTime = 0.0,
                          minF = 50,
                          maxF = 500,
                          toFile = TRUE,
                          explicitExt = "f0v",
                          outputDirectory = NULL,
                          verbose = TRUE) {

  if (FALSE) {
    cli::cli_abort(c(
      "Package {.pkg voiceanalysis} is required.",
      "i" = "Install with {.code pak::pkg_install('jckane/Voice_Analysis_Toolkit/voiceanalysis')}"
    ))
  }
  if (is.null(listOfFiles) || length(listOfFiles) == 0)
    cli::cli_abort("No input files specified in {.arg listOfFiles}")

  listOfFiles <- fast_strip_file_protocol(listOfFiles)
  listOfFiles <- normalizePath(path.expand(listOfFiles), mustWork = FALSE)
  if (!all(file.exists(listOfFiles)))
    cli::cli_abort(c("!" = "Some files do not exist:",
                     "x" = "{.file {fast_basename(listOfFiles[!file.exists(listOfFiles)])}}"))

  n_files <- length(listOfFiles)
  beginTime <- if (length(beginTime) == 1) rep(beginTime, n_files) else beginTime
  endTime   <- if (length(endTime)   == 1) rep(endTime,   n_files) else endTime

  makeOutputDirectory(outputDirectory, FALSE, "trk_pitch_vat")
  if (verbose) format_apply_msg("trk_pitch_vat", n_files, beginTime, endTime)

  results <- vector("list", n_files)
  target_sr <- 16000L

  for (i in seq_len(n_files)) {
    file_path <- listOfFiles[i]
    bt <- beginTime[i]; et <- endTime[i]
    tryCatch({
      invisible(utils::capture.output(
        audio_data <- av::read_audio_bin(
          audio = file_path,
          start_time = if (bt > 0) bt else NULL,
          end_time   = if (et > 0) et else NULL,
          channels = 1, sample_rate = target_sr
        ), type = "message"))

      vat_res <- .vat_pitch(as.numeric(audio_data),
                                          target_sr, minF, maxF,
                                          method = "srh")
      n_frames <- length(vat_res$f0)
      if (n_frames == 0) {
        cli::cli_warn("SRH returned empty result for {.file {basename(file_path)}}")
        results[[i]] <- if (toFile) FALSE else NULL; next
      }

      f0_v  <- ifelse(as.integer(vat_res$VUV) == 1L,
                       as.numeric(vat_res$f0), 0.0)
      out_obj <- list(
        f0      = matrix(f0_v, ncol = 1),
        vad     = matrix(as.numeric(vat_res$VUV), ncol = 1),
        srh_val = matrix(as.numeric(vat_res$SRHVal), ncol = 1)
      )
      attr(out_obj, "trackFormats") <- c("REAL32", "REAL32", "REAL32")
      attr(out_obj, "sampleRate")   <- 100  # 10 ms hop
      attr(out_obj, "origFreq")     <- as.numeric(target_sr)
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

attr(trk_pitch_vat, "ext")              <- "f0v"
attr(trk_pitch_vat, "tracks")           <- c("f0", "vad", "srh_val")
attr(trk_pitch_vat, "outputType")       <- "SSFF"
attr(trk_pitch_vat, "nativeFiletypes")  <- c("wav", "flac", "mp3", "mp4", "mkv", "avi")
attr(trk_pitch_vat, "suggestCaching")   <- FALSE
