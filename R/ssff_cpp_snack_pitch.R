##' Track fundamental frequency using the Snack/ESPS dp_f0 algorithm
##'
##' Extracts F0, voicing probability, RMS energy, and autocorrelation peak using
##' the Snack Sound Toolkit normalized cross-correlation + dynamic-programming pitch
##' tracker (\code{dp_f0}, Talkin 1995). Unlike \code{\link{trk_pitch_rapt}}, which
##' returns only F0, this function exposes all four tracks for downstream signal
##' quality assessment.
##'
##' @param listOfFiles Character vector of audio file paths. Any format supported by
##'   \pkg{av} is accepted; non-native inputs are transcoded automatically.
##' @param beginTime Numeric. Start of analysis window in seconds. Default 0 (file start).
##' @param endTime Numeric. End of analysis window in seconds. Default 0 (file end).
##' @param windowShift Numeric. Frame shift in milliseconds; sets output frame rate
##'   (1000 / windowShift Hz). Default 10.0 ms.
##' @param minF Numeric. Minimum F0 in Hz. Default 50.0 Hz.
##' @param maxF Numeric. Maximum F0 in Hz. Default 550.0 Hz.
##' @param voiceBias Numeric. Bias toward the voiced hypothesis in the DP cost function
##'   (range approximately −0.5 to 0.5; positive = more voiced frames). Default 0.0.
##' @param toFile Logical. If \code{TRUE}, write SSFF output files and return the
##'   count written invisibly. If \code{FALSE}, return an \code{AsspDataObj}.
##'   Default \code{TRUE}.
##' @param explicitExt Character. Output file extension. Default \code{"snackpitch"}.
##' @param outputDirectory Character. Directory for output files. \code{NULL} (default)
##'   writes alongside the input file.
##' @param verbose Logical. Print per-file progress. Default \code{TRUE}.
##'
##' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with tracks:
##'   \describe{
##'     \item{\code{f0}}{REAL32, fundamental frequency in Hz, n_frames × 1.
##'       Zero indicates unvoiced frames.}
##'     \item{\code{voicing}}{REAL32, voicing probability, 0–1, n_frames × 1.}
##'     \item{\code{rms}}{REAL32, RMS energy (linear), dimensionless, n_frames × 1.}
##'     \item{\code{acpeak}}{REAL32, autocorrelation peak magnitude, 0–1, n_frames × 1.}
##'   }
##'   Frame rate: \code{1000 / windowShift} Hz (default 100 Hz).
##'   If \code{toFile = TRUE}: integer count of files written, returned invisibly.
##'
##' @export
##' @examples
##' \dontrun{
##' # Full 4-track pitch analysis
##' res <- trk_pitch_snack("recording.wav", toFile = FALSE)
##' names(res)  # "f0" "voicing" "rms" "acpeak"
##'
##' # Custom F0 range
##' trk_pitch_snack("speech.mp3", minF = 75, maxF = 300)
##' }
trk_pitch_snack <- function(listOfFiles,
                       beginTime = 0.0,
                       endTime = 0.0,
                       windowShift = 10.0,
                       minF = 50.0,
                       maxF = 550.0,
                       voiceBias = 0.0,
                       toFile = TRUE,
                       explicitExt = "snackpitch",
                       outputDirectory = NULL,
                       verbose = TRUE) {

  if (is.null(listOfFiles) || length(listOfFiles) == 0)
    cli::cli_abort("No input files specified in {.arg listOfFiles}")

  listOfFiles <- fast_strip_file_protocol(listOfFiles)
  listOfFiles <- normalizePath(path.expand(listOfFiles), mustWork = FALSE)

  files_exist <- file.exists(listOfFiles)
  if (!all(files_exist)) {
    missing <- listOfFiles[!files_exist]
    cli::cli_abort(c("!" = "Files not found:", "x" = "{.file {fast_basename(missing)}}"))
  }

  n_files <- length(listOfFiles)
  beginTime <- if (is.null(beginTime)) 0.0 else beginTime
  endTime   <- if (is.null(endTime))   0.0 else endTime
  if (length(beginTime) == 1) beginTime <- rep(beginTime, n_files)
  if (length(endTime)   == 1) endTime   <- rep(endTime,   n_files)

  makeOutputDirectory(outputDirectory, FALSE, "trk_pitch_snack")
  if (verbose) format_apply_msg("trk_pitch_snack", n_files, beginTime, endTime)

  results <- vector("list", n_files)

  if (verbose && n_files > 1)
    cli::cli_progress_bar("Processing files", total = n_files,
      format = "{cli::pb_spin} {cli::pb_current}/{cli::pb_total} | ETA: {cli::pb_eta}")

  for (i in seq_len(n_files)) {
    file_path <- listOfFiles[i]
    bt <- beginTime[i]; et <- endTime[i]

    tryCatch({
      audio_obj <- read_audio(file_path, begin = bt, end = et)

      res <- snackp_cpp(
        audio_obj   = audio_obj,
        minF        = minF,
        maxF        = maxF,
        windowShift = windowShift,
        voiceBias   = voiceBias,
        verbose     = FALSE
      )

      out_obj <- create_snackp_asspobj(res, windowShift)

      if (toFile) {
        out_file <- generate_output_path(file_path, explicitExt, outputDirectory)
        write.AsspDataObj(out_obj, out_file)
        results[[i]] <- TRUE
      } else {
        results[[i]] <- out_obj
      }
    }, error = function(e) {
      cli::cli_warn("Error processing {.file {basename(file_path)}}: {conditionMessage(e)}")
      results[[i]] <<- if (toFile) FALSE else NULL
    })

    if (verbose && n_files > 1) cli::cli_progress_update()
  }

  if (verbose && n_files > 1) cli::cli_progress_done()

  if (toFile) {
    n_ok <- sum(unlist(results), na.rm = TRUE)
    if (verbose) cli::cli_inform("Processed {n_ok} of {n_files} file{?s}")
    return(invisible(n_ok))
  } else {
    results <- results[!sapply(results, is.null)]
    if (length(results) == 1) results[[1]] else results
  }
}

attr(trk_pitch_snack, "ext")             <- "snackpitch"
attr(trk_pitch_snack, "tracks")          <- c("f0", "voicing", "rms", "acpeak")
attr(trk_pitch_snack, "outputType")      <- "SSFF"
attr(trk_pitch_snack, "nativeFiletypes") <- c("wav", "flac", "mp3", "mp4", "mkv", "avi")
