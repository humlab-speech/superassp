##' Track fundamental frequency using SWIPE (Sawtooth Waveform Inspired Pitch Estimator)
##'
##' Extracts F0 by matching the speech spectrum against sawtooth waveform templates
##' across candidate F0 values via SPTK. SWIPE is particularly effective for noisy
##' or challenging recording conditions where cross-correlation trackers (RAPT, Snack)
##' are less reliable. Its default voicing threshold (0.3) is lower than RAPT's (0.6),
##' yielding more voiced frames.
##'
##' @inheritParams trk_pitch_rapt
##' @param voicing_threshold Numeric. Voicing decision threshold (0–1). Default 0.3
##'   (more permissive than RAPT). Increase toward 0.5 to reduce false voiced frames.
##'
##' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with track:
##'   \describe{
##'     \item{\code{f0}}{REAL32, fundamental frequency in Hz, n_frames × 1.
##'       Zero indicates unvoiced frames.}
##'   }
##'   Frame rate: \code{1000 / windowShift} Hz (default 100 Hz).
##'   If \code{toFile = TRUE}: integer count of files written, returned invisibly.
##'
##' @export
##' @examples
##' \dontrun{
##' # Extract F0 from audio file
##' trk_pitch_swipe("recording.wav")
##'
##' # Process with custom parameters
##' trk_pitch_swipe("speech.wav", minF = 100, maxF = 500, voicing_threshold = 0.4)
##' }
trk_pitch_swipe <- function(listOfFiles,
                  beginTime = 0.0,
                  endTime = 0.0,
                  windowShift = 10.0,
                  minF = 60.0,
                  maxF = 400.0,
                  voicing_threshold = 0.3,
                  toFile = TRUE,
                  explicitExt = "f0",
                  outputDirectory = NULL,
                  verbose = TRUE) {

  if (is.null(listOfFiles) || length(listOfFiles) == 0) {
    cli::cli_abort("No input files specified in {.arg listOfFiles}")
  }

  listOfFiles <- fast_strip_file_protocol(listOfFiles)
  listOfFiles <- normalizePath(path.expand(listOfFiles), mustWork = FALSE)

  files_exist <- file.exists(listOfFiles)
  if (!all(files_exist)) {
    missing_files <- listOfFiles[!files_exist]
    cli::cli_abort(c(
      "!" = "Some files do not exist:",
      "x" = "{.file {fast_basename(missing_files)}}"
    ))
  }

  n_files <- length(listOfFiles)

  beginTime <- if (is.null(beginTime)) 0.0 else beginTime
  endTime <- if (is.null(endTime)) 0.0 else endTime

  if (length(beginTime) == 1) beginTime <- rep(beginTime, n_files)
  if (length(endTime) == 1) endTime <- rep(endTime, n_files)

  makeOutputDirectory(outputDirectory, FALSE, "trk_pitch_swipe")

  if (verbose) format_apply_msg("trk_pitch_swipe", n_files, beginTime, endTime)

  results <- vector("list", n_files)

  if (verbose && n_files > 1) {
    cli::cli_progress_bar(
      "Processing files",
      total = n_files,
      format = "{cli::pb_spin} {cli::pb_current}/{cli::pb_total} | ETA: {cli::pb_eta}"
    )
  }

  for (i in seq_len(n_files)) {
    file_path <- listOfFiles[i]
    bt <- beginTime[i]
    et <- endTime[i]

    tryCatch({
      audio_obj <- read_audio(
        file_path,
        begin = bt,
        end   = et
      )

      swipe_result <- swipe_cpp(
        audio_obj = audio_obj,
        minF = minF,
        maxF = maxF,
        windowShift = windowShift,
        voicing_threshold = voicing_threshold,
        verbose = FALSE
      )

      out_obj <- create_f0_asspobj(swipe_result, windowShift)

      if (toFile) {
        out_file <- generate_output_path(file_path, explicitExt, outputDirectory)
        write.AsspDataObj(out_obj, out_file)
        results[[i]] <- TRUE
      } else {
        results[[i]] <- out_obj
      }

    }, error = function(e) {
      cli::cli_warn("Error processing {.file {basename(file_path)}}: {conditionMessage(e)}")
      results[[i]] <- if (toFile) FALSE else NULL
    })

    if (verbose && n_files > 1) {
      cli::cli_progress_update()
    }
  }

  if (verbose && n_files > 1) {
    cli::cli_progress_done()
  }

  if (toFile) {
    n_success <- sum(unlist(results), na.rm = TRUE)
    if (verbose) {
      cli::cli_inform("Successfully processed {n_success} of {n_files} file{?s}")
    }
    return(invisible(n_success))
  } else {
    results <- results[!sapply(results, is.null)]
    if (length(results) == 1) {
      return(results[[1]])
    } else {
      return(results)
    }
  }
}

attr(trk_pitch_swipe, "ext") <- "f0"
attr(trk_pitch_swipe, "tracks") <- c("f0")
attr(trk_pitch_swipe, "outputType") <- "SSFF"
attr(trk_pitch_swipe, "nativeFiletypes") <- c("wav", "flac", "mp3", "mp4", "mkv", "avi")
attr(trk_pitch_swipe, "suggestCaching") <- FALSE
