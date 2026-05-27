##' Estimate band aperiodicity using the D4C algorithm (WORLD vocoder)
##'
##' Computes per-frame, per-band aperiodicity using the D4C (Death, Destruction,
##' Diversion and Disgrace) estimator from the WORLD vocoder via SPTK. Aperiodicity
##' quantifies noise-to-harmonic energy ratio across frequency bands and is the
##' primary input to WORLD's noise excitation model. Use alongside
##' \code{trk_pitch_rapt()} or \code{trk_pitch_swipe()} for full WORLD vocoder
##' analysis/resynthesis.
##'
##' @param listOfFiles Character vector of audio file paths. Any format supported by
##'   \pkg{av} is accepted; non-native inputs are transcoded automatically.
##' @param beginTime Numeric. Start of analysis window in seconds. Default 0 (file start).
##' @param endTime Numeric. End of analysis window in seconds. Default 0 (file end).
##' @param windowShift Numeric. Frame shift in milliseconds; sets output frame rate
##'   (1000 / windowShift Hz). Default 5.0 ms.
##' @param minF Numeric. Minimum F0 in Hz used for the internal pitch estimator.
##'   Default 60.0 Hz.
##' @param maxF Numeric. Maximum F0 in Hz used for the internal pitch estimator.
##'   Default 400.0 Hz.
##' @param voicing_threshold Numeric. Voicing threshold for the internal F0 detector
##'   (0–1; higher = more conservative). Default 0.85.
##' @param threshold Numeric. D4C aperiodicity clipping threshold (0–1). Default 0.85.
##' @param toFile Logical. If \code{TRUE}, write SSFF output files and return the
##'   count written invisibly. If \code{FALSE}, return an \code{AsspDataObj}.
##'   Default \code{TRUE}.
##' @param explicitExt Character. Output file extension. Default \code{"ap"}.
##' @param outputDirectory Character. Directory for output files. \code{NULL} (default)
##'   writes alongside the input file.
##' @param verbose Logical. Print per-file progress. Default \code{TRUE}.
##'
##' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with track:
##'   \describe{
##'     \item{\code{aperiodicity}}{REAL32, band aperiodicity, 0–1 (0 = fully periodic,
##'       1 = fully aperiodic), n_frames × n_bands where n_bands =
##'       floor(fs/2 / 3000) + 1.}
##'   }
##'   Frame rate: \code{1000 / windowShift} Hz (default 200 Hz).
##'   If \code{toFile = TRUE}: integer count of files written, returned invisibly.
##'
##' @export
##' @examples
##' \dontrun{
##' # Estimate aperiodicity
##' trk_d4c("recording.wav")
##'
##' # Process with custom parameters
##' trk_d4c("speech.wav", minF = 80, maxF = 350, windowShift = 10)
##'
##' # Process multiple files
##' trk_d4c(c("file1.wav", "file2.wav"))
##' }
trk_d4c <- function(listOfFiles,
                beginTime = 0.0,
                endTime = 0.0,
                windowShift = 5.0,
                minF = 60.0,
                maxF = 400.0,
                voicing_threshold = 0.1,
                threshold = 0.85,
                toFile = TRUE,
                explicitExt = "ap",
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

  makeOutputDirectory(outputDirectory, FALSE, "trk_d4c")

  if (verbose) format_apply_msg("trk_d4c", n_files, beginTime, endTime)

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

      harvest_result <- harvest_cpp(
        audio_obj = audio_obj,
        minF = minF,
        maxF = maxF,
        windowShift = windowShift,
        voicing_threshold = voicing_threshold,
        verbose = FALSE
      )

      d4c_result <- d4c_cpp(
        audio_obj = audio_obj,
        f0 = as.numeric(harvest_result$f0),
        temporal_positions = as.numeric(harvest_result$times),
        threshold = threshold,
        verbose = FALSE
      )

      out_obj <- create_aperiodicity_asspobj(d4c_result, windowShift)

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

attr(trk_d4c, "ext") <- "ap"
attr(trk_d4c, "tracks") <- c("aperiodicity")
attr(trk_d4c, "outputType") <- "SSFF"
attr(trk_d4c, "nativeFiletypes") <- c("wav", "flac", "mp3", "mp4", "mkv", "avi")
attr(trk_d4c, "suggestCaching") <- FALSE
