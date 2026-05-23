##' CheapTrick Spectral Envelope Estimation (WORLD vocoder, C++ implementation)
##'
##' @description Estimate the spectral envelope using the CheapTrick algorithm
##'   from the WORLD vocoder (Morise 2015). CheapTrick uses a Hanning-windowed,
##'   F0-adaptive spectral analysis to produce a smooth, high-quality power
##'   spectral envelope per frame.
##'
##'   Internally: F0 is extracted via Harvest (WORLD), then fed to CheapTrick.
##'   The output is a multi-column SSFF track (\code{sp}) with
##'   \code{fft_size/2 + 1} spectral coefficients per frame, where
##'   \code{fft_size} is determined by \code{fs} and \code{f0_floor} following
##'   the WORLD formula.
##'
##'   All input media formats supported by \pkg{av} are accepted.
##'
##' @param listOfFiles Character vector of audio file paths.
##' @param beginTime Numeric. Start of analysis window in seconds. Default 0.
##' @param endTime Numeric. End of analysis window in seconds. Default 0 (file end).
##' @param windowShift Numeric. Frame shift in milliseconds. Default 5.0 ms.
##' @param minF Numeric. Minimum F0 in Hz for the Harvest pitch extractor.
##'   Default 60.0 Hz.
##' @param maxF Numeric. Maximum F0 in Hz for the Harvest pitch extractor.
##'   Default 400.0 Hz.
##' @param voicing_threshold Numeric. Voicing threshold for Harvest (default 0.1).
##' @param q1 Numeric. CheapTrick spectral regularization parameter (default -0.15).
##' @param f0_floor Numeric. Lower F0 bound used to determine FFT size (default 71.0).
##'   For 16 kHz audio, 71 Hz gives fft_size = 2048 (sp_length = 1025).
##' @param toFile Logical. If \code{TRUE}, write SSFF output and return count
##'   invisibly. If \code{FALSE}, return \code{AsspDataObj}. Default \code{TRUE}.
##' @param explicitExt Character. Output file extension. Default \code{"sp"}.
##' @param outputDirectory Character. Output directory. \code{NULL} (default)
##'   writes alongside the input file.
##' @param verbose Logical. Print per-file progress. Default \code{TRUE}.
##'
##' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with track:
##'   \describe{
##'     \item{\code{sp}}{REAL32, power spectral envelope, n_frames x (fft_size/2+1).}
##'   }
##'   Frame rate: \code{1000 / windowShift} Hz.
##'   If \code{toFile = TRUE}: integer count of files written, returned invisibly.
##'
##' @export
##' @examples
##' \dontrun{
##' # Estimate spectral envelope
##' trk_cheap_trick("recording.wav")
##'
##' # Return data without writing file
##' sp <- trk_cheap_trick("recording.wav", toFile = FALSE)
##' dim(sp$sp)  # n_frames x 1025 for 16 kHz audio
##'
##' # Process multiple files
##' trk_cheap_trick(c("file1.wav", "file2.wav"))
##' }
trk_cheap_trick <- function(listOfFiles,
                            beginTime = 0.0,
                            endTime = 0.0,
                            windowShift = 5.0,
                            minF = 60.0,
                            maxF = 400.0,
                            voicing_threshold = 0.1,
                            q1 = -0.15,
                            f0_floor = 71.0,
                            toFile = TRUE,
                            explicitExt = "sp",
                            outputDirectory = NULL,
                            verbose = TRUE) {

  if (is.null(listOfFiles) || length(listOfFiles) == 0)
    cli::cli_abort("No input files specified in {.arg listOfFiles}")

  listOfFiles <- fast_strip_file_protocol(listOfFiles)
  listOfFiles <- normalizePath(path.expand(listOfFiles), mustWork = FALSE)

  files_exist <- file.exists(listOfFiles)
  if (!all(files_exist)) {
    missing_files <- listOfFiles[!files_exist]
    cli::cli_abort(c("!" = "Some files do not exist:",
                     "x" = "{.file {fast_basename(missing_files)}}"))
  }

  n_files <- length(listOfFiles)
  beginTime <- if (is.null(beginTime)) 0.0 else beginTime
  endTime   <- if (is.null(endTime))   0.0 else endTime
  if (length(beginTime) == 1) beginTime <- rep(beginTime, n_files)
  if (length(endTime)   == 1) endTime   <- rep(endTime,   n_files)

  makeOutputDirectory(outputDirectory, FALSE, "trk_cheap_trick")
  if (verbose) format_apply_msg("trk_cheap_trick", n_files, beginTime, endTime)

  results <- vector("list", n_files)

  if (verbose && n_files > 1)
    cli::cli_progress_bar("Processing files", total = n_files,
      format = "{cli::pb_spin} {cli::pb_current}/{cli::pb_total} | ETA: {cli::pb_eta}")

  for (i in seq_len(n_files)) {
    file_path <- listOfFiles[i]
    bt <- beginTime[i]
    et <- endTime[i]

    tryCatch({
      audio_obj <- read_audio(file_path, begin = bt, end = et)

      harvest_result <- harvest_cpp(
        audio_obj = audio_obj,
        minF = minF,
        maxF = maxF,
        windowShift = windowShift,
        voicing_threshold = voicing_threshold,
        verbose = FALSE
      )

      ct_result <- cheap_trick_cpp(
        audio_obj = audio_obj,
        f0 = as.numeric(harvest_result$f0),
        temporal_positions = as.numeric(harvest_result$times),
        q1 = q1,
        f0_floor = f0_floor,
        verbose = FALSE
      )

      out_obj <- create_spectrogram_asspobj(ct_result, windowShift)

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

    if (verbose && n_files > 1) cli::cli_progress_update()
  }

  if (verbose && n_files > 1) cli::cli_progress_done()

  if (toFile) {
    n_success <- sum(unlist(results), na.rm = TRUE)
    if (verbose)
      cli::cli_inform("Successfully processed {n_success} of {n_files} file{?s}")
    return(invisible(n_success))
  } else {
    results <- results[!sapply(results, is.null)]
    if (length(results) == 1) results[[1]] else results
  }
}

attr(trk_cheap_trick, "ext")             <- "sp"
attr(trk_cheap_trick, "tracks")          <- c("sp")
attr(trk_cheap_trick, "outputType")      <- "SSFF"
attr(trk_cheap_trick, "nativeFiletypes") <- c("wav", "flac", "mp3", "mp4", "mkv", "avi")
attr(trk_cheap_trick, "suggestCaching")  <- FALSE
