##' Track fundamental frequency using RAPT (Robust Algorithm for Pitch Tracking)
##'
##' Extracts F0 using the RAPT dynamic-programming pitch tracker from SPTK. RAPT
##' is a normalized cross-correlation tracker with robust voiced/unvoiced decisions
##' and is a reliable general-purpose choice for speech with moderate noise. Prefer
##' SWIPE for cleaner but noisier signals, or PDA for higher temporal resolution.
##'
##' @param listOfFiles Character vector of audio file paths. Any format supported by
##'   \pkg{av} is accepted; non-native inputs are transcoded automatically.
##' @param beginTime Numeric. Start of analysis window in seconds. Default 0 (file start).
##' @param endTime Numeric. End of analysis window in seconds. Default 0 (file end).
##' @param windowShift Numeric. Frame shift in milliseconds; sets output frame rate
##'   (1000 / windowShift Hz). Default 10.0 ms.
##' @param minF Numeric. Minimum F0 in Hz. Default 60.0 Hz.
##' @param maxF Numeric. Maximum F0 in Hz. Default 400.0 Hz.
##' @param voicing_threshold Numeric. Voicing decision threshold (0–1; higher = more
##'   conservative, fewer voiced frames). Default 0.6.
##' @param toFile Logical. If \code{TRUE}, write SSFF output files and return the
##'   count written invisibly. If \code{FALSE}, return an \code{AsspDataObj}.
##'   Default \code{TRUE}.
##' @param explicitExt Character. Output file extension. Default \code{"f0"}.
##' @param outputDirectory Character. Directory for output files. \code{NULL} (default)
##'   writes alongside the input file.
##' @param verbose Logical. Print per-file progress. Default \code{TRUE}.
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
##' trk_pitch_rapt("recording.wav")
##'
##' # Process with custom F0 range
##' trk_pitch_rapt("speech.mp3", minF = 75, maxF = 300)
##'
##' # Return data without writing file
##' f0_data <- trk_pitch_rapt("audio.wav", toFile = FALSE)
##'
##' # Process video file (extracts audio)
##' trk_pitch_rapt("interview.mp4")
##' }
trk_pitch_rapt <- function(listOfFiles,
                 beginTime = 0.0,
                 endTime = 0.0,
                 windowShift = 10.0,
                 minF = 60.0,
                 maxF = 400.0,
                 voicing_threshold = 0.6,
                 toFile = TRUE,
                 explicitExt = "f0",
                 outputDirectory = NULL,
                 verbose = TRUE) {

  # Validate inputs
  if (is.null(listOfFiles) || length(listOfFiles) == 0) {
    cli::cli_abort("No input files specified in {.arg listOfFiles}")
  }

  # Normalize paths
  listOfFiles <- fast_strip_file_protocol(listOfFiles)
  listOfFiles <- normalizePath(path.expand(listOfFiles), mustWork = FALSE)

  # Check file existence
  files_exist <- file.exists(listOfFiles)
  if (!all(files_exist)) {
    missing_files <- listOfFiles[!files_exist]
    cli::cli_abort(c(
      "!" = "Some files do not exist:",
      "x" = "{.file {fast_basename(missing_files)}}"
    ))
  }

  n_files <- length(listOfFiles)

  # Normalize time parameters
  beginTime <- if (is.null(beginTime)) 0.0 else beginTime
  endTime <- if (is.null(endTime)) 0.0 else endTime

  # Recycle time parameters
  if (length(beginTime) == 1) beginTime <- rep(beginTime, n_files)
  if (length(endTime) == 1) endTime <- rep(endTime, n_files)

  # Setup output directory
  makeOutputDirectory(outputDirectory, FALSE, "trk_pitch_rapt")

  if (verbose) format_apply_msg("trk_pitch_rapt", n_files, beginTime, endTime)

  # Process each file
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
      # Load audio with av
      audio_obj <- read_audio(
        file_path,
        begin = bt,
        end   = et
      )

      # Call C++ RAPT
      rapt_result <- rapt_cpp(
        audio_obj = audio_obj,
        minF = minF,
        maxF = maxF,
        windowShift = windowShift,
        voicing_threshold = voicing_threshold,
        verbose = FALSE
      )

      # Convert to AsspDataObj
      out_obj <- create_f0_asspobj(rapt_result, windowShift)

      # Handle output
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

  # Return results
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

attr(trk_pitch_rapt, "ext") <- "f0"
attr(trk_pitch_rapt, "tracks") <- c("f0")
attr(trk_pitch_rapt, "outputType") <- "SSFF"
attr(trk_pitch_rapt, "nativeFiletypes") <- c("wav", "flac", "mp3", "mp4", "mkv", "avi")  # Via av
attr(trk_pitch_rapt, "suggestCaching") <- FALSE
