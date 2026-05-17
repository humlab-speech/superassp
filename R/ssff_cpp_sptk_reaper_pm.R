##' Detect glottal closure instants using REAPER (pitch marks)
##'
##' Extracts glottal closure instants (GCIs) as a binary indicator track using the
##' REAPER EpochTracker \insertCite{talkin2019reaper}{superassp} from SPTK. GCIs are
##' mapped from irregular epoch times onto a regular grid at \code{windowShift}
##' intervals. Use \code{\link{trk_pitch_reaper}} instead when F0 is also needed
##' (avoids re-running REAPER).
##'
##' @param listOfFiles Character vector of audio file paths. Any format supported by
##'   \pkg{av} is accepted; non-native inputs are transcoded automatically.
##' @param beginTime Numeric. Start of analysis window in seconds. Default 0 (file start).
##' @param endTime Numeric. End of analysis window in seconds. Default 0 (file end).
##' @param windowShift Numeric. Frame shift in milliseconds for the output indicator grid.
##'   Default 10.0 ms.
##' @param minF Numeric. Minimum F0 in Hz for the internal pitch estimator. Lower values
##'   allow lower-pitched voices but may increase false positives. Default 40.0 Hz.
##' @param maxF Numeric. Maximum F0 in Hz for the internal pitch estimator. Default 500.0 Hz.
##' @param voicing_threshold Numeric. Voicing decision threshold (0–1; higher = more
##'   conservative). Default 0.9.
##' @param toFile Logical. If \code{TRUE}, write SSFF output files and return the
##'   count written invisibly. If \code{FALSE}, return an \code{AsspDataObj}.
##'   Default \code{TRUE}.
##' @param explicitExt Character. Output file extension. Default \code{"rpm"}.
##' @param outputDirectory Character. Directory for output files. \code{NULL} (default)
##'   writes alongside the input file.
##' @param verbose Logical. Print per-file progress. Default \code{TRUE}.
##'
##' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with track:
##'   \describe{
##'     \item{\code{pm}}{INT16, binary pitch mark indicator (0 = no GCI, 1 = GCI),
##'       n_frames × 1. Frame spacing = windowShift ms.}
##'   }
##'   Additional attributes on the returned object: \code{epoch_times} (raw GCI times
##'   in seconds), \code{n_epochs}, \code{polarity}.
##'   Frame rate: \code{1000 / windowShift} Hz (default 100 Hz).
##'   If \code{toFile = TRUE}: integer count of files written, returned invisibly.
##'
##' @details
##' Raw epoch times at irregular intervals are available via
##' \code{attr(result, "epoch_times")} when \code{toFile = FALSE}.
##'
##' @references \insertAllCited{}
##'
##' @seealso
##' \code{\link{trk_pitch_reaper}} for F0 extraction (also extracts epochs as attributes)
##' \code{\link{trk_pitchmark_estk}} for ESTK-based pitch mark detection
##'
##' @export
##' @examples
##' \dontrun{
##' # Basic usage - extract pitch marks
##' trk_pitchmark_reaper("speech.wav")
##'
##' # Get pitch marks without writing to file
##' result <- trk_pitchmark_reaper("speech.wav", toFile = FALSE)
##' pm_track <- result$pm  # Binary indicator (0 or 1)
##'
##' # Access raw epoch times (irregular intervals)
##' epoch_times <- attr(result, "epoch_times")  # Times in seconds
##' n_epochs <- attr(result, "n_epochs")        # Number of epochs
##'
##' # Adjust F0 range for low-pitched voice
##' trk_pitchmark_reaper("bass_voice.wav", minF = 50, maxF = 300)
##'
##' # Process specific time window
##' trk_pitchmark_reaper("long_recording.wav",
##'               beginTime = 1.0,
##'               endTime = 5.0)
##'
##' # Batch processing
##' files <- list.files(pattern = "\\.wav$", full.names = TRUE)
##' n_success <- trk_pitchmark_reaper(files,
##'                             outputDirectory = "results/",
##'                             verbose = TRUE)
##' message("Processed ", n_success, " files")
##' }
trk_pitchmark_reaper <- function(listOfFiles,
                          beginTime = 0.0,
                          endTime = 0.0,
                          windowShift = 10.0,
                          minF = 40.0,
                          maxF = 500.0,
                          voicing_threshold = 0.9,
                          toFile = TRUE,
                          explicitExt = "rpm",
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

  makeOutputDirectory(outputDirectory, FALSE, "trk_pitchmark_reaper")

  if (verbose) format_apply_msg("trk_pitchmark_reaper", n_files, beginTime, endTime)

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
      # Load audio via av package (supports all media formats)
      audio_obj <- read_audio(
        file_path,
        begin = bt,
        end   = et
      )

      # Call C++ REAPER implementation
      reaper_result <- reaper_cpp(
        audio_obj = audio_obj,
        minF = minF,
        maxF = maxF,
        windowShift = windowShift,
        voicing_threshold = voicing_threshold,
        verbose = FALSE
      )

      # Extract epochs (pitch marks) - already computed by C++!
      epoch_times <- reaper_result$epochs
      n_epochs <- reaper_result$n_epochs
      sample_rate <- reaper_result$sample_rate
      polarity <- reaper_result$polarity

      # Convert epochs to binary pitch mark grid
      out_obj <- create_pitchmark_asspobj(epoch_times, sample_rate, windowShift)

      # Store polarity as attribute (from REAPER analysis)
      attr(out_obj, "polarity") <- polarity

      if (toFile) {
        out_file <- generate_output_path(file_path, explicitExt, outputDirectory)
        write.AsspDataObj(out_obj, out_file)
        results[[i]] <- TRUE

        if (verbose && n_files == 1) {
          cli::cli_inform("  {cli::symbol$tick} Written to {.file {basename(out_file)}} ({n_epochs} epochs)")
        }
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

attr(trk_pitchmark_reaper, "ext") <- "rpm"
attr(trk_pitchmark_reaper, "tracks") <- c("pm")
attr(trk_pitchmark_reaper, "outputType") <- "SSFF"
attr(trk_pitchmark_reaper, "nativeFiletypes") <- c("wav", "flac", "mp3", "mp4", "mkv", "avi")
attr(trk_pitchmark_reaper, "suggestCaching") <- FALSE
