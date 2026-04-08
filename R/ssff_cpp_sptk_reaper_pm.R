##' Extract Pitch Marks using REAPER (C++ implementation)
##'
##' @description Extract glottal closure instants (pitch marks) using the Robust
##'   Epoch And Pitch EstimatoR (REAPER) from SPTK. This C++ implementation is
##'   2-3x faster than the Python-based \code{reaper_pm()} function.
##'
##'   REAPER \insertCite{talkin2019reaper}{superassp} uses an EpochTracker to
##'   simultaneously estimate voiced-speech epochs (glottal closure instants),
##'   voicing state, and F0. This function returns pitch marks as a binary
##'   indicator track at regular intervals defined by \code{windowShift}.
##'
##' @param listOfFiles Character vector of file paths to audio files, or AVAudio
##'   S7 object. Supports all media formats (WAV, MP3, MP4, video files, etc.)
##'   via the av package.
##' @param beginTime Start time in seconds for processing (default: 0.0 = start of file)
##' @param endTime End time in seconds for processing (default: 0.0 = end of file)
##' @param windowShift Frame shift in milliseconds for the output grid (default: 10.0 ms).
##'   Pitch marks are mapped to this regular grid as binary indicators.
##' @param minF Minimum F0 in Hz (default: 40.0). Lower values allow tracking
##'   of lower-pitched voices but may increase false positives.
##' @param maxF Maximum F0 in Hz (default: 500.0). Higher values allow tracking
##'   of higher-pitched voices.
##' @param voicing_threshold Voicing decision threshold between 0 and 1 (default: 0.9).
##'   Higher values make voicing decisions more conservative. Equivalent to
##'   \code{unvoiced_cost} in Python pyreaper.
##' @param toFile If TRUE (default), write SSFF files to disk and return count of
##'   successful files. If FALSE, return AsspDataObj (single file) or list of
##'   AsspDataObj (multiple files).
##' @param explicitExt Output file extension (default: "rpm" = reaper pitch marks)
##' @param outputDirectory Output directory for SSFF files. If NULL (default),
##'   files are written to same directory as input files.
##' @param verbose If TRUE (default), print progress messages.
##'
##' @return
##'   If \code{toFile = TRUE}: Returns the number of successfully processed files.
##'   If \code{toFile = FALSE}: Returns an AsspDataObj (single file) or list of
##'   AsspDataObj objects (multiple files) with:
##'   \itemize{
##'     \item \strong{pm}: Binary pitch mark indicator (INT16). Values: 0 (no pitch mark),
##'       1 (pitch mark present at this time frame)
##'     \item \strong{Attributes}: \code{epoch_times} (raw epoch times in seconds),
##'       \code{n_epochs} (number of epochs), \code{polarity} (signal polarity)
##'   }
##'
##' @details
##' \strong{Algorithm:}
##'
##' REAPER uses a two-pass algorithm:
##' \enumerate{
##'   \item \strong{Epoch Detection}: Identifies potential glottal closure instants (GCIs)
##'     using normalized correlation coefficient peaks
##'   \item \strong{F0 Tracking}: Refines epoch selection using dynamic programming to
##'     find the most likely F0 contour
##'   \item \strong{Polarity Detection}: Determines if signal should be inverted
##'   \item \strong{Grid Mapping}: Maps irregular epoch times to regular grid at
##'     \code{windowShift} intervals as binary indicators
##' }
##'
##' \strong{Output Format:}
##'
##' The \code{pm} track contains binary values (0 or 1) at regular intervals defined
##' by \code{windowShift}:
##' \itemize{
##'   \item \strong{0}: No pitch mark detected in this frame
##'   \item \strong{1}: Pitch mark (glottal closure instant) detected in this frame
##' }
##'
##' For access to raw epoch times (irregular intervals), use:
##' \code{attr(result, "epoch_times")} when \code{toFile = FALSE}.
##'
##' \strong{Performance:}
##'
##' This C++ implementation is approximately 2-3x faster than the Python-based
##' \code{reaper_pm()} and requires no Python dependencies.
##'
##' \strong{Comparison with reaper_pm():}
##'
##' \tabular{lll}{
##'   \strong{Aspect} \tab \strong{trk_reaper_pm (C++)} \tab \strong{reaper_pm (Python)} \cr
##'   Backend \tab SPTK C++ \tab pyreaper Python \cr
##'   Speed \tab 2-3x faster \tab Baseline \cr
##'   Dependencies \tab None (built-in) \tab Python + pyreaper \cr
##'   Output \tab Identical \tab Binary pm track \cr
##'   Parameters \tab Simplified \tab More options
##' }
##'
##' @references
##' \insertAllCited{}
##'
##' @seealso
##' \code{\link{trk_reaper}} for F0 extraction (also extracts epochs as attributes)
##' \code{\link{trk_pitchmark}} for ESTK-based pitch mark detection
##'
##' @export
##' @examples
##' \dontrun{
##' # Basic usage - extract pitch marks
##' trk_reaper_pm("speech.wav")
##'
##' # Get pitch marks without writing to file
##' result <- trk_reaper_pm("speech.wav", toFile = FALSE)
##' pm_track <- result$pm  # Binary indicator (0 or 1)
##'
##' # Access raw epoch times (irregular intervals)
##' epoch_times <- attr(result, "epoch_times")  # Times in seconds
##' n_epochs <- attr(result, "n_epochs")        # Number of epochs
##'
##' # Adjust F0 range for low-pitched voice
##' trk_reaper_pm("bass_voice.wav", minF = 50, maxF = 300)
##'
##' # Process specific time window
##' trk_reaper_pm("long_recording.wav",
##'               beginTime = 1.0,
##'               endTime = 5.0)
##'
##' # Batch processing
##' files <- list.files(pattern = "\\.wav$", full.names = TRUE)
##' n_success <- trk_reaper_pm(files,
##'                             outputDirectory = "results/",
##'                             verbose = TRUE)
##' message("Processed ", n_success, " files")
##'
##' # For voice source analysis - get both F0 and epochs
##' f0_result <- trk_reaper("speech.wav", toFile = FALSE)
##' pm_result <- trk_reaper_pm("speech.wav", toFile = FALSE)
##'
##' # Compare regular F0 with pitch mark density
##' plot(f0_result$f0, type = "l", main = "F0 vs Pitch Marks")
##' points(which(pm_result$pm == 1), rep(100, sum(pm_result$pm)),
##'        col = "red", pch = "|")
##' }
trk_reaper_pm <- function(listOfFiles,
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

  makeOutputDirectory(outputDirectory, FALSE, "trk_reaper_pm")

  if (verbose) format_apply_msg("trk_reaper_pm", n_files, beginTime, endTime)

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

attr(trk_reaper_pm, "ext") <- "rpm"
attr(trk_reaper_pm, "tracks") <- c("pm")
attr(trk_reaper_pm, "outputType") <- "SSFF"
attr(trk_reaper_pm, "nativeFiletypes") <- c("wav", "flac", "mp3", "mp4", "mkv", "avi")
attr(trk_reaper_pm, "suggestCaching") <- FALSE
