##' REAPER Pitch and Epoch Tracking (C++ implementation)
##'
##' @description Extract F0 and glottal closure instants using the Robust Epoch
##'   And Pitch EstimatoR (REAPER) from SPTK. This C++ implementation is 2-3x
##'   faster than the Python version.
##'
##'   REAPER provides both F0 estimates and epoch marks (glottal closure instants),
##'   making it useful for voice source analysis.
##'
##' @inheritParams rapt
##' @param voicing_threshold Voicing threshold (default: 0.9)
##'
##' @return If toFile=TRUE, returns the number of successfully processed files.
##'   If toFile=FALSE, returns AsspDataObj or list of AsspDataObj objects with
##'   f0 track and optionally epochs.
##'
##' @export
##' @examples
##' \dontrun{
##' # Extract F0 and epochs
##' trk_reaper("recording.wav")
##'
##' # Get epochs for voice source analysis
##' result <- trk_reaper("speech.wav", toFile = FALSE)
##' epochs <- attr(result, "epochs")
##' }
trk_reaper <- function(listOfFiles,
                   beginTime = 0.0,
                   endTime = 0.0,
                   windowShift = 10.0,
                   minF = 60.0,
                   maxF = 400.0,
                   voicing_threshold = 0.9,
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

  makeOutputDirectory(outputDirectory, FALSE, "trk_reaper")

  if (verbose) {
    cli::cli_inform("Applying {.fun reaper} (C++) to {cli::no(n_files)} recording{?s}")
  }

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
      audio_obj <- av_to_asspDataObj(
        file_path,
        start_time = bt,
        end_time = if (et == 0.0) NULL else et
      )

      reaper_result <- reaper_cpp(
        audio_obj = audio_obj,
        minF = minF,
        maxF = maxF,
        windowShift = windowShift,
        voicing_threshold = voicing_threshold,
        verbose = FALSE
      )

      out_obj <- create_f0_asspobj(reaper_result, windowShift)

      # Store epochs and polarity as attributes
      attr(out_obj, "epochs") <- reaper_result$epochs
      attr(out_obj, "n_epochs") <- reaper_result$n_epochs
      attr(out_obj, "polarity") <- reaper_result$polarity

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

attr(trk_reaper, "ext") <- "f0"
attr(trk_reaper, "tracks") <- c("f0")
attr(trk_reaper, "outputType") <- "SSFF"
attr(trk_reaper, "nativeFiletypes") <- c("wav", "flac", "mp3", "mp4", "mkv", "avi")
attr(trk_reaper, "suggestCaching") <- FALSE
