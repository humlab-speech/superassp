##' D4C Aperiodicity Estimation (C++ implementation)
##'
##' @description Estimate band aperiodicity using the D4C algorithm from WORLD vocoder (via SPTK).
##'   D4C provides high-quality aperiodicity estimation for speech synthesis and analysis.
##'   This implementation uses direct C++ calls for optimal performance.
##'
##' @param listOfFiles Vector of file paths to process
##' @param beginTime Start time in seconds (default: 0.0)
##' @param endTime End time in seconds (default: 0.0, meaning end of file)
##' @param windowShift Frame shift in milliseconds (default: 5.0)
##' @param minF Minimum F0 in Hz (default: 60.0)
##' @param maxF Maximum F0 in Hz (default: 400.0)
##' @param voicing_threshold Voicing threshold for F0 detection (default: 0.85)
##' @param threshold D4C threshold parameter (default: 0.85)
##' @param toFile Write results to file (default: TRUE)
##' @param explicitExt Output file extension (default: "ap")
##' @param outputDirectory Output directory for files (default: NULL, same as input)
##' @param verbose Print progress messages (default: TRUE)
##'
##' @return If toFile=TRUE, returns the number of successfully processed files.
##'   If toFile=FALSE, returns AsspDataObj or list of AsspDataObj objects.
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
                voicing_threshold = 0.85,
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

  if (verbose) {
    cli::cli_inform("Applying {.fun d4c} (C++) to {cli::no(n_files)} recording{?s}")
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

      d4c_result <- d4c_cpp(
        audio_obj = audio_obj,
        minF = minF,
        maxF = maxF,
        windowShift = windowShift,
        voicing_threshold = voicing_threshold,
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
