##' RAPT Pitch Tracking (C++ implementation)
##'
##' @description Extract F0 (fundamental frequency) using the Robust Algorithm
##'   for Pitch Tracking (RAPT) from SPTK. This is a high-performance C++
##'   implementation that is 2-3x faster than the Python version and requires
##'   no Python dependencies.
##'
##'   The algorithm uses dynamic programming to find the optimal F0 contour,
##'   making it robust to noise and reliable for speech analysis.
##'
##'   All input media formats are supported via the av package, including video
##'   files from which audio will be automatically extracted.
##'
##' @param listOfFiles Vector of file paths to process
##' @param beginTime Start time in seconds (default: 0.0)
##' @param endTime End time in seconds (default: 0.0 = end of file)
##' @param windowShift Frame shift in milliseconds (default: 10.0)
##' @param minF Minimum F0 in Hz (default: 60.0)
##' @param maxF Maximum F0 in Hz (default: 400.0)
##' @param voicing_threshold Voicing threshold (default: 0.9)
##' @param toFile Write results to file (default: TRUE)
##' @param explicitExt Output file extension (default: "f0")
##' @param outputDirectory Output directory (default: NULL = same as input)
##' @param verbose Show progress messages (default: TRUE)
##'
##' @return If toFile=TRUE, returns the number of successfully processed files.
##'   If toFile=FALSE, returns AsspDataObj or list of AsspDataObj objects.
##'
##' @export
##' @examples
##' \dontrun{
##' # Extract F0 from audio file
##' trk_rapt("recording.wav")
##'
##' # Process with custom F0 range
##' trk_rapt("speech.mp3", minF = 75, maxF = 300)
##'
##' # Return data without writing file
##' f0_data <- trk_rapt("audio.wav", toFile = FALSE)
##'
##' # Process video file (extracts audio)
##' trk_rapt("interview.mp4")
##' }
trk_rapt <- function(listOfFiles,
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
  makeOutputDirectory(outputDirectory, FALSE, "trk_rapt")

  if (verbose) {
    cli::cli_inform("Applying {.fun rapt} (C++) to {cli::no(n_files)} recording{?s}")
  }

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

attr(trk_rapt, "ext") <- "f0"
attr(trk_rapt, "tracks") <- c("f0")
attr(trk_rapt, "outputType") <- "SSFF"
attr(trk_rapt, "nativeFiletypes") <- c("wav", "flac", "mp3", "mp4", "mkv", "avi")  # Via av
attr(trk_rapt, "suggestCaching") <- FALSE
