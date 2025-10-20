#' Estimate pitch using Kaldi-style pitch tracker
#'
#' @description
#' This function estimates pitch using normalized cross-correlation function (NCCF) and
#' median smoothing, similar to the Kaldi ASR toolkit's pitch tracker
#' \insertCite{Ghahremani.2014.10.1109/icassp.2014.6854049}{superassp}.
#'
#' The function uses `torchaudio.functional.detect_pitch_frequency` which replaced
#' the deprecated `compute_kaldi_pitch` (removed in torchaudio 2.9+). The algorithm
#' uses the same NCCF and median smoothing approach as the original Kaldi implementation.
#'
#' **Note**: This implementation only returns F0 values, not the NCCF (Normalized 
#' Cross-Correlation Function) that was available in the original Kaldi implementation
#' with `compute_kaldi_pitch`. If you need voicing probability, consider using
#' [trk_rapt()] or [trk_reaper()] which include voicing confidence measures.
#'
#' All input media formats are supported via torchaudio, including audio extraction
#' from video files.
#'
#' @param listOfFiles Vector of file paths to process
#' @param beginTime Start time in seconds (default: 0.0)
#' @param endTime End time in seconds (default: 0.0 = end of file)
#' @param windowShift Frame shift in milliseconds (default: 10.0)
#' @param windowSize Window length for median smoothing in number of frames (default: 30)
#' @param minF Minimum F0 in Hz (default: 85.0)
#' @param maxF Maximum F0 in Hz (default: 400.0)
#' @param toFile Write results to file (default: TRUE)
#' @param explicitExt Output file extension (default: "kap")
#' @param outputDirectory Output directory (default: NULL = same as input)
#' @param verbose Show progress messages (default: TRUE)
#'
#' @return If toFile=TRUE, returns the number of successfully processed files.
#'   If toFile=FALSE, returns AsspDataObj with F0 track.
#'
#' @note This function requires torchaudio >= 0.13.0.
#'
#' @seealso [trk_rapt()], [trk_swipe()], [trk_reaper()], [trk_crepe()]
#'
#' @references \insertAllCited{}
#' @export
#'
#' @examples
#' \dontrun{
#' # Extract F0 from audio file
#' trk_kaldi_pitch("recording.wav")
#'
#' # Process with custom F0 range
#' trk_kaldi_pitch("speech.mp3", minF = 75, maxF = 300)
#'
#' # Return data without writing file
#' f0_data <- trk_kaldi_pitch("audio.wav", toFile = FALSE)
#'
#' # Process multiple files
#' trk_kaldi_pitch(c("file1.wav", "file2.wav"))
#' }
trk_kaldi_pitch <- function(listOfFiles,
                       beginTime = 0.0,
                       endTime = 0.0,
                       windowShift = 10.0,
                       windowSize = 30,
                       minF = 85.0,
                       maxF = 400.0,
                       toFile = TRUE,
                       explicitExt = "kap",
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
  makeOutputDirectory(outputDirectory, FALSE, "trk_kaldi_pitch")

  if (verbose) {
    cli::cli_inform("Applying {.fun kaldi_pitch} to {cli::no(n_files)} recording{?s}")
  }

  # Check that Python script exists
  python_script <- system.file("python", "kaldi_pitch.py", package = "superassp")
  if (!file.exists(python_script)) {
    cli::cli_abort(c(
      "x" = "Python script not found: {.file kaldi_pitch.py}",
      "i" = "Package may be incorrectly installed"
    ))
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
      # Use external Python script
      python_script <- system.file("python", "kaldi_pitch.py", package = "superassp")
      
      # Prepare parameters as JSON
      params <- list(
        soundFile = file_path,
        windowShift = windowShift,
        windowSize = as.integer(windowSize),
        minF = minF,
        maxF = maxF,
        beginTime = bt,
        endTime = et
      )
      params_json <- jsonlite::toJSON(params, auto_unbox = TRUE)
      
      # Call Python script
      cmd <- sprintf("python3 '%s' '%s'", python_script, params_json)
      result_json <- system(cmd, intern = TRUE, ignore.stderr = FALSE)
      
      # Parse JSON result
      result <- jsonlite::fromJSON(result_json)
      
      # Extract results
      f0_values <- result$f0
      sample_rate <- result$sample_rate
      
      # Create AsspDataObj
      inTable <- data.frame(f0 = f0_values)
      
      outDataObj <- list()
      attr(outDataObj, "trackFormats") <- c("INT16")
      
      sampleRate <- 1000.0 / windowShift
      attr(outDataObj, "sampleRate") <- sampleRate
      attr(outDataObj, "origFreq") <- sample_rate
      
      startTime <- 1 / sampleRate
      attr(outDataObj, "startTime") <- as.numeric(startTime)
      attr(outDataObj, "startRecord") <- as.integer(1)
      attr(outDataObj, "endRecord") <- as.integer(nrow(inTable))
      
      class(outDataObj) <- "AsspDataObj"
      AsspFileFormat(outDataObj) <- "SSFF"
      AsspDataFormat(outDataObj) <- as.integer(2)
      
      # Add F0 track - convert to INT16 and replace NA with 0
      f0_int <- as.integer(inTable$f0)
      f0_int[is.na(f0_int)] <- 0L
      f0_matrix <- matrix(f0_int, ncol = 1)
      
      outDataObj <- wrassp::addTrack(outDataObj, "f0", f0_matrix, "INT16")
      
      # Apply Emu-SDMS fix for missing samples at start
      if (startTime > (1 / sampleRate)) {
        nr_of_missing_samples <- as.integer(floor(startTime / (1 / sampleRate)))
        
        missing_f0_vals <- matrix(0,
                                  nrow = nr_of_missing_samples,
                                  ncol = ncol(outDataObj$f0))
        
        outDataObj$f0 <- rbind(missing_f0_vals, outDataObj$f0)
        attr(outDataObj, "startTime") <- startTime - nr_of_missing_samples * (1 / sampleRate)
      }
      
      assertthat::assert_that(
        wrassp::is.AsspDataObj(outDataObj),
        msg = "The AsspDataObj created by kaldi_pitch is invalid."
      )
      
      # Handle output
      if (toFile) {
        ssff_file <- sub("wav$", explicitExt, basename(file_path))
        if (!is.null(outputDirectory)) {
          ssff_file <- file.path(outputDirectory, ssff_file)
        } else {
          ssff_file <- file.path(dirname(file_path), ssff_file)
        }
        
        attr(outDataObj, "filePath") <- as.character(ssff_file)
        wrassp::write.AsspDataObj(dobj = outDataObj, file = ssff_file)
        results[[i]] <- TRUE
      } else {
        results[[i]] <- outDataObj
      }

    }, error = function(e) {
      if (verbose) {
        cli::cli_warn("Error processing {.file {basename(file_path)}}: {conditionMessage(e)}")
      }
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

attr(trk_kaldi_pitch, "ext") <- "kap"
attr(trk_kaldi_pitch, "tracks") <- c("f0")
attr(trk_kaldi_pitch, "outputType") <- "SSFF"
attr(trk_kaldi_pitch, "nativeFiletypes") <- c("wav", "flac", "mp3", "mp4", "mkv", "avi")
attr(trk_kaldi_pitch, "suggestCaching") <- FALSE
