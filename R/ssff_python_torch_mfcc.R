#' Extract MFCC features using torchaudio (DEPRECATED)
#'
#' @description
#' **DEPRECATED**: Please use [sptk_mfcc()] instead, which provides a faster
#' C++-based implementation using the SPTK library (2-3x faster) without
#' requiring Python dependencies.
#'
#' Compute Mel-frequency cepstral coefficients (MFCCs) from audio using torchaudio.
#' MFCCs are widely used features in speech recognition, speaker identification,
#' and other audio analysis tasks.
#'
#' The implementation follows the standard MFCC pipeline:
#' 1. Compute power spectrogram via STFT
#' 2. Apply mel-scale filterbank
#' 3. Convert to dB scale
#' 4. Apply Discrete Cosine Transform (DCT)
#'
#' All input media formats are supported via torchaudio, including audio extraction
#' from video files.
#'
#' @param listOfFiles Vector of file paths to process
#' @param beginTime Start time in seconds (default: 0.0)
#' @param endTime End time in seconds (default: 0.0 = end of file)
#' @param windowShift Frame shift in milliseconds (default: 10.0)
#' @param windowSize FFT window size in milliseconds (default: 25.0)
#' @param n_mfcc Number of MFCC coefficients to return (default: 13)
#' @param n_mels Number of mel filterbanks (default: 40)
#' @param fmin Minimum frequency in Hz (default: 0.0)
#' @param fmax Maximum frequency in Hz (default: NULL = sample_rate/2)
#' @param toFile Write results to file (default: TRUE)
#' @param explicitExt Output file extension (default: "mfcc")
#' @param outputDirectory Output directory (default: NULL = same as input)
#' @param verbose Show progress messages (default: TRUE)
#'
#' @return If toFile=TRUE, returns the number of successfully processed files.
#'   If toFile=FALSE, returns AsspDataObj with MFCC tracks (mfcc_1, mfcc_2, ..., mfcc_n).
#'
#' @note **DEPRECATED**: This function requires torch and torchaudio.
#'   Use [sptk_mfcc()] for better performance.
#'
#' @seealso [sptk_mfcc()] for the recommended C++ implementation
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Extract 13 MFCCs (default)
#' mfcc("recording.wav")
#'
#' # Extract 20 MFCCs with custom parameters
#' mfcc("speech.mp3", n_mfcc = 20, n_mels = 80)
#'
#' # Return data without writing file
#' mfcc_data <- mfcc("audio.wav", toFile = FALSE)
#'
#' # Process with specific frequency range
#' mfcc("recording.wav", fmin = 80, fmax = 8000)
#' }
mfcc <- function(listOfFiles = NULL,
                beginTime = 0.0,
                endTime = 0.0,
                windowShift = 10.0,
                windowSize = 25.0,
                n_mfcc = 13,
                n_mels = 40,
                fmin = 0.0,
                fmax = NULL,
                toFile = TRUE,
                explicitExt = "mfcc",
                outputDirectory = NULL,
                verbose = TRUE) {

  # Deprecation warning
  cli::cli_warn(c(
    "!" = "{.fn mfcc} is deprecated and will be removed in a future version.",
    "i" = "Please use {.fn sptk_mfcc} instead for better performance (2-3x faster).",
    "i" = "The C++ implementation requires no Python dependencies."
  ))
  
  # Validate inputs
  if (is.null(listOfFiles) || length(listOfFiles) == 0) {
    cli::cli_abort("No input files specified in {.arg listOfFiles}")
  }
  
  if (n_mfcc < 1 || n_mfcc > 50) {
    cli::cli_abort("{.arg n_mfcc} must be between 1 and 50")
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
  makeOutputDirectory(outputDirectory, FALSE, "mfcc")

  if (verbose) {
    cli::cli_inform("Extracting {n_mfcc} MFCCs from {cli::no(n_files)} recording{?s}")
  }

  # Check that Python script exists
  python_script <- system.file("python", "mfcc.py", package = "superassp")
  if (!file.exists(python_script)) {
    cli::cli_abort(c(
      "x" = "Python script not found: {.file mfcc.py}",
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
      python_script <- system.file("python", "mfcc.py", package = "superassp")
      
      # Prepare parameters as JSON
      params <- list(
        soundFile = file_path,
        windowShift = windowShift,
        windowSize = windowSize,
        n_mfcc = as.integer(n_mfcc),
        n_mels = as.integer(n_mels),
        fmin = fmin,
        fmax = if (is.null(fmax)) NULL else fmax,
        beginTime = bt,
        endTime = et
      )
      params_json <- jsonlite::toJSON(params, auto_unbox = TRUE, null = "null")
      
      # Call Python script
      cmd <- sprintf("python3 '%s' '%s'", python_script, params_json)
      result_json <- system(cmd, intern = TRUE, ignore.stderr = FALSE)
      
      # Parse JSON result
      result <- jsonlite::fromJSON(result_json)
      
      # Extract results
      mfcc_matrix <- matrix(unlist(result$mfcc), ncol = result$n_mfcc, byrow = TRUE)
      sample_rate <- result$sample_rate
      n_frames <- result$n_frames
      
      # Create AsspDataObj
      outDataObj <- list()
      
      sampleRate <- 1000.0 / windowShift
      attr(outDataObj, "sampleRate") <- sampleRate
      attr(outDataObj, "origFreq") <- sample_rate
      
      startTime <- 1 / sampleRate
      attr(outDataObj, "startTime") <- as.numeric(startTime)
      attr(outDataObj, "startRecord") <- as.integer(1)
      attr(outDataObj, "endRecord") <- as.integer(n_frames)
      
      class(outDataObj) <- "AsspDataObj"
      AsspFileFormat(outDataObj) <- "SSFF"
      AsspDataFormat(outDataObj) <- as.integer(2)
      attr(outDataObj, "trackFormats") <- character(0)
      
      # Add MFCC tracks
      for (coef in seq_len(n_mfcc)) {
        track_name <- sprintf("mfcc_%d", coef)
        track_data <- matrix(mfcc_matrix[, coef], ncol = 1)
        outDataObj <- wrassp::addTrack(outDataObj, track_name, track_data, "REAL32")
        attr(outDataObj, "trackFormats") <- c(attr(outDataObj, "trackFormats"), "REAL32")
      }
      
      # Apply Emu-SDMS fix for missing samples at start
      if (startTime > (1 / sampleRate)) {
        nr_of_missing_samples <- as.integer(floor(startTime / (1 / sampleRate)))
        
        for (coef in seq_len(n_mfcc)) {
          track_name <- sprintf("mfcc_%d", coef)
          missing_vals <- matrix(0,
                                nrow = nr_of_missing_samples,
                                ncol = ncol(outDataObj[[track_name]]))
          outDataObj[[track_name]] <- rbind(missing_vals, outDataObj[[track_name]])
        }
        
        attr(outDataObj, "startTime") <- startTime - nr_of_missing_samples * (1 / sampleRate)
      }
      
      assertthat::assert_that(
        wrassp::is.AsspDataObj(outDataObj),
        msg = "The AsspDataObj created by mfcc is invalid."
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

attr(mfcc, "ext") <- "mfcc"
attr(mfcc, "tracks") <- function(n_mfcc = 13) sprintf("mfcc_%d", seq_len(n_mfcc))
attr(mfcc, "outputType") <- "SSFF"
attr(mfcc, "nativeFiletypes") <- c("wav", "flac", "mp3", "mp4", "mkv", "avi")
attr(mfcc, "suggestCaching") <- FALSE
