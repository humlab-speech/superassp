#' Snack-style Pitch Tracking
#'
#' Extracts fundamental frequency (F0) using the Snack pitch tracker algorithm.
#' This implementation replicates the pitch tracking method from the Snack Sound
#' Toolkit, which uses autocorrelation with dynamic programming for robust F0 estimation.
#'
#' @param listOfFiles Vector of file paths to audio files
#' @param beginTime Start time in seconds (default: 0.0)
#' @param endTime End time in seconds (default: 0.0 = end of file)
#' @param windowShift Frame shift in milliseconds (default: 10.0)
#' @param windowLength Analysis window length in milliseconds (default: 7.5)
#' @param minF Minimum F0 in Hz (default: 50)
#' @param maxF Maximum F0 in Hz (default: 550)
#' @param threshold Correlation threshold for peak detection (default: 0.3)
#' @param toFile Write results to file (default: TRUE)
#' @param explicitExt File extension for output files (default: "snackpitch")
#' @param outputDirectory Directory for output files (default: NULL = same as input)
#' @param verbose Print progress messages (default: TRUE)
#'
#' @return If toFile=FALSE, returns AsspDataObj or list of AsspDataObjs with tracks:
#'   \itemize{
#'     \item f0: fundamental frequency in Hz (0 = unvoiced)
#'     \item voicing: voicing probability (0-1)
#'     \item rms: RMS energy per frame
#'   }
#'   If toFile=TRUE, returns number of files successfully processed.
#'
#' @details
#' The Snack pitch tracker uses normalized autocorrelation to find pitch candidates
#' at each frame, then applies dynamic programming to select the most likely F0
#' trajectory. This approach balances local pitch estimates with trajectory smoothness.
#'
#' Default parameters match Snack's defaults:
#' \itemize{
#'   \item minF = 50 Hz, maxF = 550 Hz
#'   \item windowShift = 10 ms
#'   \item windowLength = 7.5 ms
#'   \item threshold = 0.3 (correlation threshold)
#' }
#'
#' This implementation provides compatibility with analyses using Snack as a reference.
#'
#' @references
#' Sjölander, K. & Beskow, J. (2000). "Wavesurfer - an open source speech tool."
#' In Proc. ICSLP 2000, Beijing, China.
#'
#' @seealso \code{\link{rapt}}, \code{\link{swipe}}, \code{\link{kaldi_pitch}}
#'   for alternative pitch tracking methods
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' snack_pitch("audio.wav")
#'
#' # Return data without file
#' pitch_data <- snack_pitch("audio.wav", toFile = FALSE)
#' plot(pitch_data$f0, type = "l")
#'
#' # Custom parameters
#' snack_pitch("speech.wav", minF = 75, maxF = 400, windowShift = 5)
#'
#' # Batch processing
#' files <- list.files("audio_dir", pattern = ".wav$", full.names = TRUE)
#' snack_pitch(files, verbose = TRUE)
#' }
#'
#' @export
snack_pitch <- function(listOfFiles = NULL,
                       beginTime = 0.0,
                       endTime = 0.0,
                       windowShift = 10.0,
                       windowLength = 7.5,
                       minF = 50,
                       maxF = 550,
                       threshold = 0.3,
                       toFile = TRUE,
                       explicitExt = "snackpitch",
                       outputDirectory = NULL,
                       verbose = TRUE) {
  
  # Input validation
  if (is.null(listOfFiles)) {
    cli::cli_abort("No files provided")
  }
  
  listOfFiles <- as.vector(listOfFiles)
  n_files <- length(listOfFiles)
  
  # Check files exist
  missing <- !file.exists(listOfFiles)
  if (any(missing)) {
    cli::cli_abort(c(
      "x" = "{sum(missing)} file{?s} not found",
      "i" = "First missing: {.file {listOfFiles[which(missing)[1]]}}"
    ))
  }
  
  # Normalize time parameters
  beginTime <- if (is.null(beginTime)) rep(0.0, n_files) else beginTime
  endTime <- if (is.null(endTime)) rep(0.0, n_files) else endTime
  
  if (length(beginTime) == 1) beginTime <- rep(beginTime, n_files)
  if (length(endTime) == 1) endTime <- rep(endTime, n_files)
  
  if (verbose) {
    cli::cli_inform("Applying {.fun snack_pitch} to {cli::no(n_files)} recording{?s}")
  }
  
  # Check that Python script exists
  python_script <- system.file("python", "snack_pitch.py", package = "superassp")
  if (!file.exists(python_script)) {
    cli::cli_abort(c(
      "x" = "Python script not found: {.file snack_pitch.py}",
      "i" = "Package may be incorrectly installed"
    ))
  }
  
  # Process each file
  results <- vector("list", n_files)
  n_success <- 0
  
  if (verbose && n_files > 1) {
    cli::cli_progress_bar("Processing files", total = n_files)
  }
  
  for (i in seq_len(n_files)) {
    file_path <- listOfFiles[i]
    bt <- beginTime[i]
    et <- endTime[i]
    
    tryCatch({
      # Use external Python script
      params <- list(
        soundFile = file_path,
        minF = minF,
        maxF = maxF,
        windowShift = windowShift,
        windowLength = windowLength,
        threshold = threshold,
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
      voicing_values <- result$voicing
      rms_values <- result$rms
      sample_rate <- result$sample_rate
      
      # Create AsspDataObj
      sampleRate <- 1000.0 / windowShift
      n_frames <- length(f0_values)
      
      outDataObj <- list()
      attr(outDataObj, "trackFormats") <- c("INT16", "REAL32", "REAL32")
      attr(outDataObj, "sampleRate") <- sampleRate
      attr(outDataObj, "origFreq") <- sample_rate
      
      startTime <- 1 / sampleRate
      attr(outDataObj, "startTime") <- as.numeric(startTime)
      attr(outDataObj, "startRecord") <- as.integer(1)
      attr(outDataObj, "endRecord") <- as.integer(n_frames)
      
      class(outDataObj) <- "AsspDataObj"
      AsspFileFormat(outDataObj) <- "SSFF"
      AsspDataFormat(outDataObj) <- as.integer(2)
      
      # Add F0 track (INT16)
      f0_int <- as.integer(f0_values)
      f0_int[is.na(f0_int)] <- 0L
      f0_matrix <- matrix(f0_int, ncol = 1)
      outDataObj <- wrassp::addTrack(outDataObj, "f0", f0_matrix, "INT16")
      
      # Add voicing track (REAL32)
      voicing_matrix <- matrix(as.numeric(voicing_values), ncol = 1)
      outDataObj <- wrassp::addTrack(outDataObj, "voicing", voicing_matrix, "REAL32")
      
      # Add RMS track (REAL32)
      rms_matrix <- matrix(as.numeric(rms_values), ncol = 1)
      outDataObj <- wrassp::addTrack(outDataObj, "rms", rms_matrix, "REAL32")
      
      # Apply Emu-SDMS fix for missing samples at start
      if (startTime > (1 / sampleRate)) {
        outDataObj <- applyWaveSampStartEmuSDMSfix(outDataObj, bt)
      }
      
      # Handle output
      if (toFile) {
        # Determine output path
        if (is.null(outputDirectory)) {
          out_path <- paste0(tools::file_path_sans_ext(file_path), ".", explicitExt)
        } else {
          base_name <- basename(tools::file_path_sans_ext(file_path))
          out_path <- file.path(outputDirectory, paste0(base_name, ".", explicitExt))
        }
        
        # Write SSFF file
        wrassp::write.AsspDataObj(outDataObj, out_path)
        n_success <- n_success + 1
        results[[i]] <- out_path
      } else {
        results[[i]] <- outDataObj
        n_success <- n_success + 1
      }
      
    }, error = function(e) {
      cli::cli_warn("Error processing {.file {basename(file_path)}}: {conditionMessage(e)}")
      results[[i]] <- NULL
    })
    
    if (verbose && n_files > 1) {
      cli::cli_progress_update()
    }
  }
  
  if (verbose && n_files > 1) {
    cli::cli_progress_done()
  }
  
  if (verbose) {
    cli::cli_inform("Successfully processed {n_success} of {n_files} file{?s}")
  }
  
  # Return results
  if (toFile) {
    return(n_success)
  } else {
    if (n_files == 1) {
      return(results[[1]])
    } else {
      return(results)
    }
  }
}

# Set function attributes
attr(snack_pitch, "ext") <- "snackpitch"
attr(snack_pitch, "tracks") <- c("f0", "voicing", "rms")
attr(snack_pitch, "outputType") <- "SSFF"
