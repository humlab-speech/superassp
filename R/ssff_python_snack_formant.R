#' Snack-style Formant Tracking
#'
#' Extracts formant frequencies and bandwidths using the Snack formant tracker algorithm.
#' This implementation replicates the formant tracking method from the Snack Sound
#' Toolkit, which uses LPC analysis with dynamic programming for robust formant estimation.
#'
#' @param listOfFiles Vector of file paths to audio files
#' @param beginTime Start time in seconds (default: 0.0)
#' @param endTime End time in seconds (default: 0.0 = end of file)
#' @param windowShift Frame shift in milliseconds (default: 5.0)
#' @param windowLength Analysis window length in milliseconds (default: 49.0)
#' @param numFormants Number of formants to track (default: 4, max: 7)
#' @param lpcOrder LPC order (default: NULL = 2*numFormants + 6)
#' @param preEmphasis Pre-emphasis factor (default: 0.7)
#' @param toFile Write results to file (default: TRUE)
#' @param explicitExt File extension for output files (default: "snackfmt")
#' @param outputDirectory Directory for output files (default: NULL = same as input)
#' @param verbose Print progress messages (default: TRUE)
#'
#' @return If toFile=FALSE, returns AsspDataObj or list of AsspDataObjs with tracks:
#'   \itemize{
#'     \item fm_1, fm_2, ..., fm_N: formant frequencies in Hz
#'     \item bw_1, bw_2, ..., bw_N: formant bandwidths in Hz
#'   }
#'   If toFile=TRUE, returns number of files successfully processed.
#'
#' @details
#' The Snack formant tracker uses Linear Predictive Coding (LPC) to estimate
#' vocal tract resonances. The algorithm:
#' \enumerate{
#'   \item Applies pre-emphasis to boost high frequencies
#'   \item Computes LPC coefficients using autocorrelation method
#'   \item Finds polynomial roots to extract pole frequencies
#'   \item Maps poles to formants using expected frequency ranges
#'   \item Applies trajectory smoothing
#' }
#'
#' Default parameters match Snack's defaults:
#' \itemize{
#'   \item numFormants = 4
#'   \item lpcOrder = 2 * numFormants + 6 = 14
#'   \item windowShift = 5 ms
#'   \item windowLength = 49 ms
#'   \item preEmphasis = 0.7
#' }
#'
#' This implementation provides compatibility with analyses using Snack as a reference.
#'
#' @references
#' Sjölander, K. & Beskow, J. (2000). "Wavesurfer - an open source speech tool."
#' In Proc. ICSLP 2000, Beijing, China.
#'
#' Talkin, D. (1987). "Speech formant trajectory estimation using dynamic programming
#' with modulated transition costs." J. Acoust. Soc. Am.
#'
#' @seealso \code{\link{praat_formant_burg}} for Praat-style formant tracking
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' trk_snackf("audio.wav")
#'
#' # Return data without file
#' formant_data <- trk_snackf("audio.wav", toFile = FALSE)
#' plot(formant_data$fm_1, type = "l", main = "F1 trajectory")
#'
#' # Track 5 formants
#' trk_snackf("speech.wav", numFormants = 5)
#'
#' # Custom parameters
#' trk_snackf("vowels.wav", numFormants = 3, windowShift = 10, preEmphasis = 0.9)
#'
#' # Batch processing
#' files <- list.files("audio_dir", pattern = ".wav$", full.names = TRUE)
#' trk_snackf(files, verbose = TRUE)
#' }
#'
#' @export
trk_snackf <- function(listOfFiles,
                         beginTime = 0.0,
                         endTime = 0.0,
                         windowShift = 5.0,
                         windowLength = 49.0,
                         numFormants = 4,
                         lpcOrder = NULL,
                         preEmphasis = 0.7,
                         toFile = TRUE,
                         explicitExt = "snackfmt",
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
  
  # Validate formant count
  if (numFormants < 1 || numFormants > 7) {
    cli::cli_abort("numFormants must be between 1 and 7")
  }
  
  # Set default LPC order if not specified
  if (is.null(lpcOrder)) {
    lpcOrder <- 2 * numFormants + 6
  }
  
  # Normalize time parameters
  beginTime <- if (is.null(beginTime)) rep(0.0, n_files) else beginTime
  endTime <- if (is.null(endTime)) rep(0.0, n_files) else endTime
  
  if (length(beginTime) == 1) beginTime <- rep(beginTime, n_files)
  if (length(endTime) == 1) endTime <- rep(endTime, n_files)
  
  if (verbose) {
    cli::cli_inform("Extracting {numFormants} formant{?s} from {cli::no(n_files)} recording{?s}")
  }
  
  # Check that Python script exists
  python_script <- system.file("python", "snack_formant.py", package = "superassp")
  if (!file.exists(python_script)) {
    cli::cli_abort(c(
      "x" = "Python script not found: {.file snack_formant.py}",
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
        numFormants = as.integer(numFormants),
        lpcOrder = as.integer(lpcOrder),
        windowShift = windowShift,
        windowLength = windowLength,
        preEmphasis = preEmphasis,
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
      formant_matrix <- matrix(unlist(result$formants), ncol = result$n_formants, byrow = TRUE)
      bandwidth_matrix <- matrix(unlist(result$bandwidths), ncol = result$n_formants, byrow = TRUE)
      sample_rate <- result$sample_rate
      n_frames <- result$n_frames
      
      # Create AsspDataObj
      sampleRate <- 1000.0 / windowShift
      
      outDataObj <- list()
      
      # Track formats: REAL32 for all formant frequencies and bandwidths
      track_formats <- rep("REAL32", 2 * numFormants)
      attr(outDataObj, "trackFormats") <- track_formats
      attr(outDataObj, "sampleRate") <- sampleRate
      attr(outDataObj, "origFreq") <- sample_rate
      
      startTime <- 1 / sampleRate
      attr(outDataObj, "startTime") <- as.numeric(startTime)
      attr(outDataObj, "startRecord") <- as.integer(1)
      attr(outDataObj, "endRecord") <- as.integer(n_frames)
      
      class(outDataObj) <- "AsspDataObj"
      AsspFileFormat(outDataObj) <- "SSFF"
      AsspDataFormat(outDataObj) <- as.integer(2)
      
      # Add formant frequency tracks
      for (j in seq_len(numFormants)) {
        track_name <- sprintf("fm_%d", j)
        track_data <- matrix(formant_matrix[, j], ncol = 1)
        outDataObj <- wrassp::addTrack(outDataObj, track_name, track_data, "REAL32")
      }
      
      # Add bandwidth tracks
      for (j in seq_len(numFormants)) {
        track_name <- sprintf("bw_%d", j)
        track_data <- matrix(bandwidth_matrix[, j], ncol = 1)
        outDataObj <- wrassp::addTrack(outDataObj, track_name, track_data, "REAL32")
      }
      
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
attr(trk_snackf, "ext") <- "snackfmt"
attr(trk_snackf, "tracks") <- function(n) {
  c(sprintf("fm_%d", seq_len(n)), sprintf("bw_%d", seq_len(n)))
}
attr(trk_snackf, "outputType") <- "SSFF"
