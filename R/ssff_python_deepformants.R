#' Track formants using DeepFormants deep learning models
#'
#' DeepFormants \insertCite{dissen2017deepformants}{superassp} uses deep neural networks trained on labeled
#' formant data to provide accurate F1-F4 tracking. The algorithm analyzes audio through LPC
#' (Linear Predictive Coding) feature extraction followed by PyTorch-based neural network prediction.
#'
#' This function tracks formants continuously across the entire audio file at **10ms intervals**,
#' returning a time-series of F1, F2, F3, and F4 values.
#'
#' The implementation uses Numba JIT compilation for 2-3x performance improvement over the original code.
#'
#' @param listOfFiles A character vector of file paths to audio files, or an AVAudio S7 object.
#' @param beginTime The start time of the section of the sound file that should be processed (in seconds).
#'   Default is 0.0 (start of file). For tracking, this parameter is used for time windowing via av package.
#' @param endTime The end time of the section of the sound file that should be processed (in seconds).
#'   Default is 0.0 (end of file). For tracking, this parameter is used for time windowing via av package.
#' @param explicitExt The file extension for the output SSFF file. Default is "dfm" (DeepFormants).
#' @param outputDirectory The directory where output files should be written.
#'   If NULL (default), files are written to the same directory as the input files.
#' @param toFile If TRUE (default), write results to SSFF files and return file count.
#'   If FALSE, return AsspDataObj for single file processing.
#' @param verbose If TRUE (default), print progress messages.
#'
#' @return
#'   If \code{toFile = TRUE}: Returns the number of successfully processed files.
#'   If \code{toFile = FALSE}: Returns an AsspDataObj with four tracks:
#'   \itemize{
#'     \item \strong{F1}: First formant frequency in Hz (REAL32 format)
#'     \item \strong{F2}: Second formant frequency in Hz (REAL32 format)
#'     \item \strong{F3}: Third formant frequency in Hz (REAL32 format)
#'     \item \strong{F4}: Fourth formant frequency in Hz (REAL32 format)
#'   }
#'
#' @details
#' \strong{DeepFormants Algorithm:}
#' \enumerate{
#'   \item \strong{Audio Loading}: File is loaded via av package (supports all formats)
#'   \item \strong{LPC Analysis}: Optimized Levinson-Durbin recursion with Numba JIT
#'   \item \strong{Feature Extraction}: LPC coefficients computed at 10ms intervals
#'   \item \strong{Neural Network}: PyTorch RNN predicts F1-F4 from LPC features
#'   \item \strong{Post-processing}: Formants returned as time-aligned tracks
#' }
#'
#' \strong{Specifications:}
#' \itemize{
#'   \item \strong{Frame Rate}: 10ms intervals (100 Hz)
#'   \item \strong{Formants}: F1, F2, F3, F4
#'   \item \strong{Model}: Pre-trained PyTorch RNN
#'   \item \strong{Performance}: ~5 seconds for 2.3s audio (2x real-time)
#'   \item \strong{Optimization}: Numba JIT provides 2-3x speedup
#' }
#'
#' \strong{Installation Requirements:}
#'
#' DeepFormants requires the R \code{torch} package (no Python needed):
#'
#' \code{install.packages("torch"); torch::install_torch()}
#'
#' \strong{Performance Notes:}
#'
#' DeepFormants is slower than traditional formant trackers like Forest but provides:
#' \itemize{
#'   \item Higher accuracy on difficult speech (creaky voice, nasalization)
#'   \item Consistent performance across different speakers
#'   \item Deep learning-based robustness
#' }
#'
#' For a 3-second audio file:
#' \itemize{
#'   \item DeepFormants tracking: ~6-7 seconds
#'   \item Forest (ASSP): ~150ms
#'   \item Trade-off: Accuracy vs. Speed
#' }
#'
#' @examples
#' \dontrun{
#' # Install R torch backend first (one-time)
#' install.packages("torch"); torch::install_torch()
#'
#' # Basic usage - track formants
#' result <- trk_deepformants("speech.wav", toFile = FALSE)
#'
#' # Extract formant tracks
#' f1 <- result$F1  # Hz
#' f2 <- result$F2  # Hz
#' f3 <- result$F3  # Hz
#' f4 <- result$F4  # Hz
#'
#' # Time windowing
#' result <- trk_deepformants("long_recording.wav",
#'                           beginTime = 1.0,
#'                           endTime = 3.0,
#'                           toFile = FALSE)
#'
#' # Batch processing to SSFF files
#' files <- c("file1.wav", "file2.wav", "file3.mp3")
#' trk_deepformants(files, toFile = TRUE)
#'
#' # S7 AVAudio dispatch (automatic)
#' audio <- read_avaudio("speech.wav", sample_rate = 16000)
#' result <- trk_deepformants(audio, toFile = FALSE)
#' }
#'
#' @references
#' Dissen, S., & Keshet, J. (2017). DeepFormants: Deep neural network for formant estimation.
#' Proceedings of Interspeech 2017.
#'
#' @export
trk_deepformants <- function(listOfFiles,
                             beginTime = 0.0,
                             endTime = 0.0,
                             explicitExt = "dfm",
                             outputDirectory = NULL,
                             toFile = TRUE,
                             verbose = TRUE) {

  # Validate inputs
  if (length(listOfFiles) > 1 && !toFile) {
    stop("length(listOfFiles) is > 1 and toFile=FALSE! toFile=FALSE only permitted for single files.",
         call. = FALSE)
  }

  # Check DeepFormants availability
  if (!deepformants_available()) {
    stop(
      "trk_deepformants() requires the R 'torch' package.\n",
      "Install with: install.packages('torch'); torch::install_torch()\n",
      call. = FALSE
    )
  }

  # Create file/time dataframe
  tryCatch({
    fileBeginEnd <- data.frame(
      listOfFiles = listOfFiles,
      beginTime = beginTime,
      endTime = endTime,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    stop("The beginTime and endTime must either be a single value or the same length as listOfFiles",
         call. = FALSE)
  })

  # Check that all files exist
  filesEx <- file.exists(listOfFiles)
  if (!all(filesEx)) {
    filesNotExists <- listOfFiles[!filesEx]
    stop("Unable to find the sound file(s): ", paste(filesNotExists, collapse = ", "),
         call. = FALSE)
  }

  # Vector of successfully processed files
  outListOfFiles <- c()

  # Process each file
  for (i in seq_len(nrow(fileBeginEnd))) {
    origSoundFile <- normalizePath(fileBeginEnd[i, "listOfFiles"], mustWork = TRUE)
    bt <- fileBeginEnd[i, "beginTime"]
    et <- fileBeginEnd[i, "endTime"]

    if (verbose) {
      message(sprintf("Processing file %d/%d: %s", i, nrow(fileBeginEnd), basename(origSoundFile)))
    }

    tryCatch({
      # Get audio info for metadata
      audio_info <- av::av_media_info(origSoundFile)
      sample_rate <- as.numeric(audio_info$audio$sample_rate)

      # Extract features and run LSTM tracker (pure R torch)
      feat_mat <- extract_deepformants_features(
        origSoundFile,
        begin = if (bt > 0) bt else NULL,
        end   = if (et > 0) et else NULL
      )
      formants <- run_deepformants_tracker(feat_mat)  # (n_frames, 4)

      formant_values <- data.frame(
        F1 = formants[, 1L],
        F2 = formants[, 2L],
        F3 = formants[, 3L],
        F4 = formants[, 4L]
      )

      # Calculate SSFF sample rate (frame rate in Hz)
      # DeepFormants uses 10ms frames = 100 Hz
      ssff_sample_rate <- 100.0

      # Create AsspDataObj
      outDataObj <- list()
      attr(outDataObj, "trackFormats") <- c("REAL32", "REAL32", "REAL32", "REAL32")
      attr(outDataObj, "sampleRate") <- as.numeric(ssff_sample_rate)
      attr(outDataObj, "origFreq") <- as.numeric(sample_rate)
      attr(outDataObj, "startTime") <- as.numeric(if (bt > 0) bt else 0.0)
      attr(outDataObj, "startRecord") <- as.integer(1)
      attr(outDataObj, "endRecord") <- as.integer(nrow(formant_values))
      class(outDataObj) <- "AsspDataObj"

      AsspFileFormat(outDataObj) <- "SSFF"
      AsspDataFormat(outDataObj) <- as.integer(2)  # binary

      # Add formant tracks
      outDataObj <- addTrack(outDataObj, "F1",
                            as.matrix(formant_values$F1),
                            "REAL32")
      outDataObj <- addTrack(outDataObj, "F2",
                            as.matrix(formant_values$F2),
                            "REAL32")
      outDataObj <- addTrack(outDataObj, "F3",
                            as.matrix(formant_values$F3),
                            "REAL32")
      outDataObj <- addTrack(outDataObj, "F4",
                            as.matrix(formant_values$F4),
                            "REAL32")

      # Validate AsspDataObj
      if (!is.AsspDataObj(outDataObj)) {
        stop("The AsspDataObj created by trk_deepformants is invalid.", call. = FALSE)
      }

      # Prepare output file path
      ssff_file <- sub("\\.[^.]+$", paste0(".", explicitExt), origSoundFile)
      if (!is.null(outputDirectory)) {
        ssff_file <- file.path(outputDirectory, basename(ssff_file))
      }

      attr(outDataObj, "filePath") <- as.character(ssff_file)

      # Write to file or accumulate
      if (toFile) {
        write.AsspDataObj(dobj = outDataObj, file = ssff_file)
        outListOfFiles <- c(outListOfFiles, TRUE)

        if (verbose) {
          message(sprintf("  → Written to: %s", ssff_file))
        }
      }

    }, error = function(e) {
      warning(sprintf("Error processing %s: %s", basename(origSoundFile), e$message))
      if (toFile) {
        outListOfFiles <- c(outListOfFiles, FALSE)
      }
    })
  }

  # Return based on toFile flag
  if (toFile) {
    return(sum(outListOfFiles))  # Count of successful files
  } else {
    return(outDataObj)
  }
}

# Set function attributes (for compatibility with existing infrastructure)
attr(trk_deepformants, "ext") <- "dfm"
attr(trk_deepformants, "tracks") <- c("F1", "F2", "F3", "F4")
attr(trk_deepformants, "outputType") <- "SSFF"


#' Estimate formants using DeepFormants deep learning models
#'
#' DeepFormants \insertCite{dissen2017deepformants}{superassp} uses deep neural networks trained on labeled
#' formant data to provide accurate F1-F4 estimation. This function estimates formants within a
#' **specific time window**, returning a single set of F1, F2, F3, and F4 values.
#'
#' For continuous tracking across an entire file, use \code{\link{trk_deepformants}} instead.
#'
#' @param listOfFiles A character vector of file paths to audio files, or an AVAudio S7 object.
#' @param beginTime The start time of the time window in seconds. Required for estimation.
#' @param endTime The end time of the time window in seconds. Required for estimation.
#' @param verbose If TRUE (default), print progress messages.
#'
#' @return
#'   A list (or list of lists for multiple files) with the following elements:
#'   \itemize{
#'     \item \strong{F1}: First formant frequency in Hz
#'     \item \strong{F2}: Second formant frequency in Hz
#'     \item \strong{F3}: Third formant frequency in Hz
#'     \item \strong{F4}: Fourth formant frequency in Hz
#'     \item \strong{file}: Original file path
#'     \item \strong{beginTime}: Start time used
#'     \item \strong{endTime}: End time used
#'   }
#'
#' @details
#' \strong{DeepFormants Estimation Algorithm:}
#' \enumerate{
#'   \item \strong{Audio Loading}: Time window is extracted via av package
#'   \item \strong{LPC Analysis}: Optimized Levinson-Durbin recursion
#'   \item \strong{Feature Extraction}: LPC coefficients from the specified window
#'   \item \strong{Neural Network}: PyTorch model predicts F1-F4 from features
#'   \item \strong{Single Estimate}: Returns one value per formant
#' }
#'
#' \strong{Use Cases:}
#' \itemize{
#'   \item Vowel quality analysis at specific time points
#'   \item Formant measurements from labeled TextGrids
#'   \item Batch processing of vowel tokens
#'   \item Comparative phonetic studies
#' }
#'
#' \strong{Performance:}
#' \itemize{
#'   \item Single estimation: ~2 seconds per file
#'   \item Batch processing: Efficient with multiple time windows
#' }
#'
#' \strong{Installation Requirements:}
#'
#' \code{install.packages("torch"); torch::install_torch()}
#'
#' @examples
#' \dontrun{
#' # Install R torch backend first (one-time)
#' install.packages("torch"); torch::install_torch()
#'
#' # Estimate formants in a time window (e.g., vowel midpoint)
#' result <- lst_deepformants("speech.wav",
#'                            beginTime = 1.2,
#'                            endTime = 1.3)
#' print(result$F1)  # 458.46 Hz
#' print(result$F2)  # 1642.70 Hz
#'
#' # Process multiple files with different time windows
#' files <- c("vowel1.wav", "vowel2.wav", "vowel3.wav")
#' begins <- c(1.0, 1.5, 2.0)
#' ends <- c(1.5, 2.0, 2.5)
#'
#' results <- lst_deepformants(files, begins, ends)
#'
#' # Extract to data frame
#' df <- do.call(rbind, lapply(results, as.data.frame))
#' }
#'
#' @references
#' Dissen, S., & Keshet, J. (2017). DeepFormants: Deep neural network for formant estimation.
#' Proceedings of Interspeech 2017.
#'
#' @seealso
#' \code{\link{trk_deepformants}} for continuous formant tracking
#' \code{\link{trk_deepformants}} for continuous formant tracking
#'
#' @export
lst_deepformants <- function(listOfFiles,
                             beginTime,
                             endTime,
                             verbose = TRUE) {

  # Validate that time windows are provided
  if (missing(beginTime) || missing(endTime)) {
    stop("lst_deepformants() requires beginTime and endTime for estimation mode.\n",
         "Use trk_deepformants() for continuous tracking without time windows.",
         call. = FALSE)
  }

  # Check DeepFormants availability
  if (!deepformants_available()) {
    stop(
      "lst_deepformants() requires the R 'torch' package.\n",
      "Install with: install.packages('torch'); torch::install_torch()\n",
      call. = FALSE)
  }

  # Create file/time dataframe
  tryCatch({
    fileBeginEnd <- data.frame(
      listOfFiles = listOfFiles,
      beginTime = beginTime,
      endTime = endTime,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    stop("The beginTime and endTime must either be a single value or the same length as listOfFiles",
         call. = FALSE)
  })

  # Check that all files exist
  filesEx <- file.exists(listOfFiles)
  if (!all(filesEx)) {
    filesNotExists <- listOfFiles[!filesEx]
    stop("Unable to find the sound file(s): ", paste(filesNotExists, collapse = ", "),
         call. = FALSE)
  }

  # Results list
  results <- list()

  # Process each file
  for (i in seq_len(nrow(fileBeginEnd))) {
    origSoundFile <- normalizePath(fileBeginEnd[i, "listOfFiles"], mustWork = TRUE)
    bt <- fileBeginEnd[i, "beginTime"]
    et <- fileBeginEnd[i, "endTime"]

    if (verbose) {
      message(sprintf("Processing file %d/%d: %s (%.3f-%.3fs)",
                     i, nrow(fileBeginEnd), basename(origSoundFile), bt, et))
    }

    tryCatch({
      # Extract features and run MLP estimator (pure R torch)
      feat_mat <- extract_deepformants_features(origSoundFile, begin = bt, end = et)
      formants  <- run_deepformants_estimator(feat_mat)  # length-4 Hz vector

      # Create result list
      result <- list(
        F1 = formants[1L],
        F2 = formants[2L],
        F3 = formants[3L],
        F4 = formants[4L],
        file = origSoundFile,
        beginTime = bt,
        endTime = et
      )

      results[[i]] <- result

      if (verbose) {
        message(sprintf("  → F1: %.0f Hz, F2: %.0f Hz, F3: %.0f Hz, F4: %.0f Hz",
                       result$F1, result$F2, result$F3, result$F4))
      }

    }, error = function(e) {
      warning(sprintf("Error processing %s (%.3f-%.3fs): %s",
                     basename(origSoundFile), bt, et, e$message))
      results[[i]] <- list(
        F1 = NA_real_,
        F2 = NA_real_,
        F3 = NA_real_,
        F4 = NA_real_,
        file = origSoundFile,
        beginTime = bt,
        endTime = et,
        error = e$message
      )
    })
  }

  # Return single result for single file, otherwise list
  if (length(results) == 1) {
    return(results[[1]])
  } else {
    return(results)
  }
}

