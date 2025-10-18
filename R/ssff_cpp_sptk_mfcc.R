##' SPTK MFCC Extraction (C++ implementation)
##'
##' @description Extract Mel-Frequency Cepstral Coefficients (MFCCs) using the
##'   SPTK library. This is a high-performance C++ implementation that is
##'   significantly faster than Python-based implementations and requires
##'   no Python dependencies.
##'
##'   MFCCs are widely used acoustic features in speech recognition, speaker
##'   identification, and audio analysis. The implementation follows the
##'   standard HTK-style MFCC pipeline with mel-scale filterbanks and
##'   discrete cosine transform.
##'
##'   All input media formats are supported via the av package, including video
##'   files from which audio will be automatically extracted.
##'
##' @param listOfFiles Vector of file paths to process
##' @param beginTime Start time in seconds (default: 0.0)
##' @param endTime End time in seconds (default: 0.0 = end of file)
##' @param windowShift Frame shift in milliseconds (default: 10.0)
##' @param windowSize Window size in milliseconds (default: 25.0)
##' @param n_mfcc Number of MFCC coefficients to extract (default: 13)
##' @param n_mels Number of mel filterbanks (default: 40)
##' @param fmin Minimum frequency in Hz (default: 0.0)
##' @param fmax Maximum frequency in Hz (default: NULL = sample_rate/2)
##' @param lifter Liftering coefficient (default: 22)
##' @param floor Floor value for mel filterbank output (default: 1.0)
##' @param toFile Write results to file (default: TRUE)
##' @param explicitExt Output file extension (default: "mfcc")
##' @param outputDirectory Output directory (default: NULL = same as input)
##' @param verbose Show progress messages (default: TRUE)
##'
##' @return If toFile=TRUE, returns the number of successfully processed files.
##'   If toFile=FALSE, returns AsspDataObj or list of AsspDataObj objects
##'   with MFCC tracks (mfcc_0, mfcc_1, ..., mfcc_n).
##'
##' @note This function uses the SPTK C++ library for high-performance MFCC
##'   extraction. It is 2-3x faster than the Python torchaudio implementation.
##'
##' @seealso [mfcc()] for the deprecated torchaudio implementation
##'
##' @export
##' @examples
##' \dontrun{
##' # Extract 13 MFCCs (default)
##' trk_mfcc("recording.wav")
##'
##' # Extract 20 MFCCs with custom parameters
##' trk_mfcc("speech.mp3", n_mfcc = 20, n_mels = 80)
##'
##' # Return data without writing file
##' mfcc_data <- trk_mfcc("audio.wav", toFile = FALSE)
##'
##' # Process with specific frequency range
##' trk_mfcc("recording.wav", fmin = 80, fmax = 8000)
##'
##' # Process video file (extracts audio)
##' trk_mfcc("interview.mp4")
##' }
trk_mfcc <- function(listOfFiles = NULL,
                      beginTime = 0.0,
                      endTime = 0.0,
                      windowShift = 10.0,
                      windowSize = 25.0,
                      n_mfcc = 13,
                      n_mels = 40,
                      fmin = 0.0,
                      fmax = NULL,
                      lifter = 22,
                      floor = 1.0,
                      toFile = TRUE,
                      explicitExt = "mfcc",
                      outputDirectory = NULL,
                      verbose = TRUE) {

  # Validate inputs
  if (is.null(listOfFiles) || length(listOfFiles) == 0) {
    cli::cli_abort("No input files specified in {.arg listOfFiles}")
  }
  
  if (n_mfcc < 1 || n_mfcc >= n_mels) {
    cli::cli_abort("{.arg n_mfcc} must be positive and less than {.arg n_mels}")
  }
  
  if (n_mels < 1) {
    cli::cli_abort("{.arg n_mels} must be positive")
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
  makeOutputDirectory(outputDirectory, FALSE, "trk_mfcc")

  if (verbose) {
    cli::cli_inform("Extracting {n_mfcc} MFCCs (SPTK C++) from {cli::no(n_files)} recording{?s}")
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
      # Load audio using av package
      audio_obj <- av_to_asspDataObj(
        file_path,
        start_time = bt,
        end_time = if (et > 0) et else NULL
      )
      
      sample_rate <- attr(audio_obj, "sampleRate")
      
      # Set fmax to Nyquist if not specified
      fmax_use <- if (is.null(fmax)) sample_rate / 2 else fmax
      
      # Call C++ function
      result <- sptk_mfcc_cpp(
        audio_obj,
        n_mfcc = as.integer(n_mfcc),
        n_mels = as.integer(n_mels),
        windowShift = windowShift,
        windowSize = windowSize,
        fmin = fmin,
        fmax = fmax_use,
        lifter = as.integer(lifter),
        floor = floor,
        verbose = FALSE
      )
      
      # Create AsspDataObj with results
      outDataObj <- list()
      
      # Calculate frame rate
      frame_rate <- 1000.0 / windowShift
      attr(outDataObj, "sampleRate") <- frame_rate
      attr(outDataObj, "origFreq") <- sample_rate
      
      # Set timing attributes
      start_time <- result$times[1]
      attr(outDataObj, "startTime") <- as.numeric(start_time)
      attr(outDataObj, "startRecord") <- as.integer(1)
      attr(outDataObj, "endRecord") <- as.integer(result$n_frames)
      
      class(outDataObj) <- "AsspDataObj"
      AsspFileFormat(outDataObj) <- "SSFF"
      AsspDataFormat(outDataObj) <- as.integer(2)
      attr(outDataObj, "trackFormats") <- character(0)
      
      # Add MFCC tracks (including c0)
      mfcc_matrix <- result$mfcc
      n_coef <- result$n_mfcc
      
      for (coef in seq_len(n_coef)) {
        track_name <- sprintf("mfcc_%d", coef - 1)  # 0-indexed: c0, c1, c2, ...
        track_data <- matrix(mfcc_matrix[, coef], ncol = 1)
        outDataObj <- wrassp::addTrack(outDataObj, track_name, track_data, "REAL32")
        attr(outDataObj, "trackFormats") <- c(attr(outDataObj, "trackFormats"), "REAL32")
      }
      
      assertthat::assert_that(
        wrassp::is.AsspDataObj(outDataObj),
        msg = "The AsspDataObj created by sptk_mfcc is invalid."
      )
      
      # Handle output
      if (toFile) {
        # Generate output filename
        base_name <- tools::file_path_sans_ext(basename(file_path))
        ssff_file <- paste0(base_name, ".", explicitExt)
        
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

attr(trk_mfcc, "ext") <- "mfcc"
attr(trk_mfcc, "tracks") <- function(n_mfcc = 13) sprintf("mfcc_%d", seq(0, n_mfcc))
attr(trk_mfcc, "outputType") <- "SSFF"
attr(trk_mfcc, "nativeFiletypes") <- c("wav", "mp3", "flac", "ogg", "mp4", "mkv", "avi")
attr(trk_mfcc, "suggestCaching") <- TRUE
