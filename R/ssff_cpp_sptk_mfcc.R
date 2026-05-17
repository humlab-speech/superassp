##' Extract Mel-Frequency Cepstral Coefficients (MFCCs) via SPTK
##'
##' Computes HTK-style MFCCs using the SPTK C++ library. MFCCs are the
##' standard frame-level feature for speech recognition, speaker identification,
##' and general audio classification. Covers the full filterbank-DCT pipeline
##' with optional cepstral liftering.
##'
##' @param listOfFiles Character vector of audio file paths. Any format supported by
##'   \pkg{av} is accepted; non-native inputs are transcoded automatically.
##' @param beginTime Numeric. Start of analysis window in seconds. Default 0 (file start).
##' @param endTime Numeric. End of analysis window in seconds. Default 0 (file end).
##' @param windowShift Numeric. Frame shift in milliseconds; sets output frame rate
##'   (1000 / windowShift Hz). Default 10.0 ms.
##' @param windowSize Numeric. Analysis window length in milliseconds. Default 25.0 ms.
##' @param n_mfcc Integer. Number of MFCC coefficients to extract (must be < n_mels).
##'   Default 13.
##' @param n_mels Integer. Number of Mel filterbank channels. Default 40.
##' @param fmin Numeric. Lowest filterbank center frequency in Hz. Default 0.0 Hz.
##' @param fmax Numeric. Highest filterbank center frequency in Hz.
##'   \code{NULL} (default) uses the Nyquist frequency.
##' @param lifter Integer. Cepstral liftering exponent (HTK default is 22). Set to 0
##'   to disable liftering. Default 22.
##' @param floor Numeric. Minimum energy floor for Mel filterbank outputs (prevents
##'   log(0)). Default 1.0.
##' @param toFile Logical. If \code{TRUE}, write SSFF output files and return the
##'   count written invisibly. If \code{FALSE}, return an \code{AsspDataObj}.
##'   Default \code{TRUE}.
##' @param explicitExt Character. Output file extension. Default \code{"mfcc"}.
##' @param outputDirectory Character. Directory for output files. \code{NULL} (default)
##'   writes alongside the input file.
##' @param verbose Logical. Print per-file progress. Default \code{TRUE}.
##'
##' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with tracks:
##'   \describe{
##'     \item{\code{mfcc_0} … \code{mfcc_\{n_mfcc-1\}}}{REAL32, cepstral coefficients
##'       c0 through c\{n_mfcc-1\}, n_frames × 1 each. Dimensionless.}
##'   }
##'   Frame rate: \code{1000 / windowShift} Hz (default 100 Hz).
##'   If \code{toFile = TRUE}: integer count of files written, returned invisibly.
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
trk_mfcc <- function(listOfFiles,
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

  if (verbose) format_apply_msg("trk_mfcc", n_files, beginTime, endTime)

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
      audio_obj <- read_audio(
        file_path,
        begin = bt,
        end   = et
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
        inherits(outDataObj, "AsspDataObj"),
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
