#' Sound intensity contour
#'
#' Computes a time-series intensity (dB SPL) contour using Praat's intensity
#' algorithm via pladdrr. Window length is derived from \code{minimal_f0_frequency}
#' to ensure at least one pitch period per frame.
#'
#' @inheritParams trk_acf
#' @param time_step Numeric. Frame shift in seconds; sets output frame rate
#'   (1 / time_step Hz). Set to 0 for Praat's automatic choice. Default 0.
#' @param minimal_f0_frequency Numeric. Minimum expected pitch in Hz; determines
#'   the effective analysis window length (3 / minimal_f0_frequency). Default 50 Hz.
#' @param subtract_mean Logical. Subtract mean intensity to correct DC offset.
#'   Default \code{TRUE}.
#' @param windowShape Character. Window shape applied to the extracted audio segment.
#'   Default \code{"Gaussian1"}.
#' @param relativeWidth Numeric. Relative width of the window. Default 1.0.
#' @param toFile Logical. If \code{TRUE}, write SSFF output files and return the
#'   count written (invisibly). If \code{FALSE}, return an \code{AsspDataObj}.
#'   Default \code{TRUE}.
#' @param explicitExt Character. Output file extension. Default \code{"int"}.
#'
#' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with track:
#'   \describe{
#'     \item{\code{intensity}}{REAL32, dB SPL, n_frames x 1. Sound pressure level
#'       contour. 0 encodes undefined frames.}
#'   }
#'   Frame rate: \code{1 / time_step} Hz (Praat automatic when \code{time_step = 0}).
#'   If \code{toFile = TRUE}: integer count of files written, returned invisibly.
#'
#' @details
#' Uses Praat's intensity-from-sound algorithm. Window length is set automatically
#' to 3 / \code{minimal_f0_frequency}. Increase \code{minimal_f0_frequency} for a
#' shorter, more time-resolved window; decrease for better low-frequency coverage.
#'
#' @examples
#' \dontrun{
#' # Analyze intensity for a single file
#' intensity <- trk_intensity("speech.wav", toFile = FALSE)
#'
#' # Batch process multiple files
#' files <- c("speech1.wav", "speech2.wav")
#' trk_intensity(files, toFile = TRUE, outputDirectory = "output/")
#'
#' # With time windowing
#' intensity <- trk_intensity("speech.wav",
#'                             beginTime = 1.0,
#'                             endTime = 3.0,
#'                             toFile = FALSE)
#' }
#'
#' @export
trk_intensity <- function(listOfFiles,
                           beginTime = 0.0,
                           endTime = 0.0,
                           time_step = 0.0,
                           minimal_f0_frequency = 50.0,
                           subtract_mean = TRUE,
                           windowShape = "Gaussian1",
                           relativeWidth = 1.0,
                           toFile = TRUE,
                           explicitExt = "int",
                           outputDirectory = NULL,
                           verbose = TRUE) {

  # Check pladdrr availability
  if (!pladdrr_available()) {
    cli::cli_abort(c(
      "x" = "pladdrr package not available",
      "i" = "Install with: install_pladdrr()"
    ))
  }

  # Check if multiple files with toFile=FALSE
  if(length(listOfFiles) > 1 & !toFile) {
    stop("length(listOfFiles) is > 1 and toFile=FALSE! toFile=FALSE only permitted for single files.")
  }

  # Create data frame for file processing
  tryCatch({
    fileBeginEnd <- data.frame(
      listOfFiles = listOfFiles,
      beginTime = beginTime,
      endTime = endTime
    )
  }, error = function(e) {
    stop("The beginTime and endTime must either be a single value or the same length as listOfFiles")
  })

  # Check that all files exist
  filesEx <- file.exists(listOfFiles)
  if(!all(filesEx)) {
    filesNotExist <- listOfFiles[!filesEx]
    stop("Unable to find the sound file(s) ", paste(filesNotExist, collapse = ", "))
  }

  outListOfFiles <- c()
  n_files <- nrow(fileBeginEnd)

  # Progress bar setup
  if (verbose && n_files > 1) {
    pb <- utils::txtProgressBar(min = 0, max = n_files, style = 3)
  }

  # Process each file
  for(i in 1:n_files) {
    origSoundFile <- normalizePath(fileBeginEnd[i, "listOfFiles"], mustWork = TRUE)

    bt <- fileBeginEnd[i, "beginTime"]
    et <- fileBeginEnd[i, "endTime"]

    # Load audio and create pladdrr Sound object (in-memory)
    sound <- av_load_for_pladdrr(
      file_path = origSoundFile,
      start_time = if (bt > 0) bt else NULL,
      end_time = if (et > 0) et else NULL
    )

    # Get duration for time bounds handling
    duration <- sound$get_total_duration()

    # Handle time bounds (if full file analysis needed)
    if (bt == 0.0 && et == 0.0) {
      analysis_sound <- sound
    } else {
      actual_end <- if (et > 0) et else duration
      analysis_sound <- sound$extract_part(
        from_time = bt,
        to_time = actual_end,
        window_shape = windowShape,
        relative_width = relativeWidth,
        preserve_times = TRUE
      )
    }

    # Compute intensity using direct API for speed
    intensity_ptr <- pladdrr::to_intensity_direct(
      analysis_sound,
      minimum_pitch = minimal_f0_frequency,
      time_step = time_step,
      subtract_mean = subtract_mean
    )
    intensity <- pladdrr::Intensity(.xptr = intensity_ptr)

    # Extract frame data
    frames_df <- intensity$as_data_frame()
    names(frames_df) <- c("time", "intensity")

    # Create AsspDataObj
    outDataObj <- list()

    # Calculate sample rate from intensity object
    actual_time_step <- intensity$get_sampling_period()
    sampleRate <- 1 / actual_time_step

    attr(outDataObj, "sampleRate") <- sampleRate
    attr(outDataObj, "origFreq") <- as.numeric(16000)

    startTime <- frames_df$time[1]
    attr(outDataObj, "startTime") <- as.numeric(startTime)
    attr(outDataObj, "startRecord") <- as.integer(1)
    attr(outDataObj, "endRecord") <- as.integer(nrow(frames_df))

    class(outDataObj) <- "AsspDataObj"
    AsspFileFormat(outDataObj) <- "SSFF"
    AsspDataFormat(outDataObj) <- as.integer(2)
    attr(outDataObj, "trackFormats") <- character(0)

    # Add intensity track
    intensity_data <- frames_df$intensity
    intensity_data[is.na(intensity_data)] <- 0
    outDataObj <- wrassp::addTrack(outDataObj, "intensity",
                                   as.matrix(intensity_data), "REAL32")
    # Manually fix trackFormats due to wrassp::addTrack bug
    attr(outDataObj, "trackFormats") <- c(attr(outDataObj, "trackFormats"), "REAL32")

    assertthat::assert_that(inherits(outDataObj, "AsspDataObj"),
                            msg = "The AsspDataObj created by intensity analysis is invalid.")

    # Determine output file path
    ssff_file <- sub("\\.[^.]*$", paste0(".", explicitExt), origSoundFile)
    if(!is.null(outputDirectory)) {
      ssff_file <- file.path(outputDirectory, basename(ssff_file))
    }

    attr(outDataObj, "filePath") <- as.character(ssff_file)

    if(toFile) {
      wrassp::write.AsspDataObj(dobj = outDataObj, file = ssff_file)
      outListOfFiles <- c(outListOfFiles, TRUE)
    }

    # Update progress bar
    if (verbose && n_files > 1) {
      utils::setTxtProgressBar(pb, i)
    }
  }

  # Close progress bar
  if (verbose && n_files > 1) {
    close(pb)
  }

  if(toFile) {
    if (verbose) {
      cli::cli_alert_success("Processed {length(outListOfFiles)} file{?s}")
    }
    return(length(outListOfFiles))
  } else {
    return(outDataObj)
  }
}

attr(trk_intensity, "ext") <- c("int")
attr(trk_intensity, "tracks") <- c("intensity")
attr(trk_intensity, "outputType") <- c("SSFF")
attr(trk_intensity, "nativeFiletypes") <- c("wav")
