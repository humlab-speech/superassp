#' Optimized intensity analysis using pladdrr
#'
#' Computes the intensity (loudness) contour of a sound signal using pladdrr's
#' native R bindings to Praat's C library. Audio is loaded using the av package
#' and converted directly to pladdrr Sound objects, eliminating temporary file
#' creation. Supports all media formats (WAV, MP3, MP4, etc.).
#'
#' @param listOfFiles Vector of file paths to audio files
#' @param beginTime Start time in seconds (0 for beginning of file)
#' @param endTime End time in seconds (0 for end of file)
#' @param time_step Time step between frames in seconds (0 for automatic)
#' @param minimal_f0_frequency Minimum pitch frequency in Hz (determines window length)
#' @param subtract_mean Whether to subtract the mean intensity (DC offset correction)
#' @param windowShape Window shape for time windowing
#' @param relativeWidth Relative width for windowing
#' @param toFile Write output to file (TRUE) or return object (FALSE)
#' @param explicitExt File extension for output files
#' @param outputDirectory Output directory path (NULL for same as input)
#' @param verbose Show progress messages (default: TRUE)
#'
#' @return If toFile=TRUE, returns number of files successfully processed.
#'   If toFile=FALSE, returns an AsspDataObj with intensity track.
#'
#' @details
#' This function uses pladdrr (R bindings to Praat's C library) instead of
#' Python's parselmouth. Advantages:
#' \itemize{
#'   \item Pure R/C implementation (no Python dependency)
#'   \item Native R6 object-oriented interface
#'   \item Direct C library access for performance
#'   \item No numpy conversion overhead
#' }
#'
#' The intensity contour represents the loudness of the sound over time,
#' measured in decibels (dB). The minimum pitch parameter determines the
#' effective window length used for the analysis.
#'
#' @examples
#' \dontrun{
#' # Analyze intensity for a single file
#' intensity <- trk_intensityp("speech.wav", toFile = FALSE)
#'
#' # Batch process multiple files
#' files <- c("speech1.wav", "speech2.wav")
#' trk_intensityp(files, toFile = TRUE, outputDirectory = "output/")
#'
#' # With time windowing
#' intensity <- trk_intensityp("speech.wav",
#'                             beginTime = 1.0,
#'                             endTime = 3.0,
#'                             toFile = FALSE)
#' }
#'
#' @export
trk_intensityp <- function(listOfFiles,
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
      end_time = if (et > 0) et else NULL,
      channels = 1
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

    assertthat::assert_that(wrassp::is.AsspDataObj(outDataObj),
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

attr(trk_intensityp, "ext") <- c("int")
attr(trk_intensityp, "tracks") <- c("intensity")
attr(trk_intensityp, "outputType") <- c("SSFF")
attr(trk_intensityp, "nativeFiletypes") <- c("wav")
