#' Optimized formant analysis using Python/Parselmouth (Burg method)
#'
#' This is an optimized version of \code{\link{praat_formant_burg}} that uses
#' Python's Parselmouth library instead of calling external Praat. It provides
#' significant performance improvements while maintaining the same output format.
#'
#' The function performs formant analysis using Burg's algorithm and optionally
#' tracks formants across time. It also extracts formant amplitudes from a
#' spectrogram.
#'
#' @param listOfFiles Vector of file paths to audio files
#' @param beginTime Start time in seconds (0 for beginning of file)
#' @param endTime End time in seconds (0 for end of file)
#' @param timeStep Time step between frames in seconds
#' @param number_of_formants Maximum number of formants to extract
#' @param maxHzFormant Maximum formant frequency in Hz
#' @param windowLength Analysis window length in seconds
#' @param pre_emphasis Pre-emphasis frequency in Hz
#' @param track_formants Whether to apply formant tracking
#' @param number_of_tracks Number of formant tracks (used if track_formants=TRUE)
#' @param reference_F1 Reference frequency for F1 in Hz (for tracking)
#' @param reference_F2 Reference frequency for F2 in Hz (for tracking)
#' @param reference_F3 Reference frequency for F3 in Hz (for tracking)
#' @param reference_F4 Reference frequency for F4 in Hz (for tracking)
#' @param reference_F5 Reference frequency for F5 in Hz (for tracking)
#' @param frequency_cost Weight for frequency deviation in tracking
#' @param bandwidth_cost Weight for bandwidth in tracking
#' @param transition_cost Weight for formant transitions in tracking
#' @param windowShape Window shape for time windowing (passed to Python)
#' @param relativeWidth Relative width for windowing
#' @param spectrogram_window_shape Window shape for spectrogram
#' @param spectrogram_resolution Frequency resolution for spectrogram in Hz
#' @param toFile Write output to file (TRUE) or return object (FALSE)
#' @param explicitExt File extension for output files
#' @param outputDirectory Output directory path (NULL for same as input)
#'
#' @return If toFile=TRUE, returns number of files successfully processed.
#'   If toFile=FALSE, returns an AsspDataObj with formant tracks.
#'
#' @export
praat_formant_burg <- function(listOfFiles,
                                    beginTime = 0.0,
                                    endTime = 0.0,
                                    timeStep = 0.005,
                                    number_of_formants = 5,
                                    maxHzFormant = 5500.0,
                                    windowLength = 0.025,
                                    pre_emphasis = 50.0,
                                    track_formants = FALSE,
                                    number_of_tracks = 3,
                                    reference_F1 = 550,
                                    reference_F2 = 1650,
                                    reference_F3 = 2750,
                                    reference_F4 = 3850,
                                    reference_F5 = 4950,
                                    frequency_cost = 1.0,
                                    bandwidth_cost = 1.0,
                                    transition_cost = 1.0,
                                    windowShape = "Gaussian1",
                                    relativeWidth = 1.0,
                                    spectrogram_window_shape = "Gaussian",
                                    spectrogram_resolution = 40.0,
                                    toFile = TRUE,
                                    explicitExt = "pfm",
                                    outputDirectory = NULL) {

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

  # Check that all files exist before we begin
  filesEx <- file.exists(listOfFiles)
  if(!all(filesEx)) {
    filesNotExist <- listOfFiles[!filesEx]
    stop("Unable to find the sound file(s) ", paste(filesNotExist, collapse = ", "))
  }

  # Source the Python script
  reticulate::source_python(system.file("python", "praat_formant_burg.py", package = "superassp"))

  # The empty vector of file names that should be returned
  outListOfFiles <- c()

  # Process each file
  for(i in 1:nrow(fileBeginEnd)) {
    origSoundFile <- normalizePath(fileBeginEnd[i, "listOfFiles"], mustWork = TRUE)

    bt <- fileBeginEnd[i, "beginTime"]
    et <- fileBeginEnd[i, "endTime"]

    # Convert window shape to Python enum
    py_windowShape <- if(windowShape == "Gaussian1") {
      reticulate::py_eval("pm.WindowShape.GAUSSIAN1")
    } else if(windowShape == "Rectangular") {
      reticulate::py_eval("pm.WindowShape.RECTANGULAR")
    } else {
      reticulate::py_eval("pm.WindowShape.GAUSSIAN1")
    }

    # Call Python function
    result_df <- reticulate::py$praat_formant_burg(
      origSoundFile,
      beginTime = bt,
      endTime = et,
      timeStep = timeStep,
      number_of_formants = number_of_formants,
      maxHzFormant = maxHzFormant,
      windowLength = windowLength,
      pre_emphasis = pre_emphasis,
      track_formants = track_formants,
      number_of_tracks = number_of_tracks,
      reference_F1 = reference_F1,
      reference_F2 = reference_F2,
      reference_F3 = reference_F3,
      reference_F4 = reference_F4,
      reference_F5 = reference_F5,
      frequency_cost = frequency_cost,
      bandwidth_cost = bandwidth_cost,
      transition_cost = transition_cost,
      windowShape = py_windowShape,
      relativeWidth = relativeWidth,
      spectrogram_window_shape = spectrogram_window_shape,  # Pass string directly
      spectrogram_resolution = spectrogram_resolution
    )

    # Create AsspDataObj
    outDataObj <- list()

    # Get sample rate from time step
    sampleRate <- 1 / timeStep
    attr(outDataObj, "sampleRate") <- sampleRate

    # Get original frequency (sampling frequency of audio)
    # We need to load this from the file or estimate
    origFreq <- 16000  # Default estimate
    attr(outDataObj, "origFreq") <- as.numeric(origFreq)

    startTime <- result_df$`time(s)`[1]
    attr(outDataObj, "startTime") <- as.numeric(startTime)
    attr(outDataObj, "startRecord") <- as.integer(1)
    attr(outDataObj, "endRecord") <- as.integer(nrow(result_df))

    class(outDataObj) <- "AsspDataObj"
    AsspFileFormat(outDataObj) <- "SSFF"
    AsspDataFormat(outDataObj) <- as.integer(2) # binary
    attr(outDataObj, "trackFormats") <- character(0)  # Initialize empty, will be populated by addTrack

    # Determine how many formants were actually extracted
    actual_formants <- if(track_formants) {
      min(number_of_formants, number_of_tracks)
    } else {
      number_of_formants
    }

    # Add formant frequency tracks (F1, F2, F3, ...)
    for(f in 1:actual_formants) {
      col_name <- paste0("F", f, "(Hz)")
      if(col_name %in% names(result_df)) {
        track_data <- result_df[[col_name]]
        # Convert to numeric (handles character columns with "--undefined--")
        track_data <- suppressWarnings(as.numeric(track_data))
        track_data[is.na(track_data)] <- 0  # Replace NA with 0
        outDataObj <- wrassp::addTrack(outDataObj, paste0("fm", f),
                                       as.matrix(track_data), "REAL32")
        # Manually fix trackFormats due to wrassp::addTrack bug
        attr(outDataObj, "trackFormats") <- c(attr(outDataObj, "trackFormats"), "REAL32")
      }
    }

    # Add bandwidth tracks (B1, B2, B3, ...)
    for(f in 1:actual_formants) {
      col_name <- paste0("B", f, "(Hz)")
      if(col_name %in% names(result_df)) {
        track_data <- result_df[[col_name]]
        # Convert to numeric (handles character columns with "--undefined--")
        track_data <- suppressWarnings(as.numeric(track_data))
        track_data[is.na(track_data)] <- 0
        outDataObj <- wrassp::addTrack(outDataObj, paste0("bw", f),
                                       as.matrix(track_data), "REAL32")
        # Manually fix trackFormats due to wrassp::addTrack bug
        attr(outDataObj, "trackFormats") <- c(attr(outDataObj, "trackFormats"), "REAL32")
      }
    }

    # Add amplitude/level tracks (L1, L2, L3, ...)
    for(f in 1:actual_formants) {
      col_name <- paste0("L", f, "(dB)")
      if(col_name %in% names(result_df)) {
        track_data <- result_df[[col_name]]
        # Convert to numeric (handles character columns with "?" or scientific notation as strings)
        track_data <- suppressWarnings(as.numeric(track_data))
        track_data[is.na(track_data)] <- 0
        outDataObj <- wrassp::addTrack(outDataObj, paste0("lv", f),
                                       as.matrix(track_data), "REAL32")
        # Manually fix trackFormats due to wrassp::addTrack bug
        attr(outDataObj, "trackFormats") <- c(attr(outDataObj, "trackFormats"), "REAL32")
      }
    }

    assertthat::assert_that(wrassp::is.AsspDataObj(outDataObj),
                            msg = "The AsspDataObj created by the praat_formant_burg function is invalid.")

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
  }

  if(toFile) {
    return(length(outListOfFiles))
  } else {
    return(outDataObj)
  }
}

attr(praat_formant_burg, "ext") <- c("pfm")
attr(praat_formant_burg, "tracks") <- c("fm", "bw", "lv")
attr(praat_formant_burg, "outputType") <- c("SSFF")
