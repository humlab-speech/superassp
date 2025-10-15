##' Optimized Praat functions using Python/Parselmouth
##'
##' This file contains optimized versions of Praat-based signal processing functions
##' that use Python's Parselmouth library instead of calling external Praat processes.
##' This approach provides significant performance improvements while maintaining
##' compatibility with the original functions.
##'
##' @name praat_python_optimized
##' @keywords internal
NULL

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
        track_data <- as.numeric(track_data)
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
        track_data <- as.numeric(track_data)
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
        track_data <- as.numeric(track_data)
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


#' Optimized pitch analysis using Python/Parselmouth
#'
#' This is an optimized version of pitch analysis that uses Python's Parselmouth
#' library instead of calling external Praat. It computes F0 using multiple
#' pitch tracking methods (autocorrelation and cross-correlation).
#'
#' @param listOfFiles Vector of file paths to audio files
#' @param beginTime Start time in seconds (0 for beginning of file)
#' @param endTime End time in seconds (0 for end of file)
#' @param time_step Time step between frames in seconds
#' @param window_length Analysis window length in seconds
#' @param minimum_f0 Minimum pitch in Hz
#' @param maximum_f0 Maximum pitch in Hz
#' @param very_accurate Use very accurate pitch tracking
#' @param number_of_candidates Number of pitch candidates to consider
#' @param silence_threshold Silence threshold
#' @param voicing_threshold Voicing threshold
#' @param octave_cost Cost for octave jumps
#' @param octave_jump_cost Cost for octave jumps
#' @param voiced_voiceless_cost Cost for voiced/voiceless transitions
#' @param only_correlation_methods If TRUE, only use cc and ac methods
#' @param minimum_filter_frequency Minimum filter frequency for SPINET (Hz)
#' @param maximum_filter_frequency Maximum filter frequency for SPINET (Hz)
#' @param number_of_filters Number of filters for SPINET
#' @param maximum_frequency_components Maximum frequency components for SHS (Hz)
#' @param maximum_number_of_subharmonics Maximum number of subharmonics for SHS
#' @param compression_factor Compression factor for SHS
#' @param number_of_points_per_octave Number of points per octave for SHS
#' @param windowShape Window shape for time windowing
#' @param relativeWidth Relative width for windowing
#' @param toFile Write output to file (TRUE) or return object (FALSE)
#' @param explicitExt File extension for output files
#' @param outputDirectory Output directory path (NULL for same as input)
#'
#' @return If toFile=TRUE, returns number of files successfully processed.
#'   If toFile=FALSE, returns an AsspDataObj with pitch tracks.
#'
#' @export
praat_pitch <- function(listOfFiles,
                            beginTime = 0.0,
                            endTime = 0.0,
                            time_step = 0.005,
                            window_length = 0.040,
                            minimum_f0 = 75.0,
                            maximum_f0 = 600.0,
                            very_accurate = TRUE,
                            number_of_candidates = 15,
                            silence_threshold = 0.03,
                            voicing_threshold = 0.45,
                            octave_cost = 0.01,
                            octave_jump_cost = 0.35,
                            voiced_voiceless_cost = 0.14,
                            only_correlation_methods = TRUE,
                            minimum_filter_frequency = 70.0,
                            maximum_filter_frequency = 5000.0,
                            number_of_filters = 250,
                            maximum_frequency_components = 1250.0,
                            maximum_number_of_subharmonics = 15,
                            compression_factor = 0.84,
                            number_of_points_per_octave = 48,
                            windowShape = "Gaussian1",
                            relativeWidth = 1.0,
                            toFile = TRUE,
                            explicitExt = "pit",
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

  # Check that all files exist
  filesEx <- file.exists(listOfFiles)
  if(!all(filesEx)) {
    filesNotExist <- listOfFiles[!filesEx]
    stop("Unable to find the sound file(s) ", paste(filesNotExist, collapse = ", "))
  }

  # Source the Python script
  reticulate::source_python(system.file("python", "praat_pitch.py", package = "superassp"))

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
    result_df <- reticulate::py$praat_pitch(
      origSoundFile,
      beginTime = bt,
      endTime = et,
      time_step = time_step,
      window_length = window_length,
      minimum_f0 = minimum_f0,
      maximum_f0 = maximum_f0,
      very_accurate = very_accurate,
      number_of_candidates = number_of_candidates,
      silence_threshold = silence_threshold,
      voicing_threshold = voicing_threshold,
      octave_cost = octave_cost,
      octave_jump_cost = octave_jump_cost,
      voiced_voiceless_cost = voiced_voiceless_cost,
      only_correlation_methods = only_correlation_methods,
      minimum_filter_frequency = minimum_filter_frequency,
      maximum_filter_frequency = maximum_filter_frequency,
      number_of_filters = number_of_filters,
      maximum_frequency_components = maximum_frequency_components,
      maximum_number_of_subharmonics = maximum_number_of_subharmonics,
      compression_factor = compression_factor,
      number_of_points_per_octave = number_of_points_per_octave,
      windowShape = py_windowShape,
      relativeWidth = relativeWidth
    )

    # Create AsspDataObj
    outDataObj <- list()

    sampleRate <- 1 / time_step
    attr(outDataObj, "sampleRate") <- sampleRate
    attr(outDataObj, "origFreq") <- as.numeric(16000)

    startTime <- result_df$time[1]
    attr(outDataObj, "startTime") <- as.numeric(startTime)
    attr(outDataObj, "startRecord") <- as.integer(1)
    attr(outDataObj, "endRecord") <- as.integer(nrow(result_df))

    class(outDataObj) <- "AsspDataObj"
    AsspFileFormat(outDataObj) <- "SSFF"
    AsspDataFormat(outDataObj) <- as.integer(2)
    attr(outDataObj, "trackFormats") <- character(0)  # Initialize empty, will be populated by addTrack

    # Add cross-correlation track
    if("cc" %in% names(result_df)) {
      cc_data <- result_df$cc
      cc_data[is.na(cc_data)] <- 0
      outDataObj <- wrassp::addTrack(outDataObj, "cc",
                                     as.matrix(cc_data), "REAL32")
      # Manually fix trackFormats due to wrassp::addTrack bug
      attr(outDataObj, "trackFormats") <- c(attr(outDataObj, "trackFormats"), "REAL32")
    }

    # Add auto-correlation track
    if("ac" %in% names(result_df)) {
      ac_data <- result_df$ac
      ac_data[is.na(ac_data)] <- 0
      outDataObj <- wrassp::addTrack(outDataObj, "ac",
                                     as.matrix(ac_data), "REAL32")
      # Manually fix trackFormats due to wrassp::addTrack bug
      attr(outDataObj, "trackFormats") <- c(attr(outDataObj, "trackFormats"), "REAL32")
    }

    # Add optional tracks if present
    if("spinet" %in% names(result_df)) {
      spinet_data <- result_df$spinet
      spinet_data[is.na(spinet_data)] <- 0
      outDataObj <- wrassp::addTrack(outDataObj, "spinet",
                                     as.matrix(spinet_data), "REAL32")
      # Manually fix trackFormats due to wrassp::addTrack bug
      attr(outDataObj, "trackFormats") <- c(attr(outDataObj, "trackFormats"), "REAL32")
    }

    if("shs" %in% names(result_df)) {
      shs_data <- result_df$shs
      shs_data[is.na(shs_data)] <- 0
      outDataObj <- wrassp::addTrack(outDataObj, "shs",
                                     as.matrix(shs_data), "REAL32")
      # Manually fix trackFormats due to wrassp::addTrack bug
      attr(outDataObj, "trackFormats") <- c(attr(outDataObj, "trackFormats"), "REAL32")
    }

    assertthat::assert_that(wrassp::is.AsspDataObj(outDataObj),
                            msg = "The AsspDataObj created by the praat_pitch function is invalid.")

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

attr(praat_pitch, "ext") <- c("pit")
attr(praat_pitch, "tracks") <- c("cc", "ac")
attr(praat_pitch, "outputType") <- c("SSFF")


#' Optimized intensity analysis using Python/Parselmouth
#'
#' This is an optimized version of intensity analysis that uses Python's
#' Parselmouth library instead of calling external Praat. It computes the
#' intensity (loudness) contour of a sound signal.
#'
#' @param listOfFiles Vector of file paths to audio files
#' @param beginTime Start time in seconds (0 for beginning of file)
#' @param endTime End time in seconds (0 for end of file)
#' @param time_step Time step between frames in seconds (0 for automatic)
#' @param minimal_f0_frequency Minimum pitch frequency in Hz
#' @param subtract_mean Whether to subtract the mean intensity
#' @param windowShape Window shape for time windowing
#' @param relativeWidth Relative width for windowing
#' @param toFile Write output to file (TRUE) or return object (FALSE)
#' @param explicitExt File extension for output files
#' @param outputDirectory Output directory path (NULL for same as input)
#'
#' @return If toFile=TRUE, returns number of files successfully processed.
#'   If toFile=FALSE, returns an AsspDataObj with intensity track.
#'
#' @export
praat_intensity <- function(listOfFiles,
                                beginTime = 0.0,
                                endTime = 0.0,
                                time_step = 0.0,
                                minimal_f0_frequency = 50.0,
                                subtract_mean = TRUE,
                                windowShape = "Gaussian1",
                                relativeWidth = 1.0,
                                toFile = TRUE,
                                explicitExt = "int",
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

  # Check that all files exist
  filesEx <- file.exists(listOfFiles)
  if(!all(filesEx)) {
    filesNotExist <- listOfFiles[!filesEx]
    stop("Unable to find the sound file(s) ", paste(filesNotExist, collapse = ", "))
  }

  # Source the Python script
  reticulate::source_python(system.file("python", "praat_intensity.py", package = "superassp"))

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
    result_df <- reticulate::py$praat_intensity(
      origSoundFile,
      beginTime = bt,
      endTime = et,
      time_step = time_step,
      minimal_f0_frequency = minimal_f0_frequency,
      subtract_mean = subtract_mean,
      windowShape = py_windowShape,
      relativeWidth = relativeWidth
    )

    # Create AsspDataObj
    outDataObj <- list()

    # Calculate sample rate from time differences
    # Column name has a space: "Time (s)"
    if(nrow(result_df) > 1) {
      time_diff <- result_df[["Time (s)"]][2] - result_df[["Time (s)"]][1]
      sampleRate <- 1 / time_diff
    } else {
      sampleRate <- 1 / 0.01  # Default
    }

    attr(outDataObj, "sampleRate") <- sampleRate
    attr(outDataObj, "origFreq") <- as.numeric(16000)

    startTime <- result_df[["Time (s)"]][1]
    attr(outDataObj, "startTime") <- as.numeric(startTime)
    attr(outDataObj, "startRecord") <- as.integer(1)
    attr(outDataObj, "endRecord") <- as.integer(nrow(result_df))

    class(outDataObj) <- "AsspDataObj"
    AsspFileFormat(outDataObj) <- "SSFF"
    AsspDataFormat(outDataObj) <- as.integer(2)
    attr(outDataObj, "trackFormats") <- character(0)  # Initialize empty, will be populated by addTrack

    # Add intensity track
    # Column name has a space: "Intensity (dB)"
    intensity_data <- result_df[["Intensity (dB)"]]
    intensity_data[is.na(intensity_data)] <- 0
    outDataObj <- wrassp::addTrack(outDataObj, "intensity",
                                   as.matrix(intensity_data), "REAL32")
    # Manually fix trackFormats due to wrassp::addTrack bug
    attr(outDataObj, "trackFormats") <- c(attr(outDataObj, "trackFormats"), "REAL32")

    assertthat::assert_that(wrassp::is.AsspDataObj(outDataObj),
                            msg = "The AsspDataObj created by the praat_intensity function is invalid.")

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

attr(praat_intensity, "ext") <- c("int")
attr(praat_intensity, "tracks") <- c("intensity")
attr(praat_intensity, "outputType") <- c("SSFF")


#' Optimized spectral moments analysis using Python/Parselmouth
#'
#' This is an optimized version that uses Python's Parselmouth library to
#' compute spectral moments (center of gravity, standard deviation, skewness,
#' and kurtosis) of the spectrum.
#'
#' The spectral moments characterize the shape of the spectrum and are useful
#' for acoustic analysis of speech sounds.
#'
#' @param listOfFiles Vector of file paths to audio files
#' @param beginTime Start time in seconds (0 for beginning of file)
#' @param endTime End time in seconds (0 for end of file)
#' @param windowLength Analysis window length in seconds
#' @param maximum_frequency Maximum frequency to analyze in Hz (0 for Nyquist)
#' @param time_step Time step between frames in seconds
#' @param frequency_step Frequency resolution in Hz
#' @param power Power for moment calculation
#' @param windowShape Window shape for time windowing
#' @param relativeWidth Relative width for windowing
#' @param toFile Write output to file (TRUE) or return object (FALSE)
#' @param explicitExt File extension for output files
#' @param outputDirectory Output directory path (NULL for same as input)
#'
#' @return If toFile=TRUE, returns number of files successfully processed.
#'   If toFile=FALSE, returns an AsspDataObj with spectral moment tracks.
#'
#' @export
praat_spectral_moments <- function(listOfFiles,
                                       beginTime = 0.0,
                                       endTime = 0.0,
                                       windowLength = 0.005,
                                       maximum_frequency = 0.0,
                                       time_step = 0.005,
                                       frequency_step = 20.0,
                                       power = 2.0,
                                       windowShape = "Gaussian1",
                                       relativeWidth = 1.0,
                                       toFile = TRUE,
                                       explicitExt = "spm",
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

  # Check that all files exist
  filesEx <- file.exists(listOfFiles)
  if(!all(filesEx)) {
    filesNotExist <- listOfFiles[!filesEx]
    stop("Unable to find the sound file(s) ", paste(filesNotExist, collapse = ", "))
  }

  # Source the Python script
  reticulate::source_python(system.file("python", "praat_spectral_moments.py", package = "superassp"))

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
    result_df <- reticulate::py$praat_spectral_moments(
      origSoundFile,
      beginTime = bt,
      endTime = et,
      windowLength = windowLength,
      maximum_frequency = maximum_frequency,
      time_step = time_step,
      frequency_step = frequency_step,
      power = power,
      windowShape = py_windowShape,
      relativeWidth = relativeWidth
    )

    # Create AsspDataObj
    outDataObj <- list()

    sampleRate <- 1 / time_step
    attr(outDataObj, "sampleRate") <- sampleRate
    attr(outDataObj, "origFreq") <- as.numeric(16000)

    startTime <- result_df$Time[1]
    attr(outDataObj, "startTime") <- as.numeric(startTime)
    attr(outDataObj, "startRecord") <- as.integer(1)
    attr(outDataObj, "endRecord") <- as.integer(nrow(result_df))

    class(outDataObj) <- "AsspDataObj"
    AsspFileFormat(outDataObj) <- "SSFF"
    AsspDataFormat(outDataObj) <- as.integer(2)
    attr(outDataObj, "trackFormats") <- character(0)  # Initialize empty, will be populated by addTrack

    # Add spectral moment tracks
    if("CenterOfGravity" %in% names(result_df)) {
      cog_data <- result_df$CenterOfGravity
      cog_data[is.na(cog_data)] <- 0
      outDataObj <- wrassp::addTrack(outDataObj, "cog",
                                     as.matrix(cog_data), "REAL32")
      # Manually fix trackFormats due to wrassp::addTrack bug
      attr(outDataObj, "trackFormats") <- c(attr(outDataObj, "trackFormats"), "REAL32")
    }

    if("SD" %in% names(result_df)) {
      sd_data <- result_df$SD
      sd_data[is.na(sd_data)] <- 0
      outDataObj <- wrassp::addTrack(outDataObj, "sd",
                                     as.matrix(sd_data), "REAL32")
      # Manually fix trackFormats due to wrassp::addTrack bug
      attr(outDataObj, "trackFormats") <- c(attr(outDataObj, "trackFormats"), "REAL32")
    }

    if("Skewness" %in% names(result_df)) {
      skew_data <- result_df$Skewness
      skew_data[is.na(skew_data)] <- 0
      outDataObj <- wrassp::addTrack(outDataObj, "skewness",
                                     as.matrix(skew_data), "REAL32")
      # Manually fix trackFormats due to wrassp::addTrack bug
      attr(outDataObj, "trackFormats") <- c(attr(outDataObj, "trackFormats"), "REAL32")
    }

    if("Kurtosis" %in% names(result_df)) {
      kurt_data <- result_df$Kurtosis
      kurt_data[is.na(kurt_data)] <- 0
      outDataObj <- wrassp::addTrack(outDataObj, "kurtosis",
                                     as.matrix(kurt_data), "REAL32")
      # Manually fix trackFormats due to wrassp::addTrack bug
      attr(outDataObj, "trackFormats") <- c(attr(outDataObj, "trackFormats"), "REAL32")
    }

    assertthat::assert_that(wrassp::is.AsspDataObj(outDataObj),
                            msg = "The AsspDataObj created by the praat_spectral_moments function is invalid.")

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

attr(praat_spectral_moments, "ext") <- c("spm")
attr(praat_spectral_moments, "tracks") <- c("cog", "sd", "skewness", "kurtosis")
attr(praat_spectral_moments, "outputType") <- c("SSFF")


#' Optimized formant path analysis using Python/Parselmouth (Burg method)
#'
#' This is an optimized version that uses Python's Parselmouth library to
#' compute formant tracks using Praat's FormantPath (Burg) algorithm. This
#' method automatically finds the optimal formant ceiling and tracks formants
#' across time, making it more robust than simple Burg analysis.
#'
#' @param listOfFiles Vector of file paths to audio files
#' @param beginTime Start time in seconds (0 for beginning of file)
#' @param endTime End time in seconds (0 for end of file)
#' @param time_step Time step between frames in seconds
#' @param number_of_formants Maximum number of formants to extract
#' @param maxHzFormant Maximum formant frequency in Hz
#' @param windowLength Analysis window length in seconds
#' @param pre_emphasis Pre-emphasis frequency in Hz
#' @param ceiling_step_size Step size for ceiling variation
#' @param number_of_steps_each_direction Number of ceiling steps in each direction
#' @param track_formants Whether to use formant tracking
#' @param number_of_tracks Number of formant tracks to extract
#' @param reference_F1 Reference frequency for F1 in Hz
#' @param reference_F2 Reference frequency for F2 in Hz
#' @param reference_F3 Reference frequency for F3 in Hz
#' @param reference_F4 Reference frequency for F4 in Hz
#' @param reference_F5 Reference frequency for F5 in Hz
#' @param frequency_cost Weight for frequency deviation in tracking
#' @param bandwidth_cost Weight for bandwidth in tracking
#' @param transition_cost Weight for formant transitions in tracking
#' @param windowShape Window shape for time windowing
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
praat_formantpath_burg <- function(listOfFiles,
                                       beginTime = 0.0,
                                       endTime = 0.0,
                                       time_step = 0.005,
                                       number_of_formants = 5,
                                       maxHzFormant = 5500.0,
                                       windowLength = 0.025,
                                       pre_emphasis = 50.0,
                                       ceiling_step_size = 0.05,
                                       number_of_steps_each_direction = 4,
                                       track_formants = TRUE,
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
                                       explicitExt = "fpb",
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

  # Check that all files exist
  filesEx <- file.exists(listOfFiles)
  if(!all(filesEx)) {
    filesNotExist <- listOfFiles[!filesEx]
    stop("Unable to find the sound file(s) ", paste(filesNotExist, collapse = ", "))
  }

  # Source the Python script
  reticulate::source_python(system.file("python", "praat_formantpath_burg.py", package = "superassp"))

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

    # Pass spectrogram_window_shape as string (Python will convert to enum)
    # Call Python function
    result_df <- reticulate::py$praat_formantpath_burg(
      origSoundFile,
      beginTime = bt,
      endTime = et,
      time_step = time_step,
      number_of_formants = number_of_formants,
      maxHzFormant = maxHzFormant,
      windowLength = windowLength,
      pre_emphasis = pre_emphasis,
      ceiling_step_size = ceiling_step_size,
      number_of_steps_each_direction = number_of_steps_each_direction,
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
      spectrogram_window_shape = spectrogram_window_shape,  # Pass as string
      spectrogram_resolution = spectrogram_resolution
    )

    # Create AsspDataObj
    outDataObj <- list()

    sampleRate <- 1 / time_step
    attr(outDataObj, "sampleRate") <- sampleRate
    attr(outDataObj, "origFreq") <- as.numeric(16000)

    startTime <- result_df$`time(s)`[1]
    attr(outDataObj, "startTime") <- as.numeric(startTime)
    attr(outDataObj, "startRecord") <- as.integer(1)
    attr(outDataObj, "endRecord") <- as.integer(nrow(result_df))

    class(outDataObj) <- "AsspDataObj"
    AsspFileFormat(outDataObj) <- "SSFF"
    AsspDataFormat(outDataObj) <- as.integer(2)
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
        track_data <- as.numeric(track_data)
        track_data[is.na(track_data)] <- 0
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
        track_data <- as.numeric(track_data)
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
        track_data <- as.numeric(track_data)
        track_data[is.na(track_data)] <- 0
        outDataObj <- wrassp::addTrack(outDataObj, paste0("lv", f),
                                       as.matrix(track_data), "REAL32")
        # Manually fix trackFormats due to wrassp::addTrack bug
        attr(outDataObj, "trackFormats") <- c(attr(outDataObj, "trackFormats"), "REAL32")
      }
    }

    assertthat::assert_that(wrassp::is.AsspDataObj(outDataObj),
                            msg = "The AsspDataObj created by the praat_formantpath_burg function is invalid.")

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

attr(praat_formantpath_burg, "ext") <- c("fpb")
attr(praat_formantpath_burg, "tracks") <- c("fm", "bw", "lv")
attr(praat_formantpath_burg, "outputType") <- c("SSFF")
