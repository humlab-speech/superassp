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
trk_pitchp <- function(listOfFiles,
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
                            outputDirectory = NULL,
                            # Backwards compatibility parameters
                            windowShift = NULL,
                            corr.only = NULL) {

  # Handle backwards compatibility: windowShift -> time_step
  if (!is.null(windowShift)) {
    warning("Parameter 'windowShift' is deprecated. Please use 'time_step' instead.", call. = FALSE)
    time_step <- windowShift / 1000  # Convert from ms to seconds
  }

  # Handle backwards compatibility: corr.only -> only_correlation_methods
  if (!is.null(corr.only)) {
    warning("Parameter 'corr.only' is deprecated. Please use 'only_correlation_methods' instead.", call. = FALSE)
    only_correlation_methods <- corr.only
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

  # Source the Python script
  reticulate::source_python(system.file("python", "praat_pitch.py", package = "superassp"))

  outListOfFiles <- c()

  # Process each file
  for(i in 1:nrow(fileBeginEnd)) {
    origSoundFile <- normalizePath(fileBeginEnd[i, "listOfFiles"], mustWork = TRUE)

    bt <- fileBeginEnd[i, "beginTime"]
    et <- fileBeginEnd[i, "endTime"]

    # Load audio using av and convert to parselmouth Sound
    # This uses pure in-memory processing with NO temp files
    sound <- av_load_for_parselmouth(
      file_path = origSoundFile,
      start_time = if (bt > 0) bt else NULL,
      end_time = if (et > 0) et else NULL,
      channels = 1
    )

    # Convert window shape to Python enum
    py_windowShape <- if(windowShape == "Gaussian1") {
      reticulate::py_eval("pm.WindowShape.GAUSSIAN1")
    } else if(windowShape == "Rectangular") {
      reticulate::py_eval("pm.WindowShape.RECTANGULAR")
    } else {
      reticulate::py_eval("pm.WindowShape.GAUSSIAN1")
    }

    # Call Python function with Sound object (in-memory processing)
    result_df <- reticulate::py$trk_pitchp_from_sound(
      sound = sound,
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

attr(trk_pitchp, "ext") <- c("pit")
attr(trk_pitchp, "tracks") <- c("cc", "ac")
attr(trk_pitchp, "outputType") <- c("SSFF")
