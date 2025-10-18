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
trk_spectral_momentsp <- function(listOfFiles,
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
    result_df <- reticulate::py$trk_spectral_momentsp(
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

attr(trk_spectral_momentsp, "ext") <- c("spm")
attr(trk_spectral_momentsp, "tracks") <- c("cog", "sd", "skewness", "kurtosis")
attr(trk_spectral_momentsp, "outputType") <- c("SSFF")
