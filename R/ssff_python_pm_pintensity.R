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
