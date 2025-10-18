#' Compute f0 using the RAPT algorithm (Non-optimized Python version)
#'
#' @description
#' \code{\link[lifecycle]{lifecycle-superseded}}
#'
#' This function has been superseded by \code{\link{rapt_cpp}}, which provides a native
#' C++ implementation that is 2-3x faster and does not require Python dependencies.
#' This Python-based version is retained for compatibility but is no longer exported.
#'
#' @details
#' This function computes f0 using the Robust Algorithm for Pitch Tracking (RAPT)
#' from the Speech Signal Processing Toolkit (SPTK). It requires Python with pysptk installed.
#'
#' @inheritParams nonopt_swipe
#'
#' @return If `toFile` is `FALSE`, the function returns a list of AsspDataObj
#'   objects. If `toFile` is `TRUE`, the number (integer) of successfully
#'   processed and stored output files is returned.
#'
#' @keywords internal
nonopt_rapt <- function(listOfFiles = NULL,
                     beginTime = 0.0,
                     endTime = 0.0,
                     windowShift = 5.0,
                     minF = 70,
                     maxF = 200,
                     voicing.threshold = 0.3,
                     explicitExt = "rpt",
                     outputDirectory = NULL,
                     toFile = TRUE,
                     assertLossless = NULL,
                     logToFile = FALSE,
                     keepConverted = FALSE,
                     verbose = TRUE) {
  
  ## Initial constants
  explicitExt <- ifelse(is.null(explicitExt), "rpt", explicitExt)
  newTracknames <- c("f0", "pitch")
  nativeFiletypes <- c("wav")
  
  currCall <- rlang::current_call()
  funName <- rlang::call_name(currCall)
  preferedFiletype <- nativeFiletypes[[1]]
  
  knownLossless <- c(assertLossless, knownLossless())
  
  # Normalize time parameters
  beginTime <- if(is.null(beginTime)) 0.0 else beginTime
  endTime <- if(is.null(endTime)) 0.0 else endTime
  
  n_files <- length(listOfFiles)
  
  # Validate and recycle time parameters
  if(length(beginTime) > 1 && length(beginTime) != n_files) {
    cli::cli_abort("The {.par beginTime} must be length 1 or match {.par listOfFiles} length.")
  }
  if(length(endTime) > 1 && length(endTime) != n_files) {
    cli::cli_abort("The {.par endTime} must be length 1 or match {.par listOfFiles} length.")
  }
  
  beginTime <- fast_recycle_times(beginTime, n_files)
  endTime <- fast_recycle_times(endTime, n_files)
  
  #### Setup ####
  makeOutputDirectory(outputDirectory, logToFile, funName)

  #### Memory-based processing: av handles ALL formats directly ####
  listOfFiles <- fast_strip_file_protocol(listOfFiles)
  listOfFiles <- normalizePath(path.expand(listOfFiles))

  # With av_load_for_python(), we don't need file conversion anymore!
  # av can read any format and handle time windows natively
  listOfFilesDF <- data.frame(
    audio = listOfFiles,
    dsp_input = listOfFiles,  # Use original files directly
    beginTime = beginTime,
    endTime = endTime,
    stringsAsFactors = FALSE
  )

  # No files to clean up - everything stays in memory!
  toClear <- character(0)
  
  #### DSP Processing ####
  if(verbose) {
    cli::cli_inform("Applying {.fun {funName}} to {cli::no(n_files)} recording{?s}")
  }
  
  if(!reticulate::py_module_available("pysptk")) {
    cli::cli_abort(c(
      "Python module {.pkg pysptk} is not available.",
      "i" = "Install with: pip install pysptk"
    ))
  }
  
  externalRes <- vector("list", nrow(listOfFilesDF))
  
  for(i in seq_len(nrow(listOfFilesDF))) {
    origSoundFile <- listOfFilesDF$dsp_input[i]
    bt <- listOfFilesDF$beginTime[i]
    et <- listOfFilesDF$endTime[i]
    
    tryCatch({
      externalRes[[i]] <- process_rapt_single(
        origSoundFile, bt, et, windowShift, minF, maxF,
        voicing.threshold, explicitExt, outputDirectory, toFile
      )
    }, error = function(e) {
      cli::cli_abort(c(
        "Error processing {.file {basename(origSoundFile)}}",
        "x" = conditionMessage(e)
      ))
    })
  }
  
  if(!toFile && !is.null(newTracknames)) {
    n_tracks <- length(names(externalRes[[1]]))
    if(n_tracks != length(newTracknames)) {
      cli::cli_abort("Wrong number of track names supplied")
    }
    externalRes <- fast_rename_tracks(externalRes, newTracknames)
  }
  
  if(n_files == 1) externalRes <- externalRes[[1]]
  
  #### Cleanup ####
  cleanupConvertedInputMediaFiles(toClear, keepConverted, verbose)
  
  if(toFile) {
    return(n_files)
  } else {
    return(externalRes)
  }
}

#' Internal function to process single file with RAPT
#' @keywords internal
process_rapt_single <- function(soundFile, beginTime, endTime, windowShift,
                                 minF, maxF, voicing.threshold,
                                 explicitExt, outputDirectory, toFile) {

  # Load audio with av → convert to numpy (MEMORY-BASED, no disk I/O!)
  audio_result <- av_load_for_python(
    soundFile,
    start_time = beginTime,
    end_time = if(endTime == 0) NULL else endTime
  )

  reticulate::py_run_string("
import numpy as np
import gc
import pysptk as sp
")

  py <- reticulate::import_main()
  py$x <- audio_result$audio_np  # Audio already in memory as numpy array!
  py$fs <- audio_result$sample_rate
  py$ws <- windowShift
  py$fMax <- maxF
  py$fMin <- minF

  # Run Python computation (x is already loaded - no librosa.load()!)
  reticulate::py_run_string("
f0_rapt = sp.trk_rapt(x.astype(np.float64), fs=fs, hopsize=ws / 1000 * fs, min=fMin, max=fMax, otype='f0')
pitch_rapt = sp.trk_rapt(x.astype(np.float64), fs=fs, hopsize=ws / 1000 * fs, min=fMin, max=fMax, otype='pitch')
del x
gc.collect()
")
  
  inTable <- data.frame(
    f0 = py$f0_rapt,
    pitch = py$pitch_rapt
  )
  
  sampleRate <- 1 / windowShift * 1000
  startTime <- 1 / sampleRate
  
  outDataObj <- list()
  attr(outDataObj, "trackFormats") <- c("INT16", "INT16")
  attr(outDataObj, "sampleRate") <- sampleRate
  attr(outDataObj, "origFreq") <- as.numeric(py$fs)
  attr(outDataObj, "startTime") <- as.numeric(startTime)
  attr(outDataObj, "startRecord") <- as.integer(1)
  attr(outDataObj, "endRecord") <- as.integer(nrow(inTable))
  class(outDataObj) <- "AsspDataObj"
  
  AsspFileFormat(outDataObj) <- "SSFF"
  AsspDataFormat(outDataObj) <- as.integer(2)
  
  # Add tracks
  f0Table <- inTable %>%
    dplyr::select(f0) %>%
    replace(is.na(.), 0) %>%
    dplyr::mutate(dplyr::across(tidyselect::everything(), as.integer))
  
  names(f0Table) <- NULL
  outDataObj <- addTrack(outDataObj, "f0", as.matrix(f0Table[,1]), "INT16")
  
  pitchTable <- inTable %>%
    dplyr::select(pitch) %>%
    replace(is.na(.), 0) %>%
    dplyr::mutate(dplyr::across(tidyselect::everything(), as.integer))
  
  names(pitchTable) <- NULL
  outDataObj <- addTrack(outDataObj, "pitch", as.matrix(pitchTable[,1]), "INT16")
  
  # Fix missing values
  if(startTime > (1/sampleRate)) {
    nr_of_missing_samples <- as.integer(floor(startTime / (1/sampleRate)))
    
    missing_f0_vals <- matrix(0, nrow = nr_of_missing_samples, ncol = ncol(outDataObj$f0))
    missing_pitch_vals <- matrix(0, nrow = nr_of_missing_samples, ncol = ncol(outDataObj$pitch))
    
    outDataObj$f0 <- rbind(missing_f0_vals, outDataObj$f0)
    outDataObj$pitch <- rbind(missing_pitch_vals, outDataObj$pitch)
    
    attr(outDataObj, "startTime") <- startTime - nr_of_missing_samples * (1/sampleRate)
  }
  
  assertthat::assert_that(is.AsspDataObj(outDataObj),
                          msg = "Invalid AsspDataObj created by rapt")
  
  if(toFile) {
    ssff_file <- sub("wav$", explicitExt, soundFile)
    if(!is.null(outputDirectory)) {
      ssff_file <- file.path(outputDirectory, basename(ssff_file))
    }
    attr(outDataObj, "filePath") <- as.character(ssff_file)
    write.AsspDataObj(dobj = outDataObj, file = ssff_file)
  }
  
  return(outDataObj)
}

attr(nonopt_rapt, "ext") <- "rpt"
attr(nonopt_rapt, "tracks") <- c("f0", "pitch")
attr(nonopt_rapt, "outputType") <- "SSFF"
attr(nonopt_rapt, "nativeFiletypes") <- c("wav")
attr(nonopt_rapt, "suggestCaching") <- FALSE
