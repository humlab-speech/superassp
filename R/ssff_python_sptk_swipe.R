#' Compute f0 using the SWIPE algorithm (Non-optimized Python version)
#'
#' @description
#' \code{\link[lifecycle]{lifecycle-superseded}}
#'
#' This function has been superseded by \code{\link{swipe_cpp}}, which provides a native
#' C++ implementation that is 2-3x faster and does not require Python dependencies.
#' This Python-based version is retained for compatibility but is no longer exported.
#'
#' 
#' This function takes a sound file and computes f$_0$ and an estimate of pitch 
#' using the Sawtooth Waveform Inspired Pitch Estimator (SWIPE) algorithm \insertCite{Camacho.2008.10.1121/1.2951592}{superassp}. 
#' 
#' @details 
#' The implementation of SWIPE in the Speech Signal Processing Toolkit (SPTK) \insertCite{sptkspeech}{superassp} is used, and called via its Python interface and the [retiulate] R package to compute the signal track.
#' Therefore, the user will have to make sure that a python environment is present and can be attached by the [reticulate]. An anaconda environment is recommended, and can set up by the user by a setup procedure that involve at least these commands:
#' 
#' ```
#' conda create conda create --prefix -n pysuperassp python=3.8 
#' conda activate pysuperassp
#' pip install librosa
#' pip install pysptk
#' ```
#' to make the functionality that this function requires available. 
#'  
#'
#' @param listOfFiles A vector of file paths to wav files.
#' @param beginTime The start time of the section of the sound file that should be processed.
#' @param endTime The end time of the section of the sound file that should be processed.
#' @param windowShift  The measurement interval (frame duration), in seconds.
#' @param minF Candidate f0 frequencies below this frequency will not be considered. 
#' @param maxF Candidates above this frequency will be ignored.
#' @param voicing.threshold Voice/unvoiced threshold. Default is 0.3.
#' @param conda.env The name of the conda environment in which Python and its
#'   required packages are stored. Please make sure that you know what you are
#'   doing if you change this. Defaults to `NULL`, which means that the default enviroment or the environment set in the 
#'   `RETICULATE_PYTHON` environment variable will be used.
#' @inheritParams praat_formant_burg
#'
#' @return
#'  An SSFF track object containing two tracks (f0 and pitch) that are either returned (toFile == FALSE) or stored on disk.
#' 
#' @keywords internal
#'
#' @references
#'   \insertAllCited{}
nonopt_swipe <- function(listOfFiles = NULL,
                      beginTime = 0.0,
                      endTime = 0.0,
                      windowShift = 5.0, 
                      minF = 70, 
                      maxF = 200, 
                      voicing.threshold = 0.3,
                      explicitExt = "swi",
                      outputDirectory = NULL,
                      toFile = TRUE,
                      assertLossless = NULL,
                      logToFile = FALSE,
                      keepConverted = FALSE,
                      verbose = TRUE) {
  
  ## Initial constants -- specific to this function
  explicitExt <- ifelse(is.null(explicitExt), "swi", explicitExt)
  newTracknames <- c("f0", "pitch")
  nativeFiletypes <- c("wav")  # Python DSP only needs wav
  
  ## Initial constants -- generics
  currCall <- rlang::current_call()
  funName <- rlang::call_name(currCall)
  preferedFiletype <- nativeFiletypes[[1]]
  
  knownLossless <- c(assertLossless, knownLossless())
  
  # Normalize time parameters
  beginTime <- if(is.null(beginTime)) 0.0 else beginTime
  endTime <- if(is.null(endTime)) 0.0 else endTime
  
  n_files <- length(listOfFiles)
  
  # Validate time parameter lengths
  if(length(beginTime) > 1 && length(beginTime) != n_files) {
    cli::cli_abort("The {.par beginTime} must be length 1 or match {.par listOfFiles} length.")
  }
  if(length(endTime) > 1 && length(endTime) != n_files) {
    cli::cli_abort("The {.par endTime} must be length 1 or match {.par listOfFiles} length.")
  }
  
  # Use Rcpp for efficient time parameter recycling
  beginTime <- fast_recycle_times(beginTime, n_files)
  endTime <- fast_recycle_times(endTime, n_files)
  
  #### Setup logging ####
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

  #### Application of Python DSP function ####
  if(verbose) {
    cli::cli_inform("Applying {.fun {funName}} to {cli::no(n_files)} recording{?s}")
  }
  
  # Initialize Python environment (once)
  if(!reticulate::py_module_available("pysptk")) {
    cli::cli_abort(c(
      "Python module {.pkg pysptk} is not available.",
      "i" = "Install with: pip install pysptk"
    ))
  }
  
  # Process files with vectorized approach
  externalRes <- vector("list", nrow(listOfFilesDF))
  
  for(i in seq_len(nrow(listOfFilesDF))) {
    origSoundFile <- listOfFilesDF$dsp_input[i]
    bt <- listOfFilesDF$beginTime[i]
    et <- listOfFilesDF$endTime[i]
    
    tryCatch({
      externalRes[[i]] <- process_swipe_single(
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
  
  # Rename tracks if returning data
  if(!toFile && !is.null(newTracknames)) {
    n_tracks <- length(names(externalRes[[1]]))
    if(n_tracks != length(newTracknames)) {
      cli::cli_abort(c(
        "Wrong number of track names supplied:",
        "i" = "Track{?s} named: {.field {names(externalRes[[1]])}}"
      ))
    }
    externalRes <- fast_rename_tracks(externalRes, newTracknames)
  }
  
  # Simplify output for single file
  if(n_files == 1) externalRes <- externalRes[[1]]
  
  #### Cleanup ####
  cleanupConvertedInputMediaFiles(toClear, keepConverted, verbose)
  
  if(toFile) {
    return(n_files)  # Return count of successfully processed files
  } else {
    return(externalRes)
  }
}

#' Internal function to process single file with SWIPE
#'
#' @keywords internal
process_swipe_single <- function(soundFile, beginTime, endTime, windowShift,
                                  minF, maxF, voicing.threshold,
                                  explicitExt, outputDirectory, toFile) {

  # Load audio with av → convert to numpy (MEMORY-BASED, no disk I/O!)
  audio_result <- av_load_for_python(
    soundFile,
    start_time = beginTime,
    end_time = if(endTime == 0) NULL else endTime
  )

  # Pass parameters to Python
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
  py$vt <- voicing.threshold

  # Run Python computation (x is already loaded - no librosa.load()!)
  reticulate::py_run_string("
f0_swipe = sp.swipe(x.astype(np.float64), fs=fs, hopsize=ws / 1000 * fs, min=fMin, max=fMax, otype='f0')
pitch_swipe = sp.swipe(x.astype(np.float64), fs=fs, hopsize=ws / 1000 * fs, min=fMin, max=fMax, otype='pitch')
del x
gc.collect()
")
  
  # Build AsspDataObj
  inTable <- data.frame(
    f0 = py$f0_swipe,
    pitch = py$pitch_swipe
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
  AsspDataFormat(outDataObj) <- as.integer(2)  # binary
  
  # Add f0 track
  f0Table <- inTable %>%
    dplyr::select(f0) %>%
    replace(is.na(.), 0) %>%
    dplyr::mutate(dplyr::across(tidyselect::everything(), as.integer))
  
  names(f0Table) <- NULL
  outDataObj <- addTrack(outDataObj, "f0", as.matrix(f0Table[,1]), "INT16")
  
  # Add pitch track
  pitchTable <- inTable %>%
    dplyr::select(pitch) %>%
    replace(is.na(.), 0) %>%
    dplyr::mutate(dplyr::across(tidyselect::everything(), as.integer))
  
  names(pitchTable) <- NULL
  outDataObj <- addTrack(outDataObj, "pitch", as.matrix(pitchTable[,1]), "INT16")
  
  # Fix missing values at start (Emu-SDMS fix)
  if(startTime > (1/sampleRate)) {
    nr_of_missing_samples <- as.integer(floor(startTime / (1/sampleRate)))
    
    missing_f0_vals <- matrix(0, nrow = nr_of_missing_samples, ncol = ncol(outDataObj$f0))
    missing_pitch_vals <- matrix(0, nrow = nr_of_missing_samples, ncol = ncol(outDataObj$pitch))
    
    outDataObj$f0 <- rbind(missing_f0_vals, outDataObj$f0)
    outDataObj$pitch <- rbind(missing_pitch_vals, outDataObj$pitch)
    
    attr(outDataObj, "startTime") <- startTime - nr_of_missing_samples * (1/sampleRate)
  }
  
  assertthat::assert_that(is.AsspDataObj(outDataObj),
                          msg = "Invalid AsspDataObj created by swipe")
  
  # Write to file if requested
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

attr(nonopt_swipe, "ext") <- "swi"
attr(nonopt_swipe, "tracks") <- c("f0", "pitch")
attr(nonopt_swipe, "outputType") <- "SSFF"
attr(nonopt_swipe, "nativeFiletypes") <- c("wav")
attr(nonopt_swipe, "suggestCaching") <- FALSE


