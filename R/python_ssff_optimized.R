#' Compute f0 using the SWIPE algorithm (Optimized)
#' 
#' Optimized version of SWIPE f0 computation following the forest function template
#' with batch processing, efficient file handling, and proper logging.
#' 
#' @inheritParams forest
#' @param minF Candidate f0 frequencies below this frequency will not be considered. 
#' @param maxF Candidates above this frequency will be ignored.
#' @param voicing.threshold Voice/unvoiced threshold. Default is 0.3.
#' 
#' @return If `toFile` is `FALSE`, the function returns a list of AsspDataObj
#'   objects. If `toFile` is `TRUE`, the number (integer) of successfully
#'   processed and stored output files is returned.
#'   
#' @export
swipe_opt <- function(listOfFiles = NULL,
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
  
  #### Fast-path: check if all files are native ####
  listOfFiles <- fast_strip_file_protocol(listOfFiles)
  listOfFiles <- normalizePath(path.expand(listOfFiles))
  
  file_exts <- fast_file_ext(listOfFiles)
  is_native <- fast_is_native(file_exts, nativeFiletypes)
  needs_timewindow <- (beginTime != 0.0 | endTime != 0.0)
  
  # Fast path: all files native, no time windows
  if(all(is_native) && !any(needs_timewindow)) {
    if(verbose) {
      cli::cli_inform("Applying {.fun {funName}} to {cli::no(n_files)} recording{?s}")
    }
    
    listOfFilesDF <- data.frame(
      audio = listOfFiles,
      dsp_input = listOfFiles,
      beginTime = beginTime,
      endTime = endTime,
      stringsAsFactors = FALSE
    )
    toClear <- character(0)
    
  } else {
    # Slow path: needs conversion or time windowing
    listOfFiles_toClear <- convertInputMediaFiles(
      listOfFiles, beginTime, endTime, windowShift,
      nativeFiletypes, preferedFiletype, knownLossless,
      funName, keepConverted, verbose
    )
    
    listOfFilesDF <- listOfFiles_toClear[[1]]
    toClear <- listOfFiles_toClear[[3]]
    
    # Verify all files are in native format
    file_exts <- fast_file_ext(listOfFilesDF$dsp_input)
    if(!all(fast_is_native(file_exts, nativeFiletypes))) {
      cli::cli_abort("File conversion failed - non-native formats remain")
    }
  }
  
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
  
  # Pass parameters to Python
  reticulate::py_run_string("
import numpy as np
import gc
import pysptk as sp
import librosa as lr
")
  
  py <- reticulate::import_main()
  py$soundFile <- soundFile
  py$ws <- windowShift
  py$fMax <- maxF
  py$fMin <- minF
  py$bt <- beginTime
  py$et <- endTime
  py$vt <- voicing.threshold
  
  # Run Python computation
  reticulate::py_run_string("
if et > 0:
    x, fs = lr.load(soundFile, dtype=np.float64, offset=bt, duration=et - bt)
else:
    x, fs = lr.load(soundFile, dtype=np.float64, offset=bt)

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
                          msg = "Invalid AsspDataObj created by swipe_opt")
  
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

attr(swipe_opt, "ext") <- "swi"
attr(swipe_opt, "tracks") <- c("f0", "pitch")
attr(swipe_opt, "outputType") <- "SSFF"
attr(swipe_opt, "nativeFiletypes") <- c("wav")
attr(swipe_opt, "suggestCaching") <- FALSE


#' Compute f0 using the RAPT algorithm (Optimized)
#' 
#' Optimized version of RAPT f0 computation following the forest function template.
#' 
#' @inheritParams swipe_opt
#' 
#' @return If `toFile` is `FALSE`, the function returns a list of AsspDataObj
#'   objects. If `toFile` is `TRUE`, the number (integer) of successfully
#'   processed and stored output files is returned.
#'   
#' @export
rapt_opt <- function(listOfFiles = NULL,
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
  
  listOfFiles <- fast_strip_file_protocol(listOfFiles)
  listOfFiles <- normalizePath(path.expand(listOfFiles))
  
  file_exts <- fast_file_ext(listOfFiles)
  is_native <- fast_is_native(file_exts, nativeFiletypes)
  needs_timewindow <- (beginTime != 0.0 | endTime != 0.0)
  
  # Fast path check
  if(all(is_native) && !any(needs_timewindow)) {
    if(verbose) {
      cli::cli_inform("Applying {.fun {funName}} to {cli::no(n_files)} recording{?s}")
    }
    
    listOfFilesDF <- data.frame(
      audio = listOfFiles,
      dsp_input = listOfFiles,
      beginTime = beginTime,
      endTime = endTime,
      stringsAsFactors = FALSE
    )
    toClear <- character(0)
    
  } else {
    listOfFiles_toClear <- convertInputMediaFiles(
      listOfFiles, beginTime, endTime, windowShift,
      nativeFiletypes, preferedFiletype, knownLossless,
      funName, keepConverted, verbose
    )
    
    listOfFilesDF <- listOfFiles_toClear[[1]]
    toClear <- listOfFiles_toClear[[3]]
    
    file_exts <- fast_file_ext(listOfFilesDF$dsp_input)
    if(!all(fast_is_native(file_exts, nativeFiletypes))) {
      cli::cli_abort("File conversion failed")
    }
  }
  
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
  
  reticulate::py_run_string("
import numpy as np
import gc
import pysptk as sp
import librosa as lr
")
  
  py <- reticulate::import_main()
  py$soundFile <- soundFile
  py$ws <- windowShift
  py$fMax <- maxF
  py$fMin <- minF
  py$bt <- beginTime
  py$et <- endTime
  
  reticulate::py_run_string("
if et > 0:
    x, fs = lr.load(soundFile, dtype=np.float64, offset=bt, duration=et - bt)
else:
    x, fs = lr.load(soundFile, dtype=np.float64, offset=bt)

f0_rapt = sp.rapt(x.astype(np.float64), fs=fs, hopsize=ws / 1000 * fs, min=fMin, max=fMax, otype='f0')
pitch_rapt = sp.rapt(x.astype(np.float64), fs=fs, hopsize=ws / 1000 * fs, min=fMin, max=fMax, otype='pitch')
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
                          msg = "Invalid AsspDataObj created by rapt_opt")
  
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

attr(rapt_opt, "ext") <- "rpt"
attr(rapt_opt, "tracks") <- c("f0", "pitch")
attr(rapt_opt, "outputType") <- "SSFF"
attr(rapt_opt, "nativeFiletypes") <- c("wav")
attr(rapt_opt, "suggestCaching") <- FALSE


#' Extract f0 tracks using the REAPER algorithm (Optimized)
#'
#' Robust Epoch And Pitch EstimatoR (REAPER) algorithm with optimized
#' file processing following the forest function template.
#'
#' @inheritParams swipe_opt
#' @param minF Minimum f0 in Hz
#' @param maxF Maximum f0 in Hz
#' @param unvoiced_cost Cost for unvoiced segments (0-1, higher = more f0)
#' @param high.pass Perform high-pass filtering?
#' @param hilbert.transform Use Hilbert transform for phase distortion?
#'
#' @return If \code{toFile} is \code{FALSE}, returns list of AsspDataObj objects.
#'   If \code{toFile} is \code{TRUE}, returns count of successfully processed files.
#'
#' @export
reaper_opt <- function(listOfFiles = NULL,
                       beginTime = 0.0,
                       endTime = 0.0,
                       windowShift = 5.0,
                       minF = 40,
                       maxF = 500,
                       unvoiced_cost = 0.9,
                       high.pass = TRUE,
                       hilbert.transform = FALSE,
                       explicitExt = "rp0",
                       outputDirectory = NULL,
                       toFile = TRUE,
                       assertLossless = NULL,
                       logToFile = FALSE,
                       keepConverted = FALSE,
                       verbose = TRUE) {
  
  ## Initial constants
  explicitExt <- ifelse(is.null(explicitExt), "rp0", explicitExt)
  newTracknames <- c("f0", "corr")
  nativeFiletypes <- c("wav")
  
  currCall <- rlang::current_call()
  funName <- rlang::call_name(currCall)
  preferedFiletype <- nativeFiletypes[[1]]
  
  knownLossless <- c(assertLossless, knownLossless())
  
  # Normalize time parameters
  beginTime <- if(is.null(beginTime)) 0.0 else beginTime
  endTime <- if(is.null(endTime)) 0.0 else endTime
  
  n_files <- length(listOfFiles)
  
  # Validate and recycle
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
  
  listOfFiles <- fast_strip_file_protocol(listOfFiles)
  listOfFiles <- normalizePath(path.expand(listOfFiles))
  
  file_exts <- fast_file_ext(listOfFiles)
  is_native <- fast_is_native(file_exts, nativeFiletypes)
  needs_timewindow <- (beginTime != 0.0 | endTime != 0.0)
  
  # Fast path
  if(all(is_native) && !any(needs_timewindow)) {
    if(verbose) {
      cli::cli_inform("Applying {.fun {funName}} to {cli::no(n_files)} recording{?s}")
    }
    
    listOfFilesDF <- data.frame(
      audio = listOfFiles,
      dsp_input = listOfFiles,
      beginTime = beginTime,
      endTime = endTime,
      stringsAsFactors = FALSE
    )
    toClear <- character(0)
    
  } else {
    listOfFiles_toClear <- convertInputMediaFiles(
      listOfFiles, beginTime, endTime, windowShift,
      nativeFiletypes, preferedFiletype, knownLossless,
      funName, keepConverted, verbose
    )
    
    listOfFilesDF <- listOfFiles_toClear[[1]]
    toClear <- listOfFiles_toClear[[3]]
    
    file_exts <- fast_file_ext(listOfFilesDF$dsp_input)
    if(!all(fast_is_native(file_exts, nativeFiletypes))) {
      cli::cli_abort("File conversion failed")
    }
  }
  
  #### DSP Processing ####
  if(verbose) {
    cli::cli_inform("Applying {.fun {funName}} to {cli::no(n_files)} recording{?s}")
  }
  
  if(!reticulate::py_module_available("pyreaper")) {
    cli::cli_abort(c(
      "Python module {.pkg pyreaper} is not available.",
      "i" = "Install with: pip install pyreaper"
    ))
  }
  
  externalRes <- vector("list", nrow(listOfFilesDF))
  
  for(i in seq_len(nrow(listOfFilesDF))) {
    origSoundFile <- listOfFilesDF$dsp_input[i]
    bt <- listOfFilesDF$beginTime[i]
    et <- listOfFilesDF$endTime[i]
    
    tryCatch({
      externalRes[[i]] <- process_reaper_single(
        origSoundFile, bt, et, windowShift, minF, maxF,
        unvoiced_cost, high.pass, hilbert.transform,
        explicitExt, outputDirectory, toFile
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

#' Internal function to process single file with REAPER
#' @keywords internal
process_reaper_single <- function(soundFile, beginTime, endTime, windowShift,
                                   minF, maxF, unvoiced_cost, high.pass,
                                   hilbert.transform, explicitExt,
                                   outputDirectory, toFile) {
  
  # Initialize Python modules (reuse if already loaded)
  if(!exists(".py_reaper_initialized", envir = .GlobalEnv)) {
    reticulate::py_run_string("
import numpy as np
import gc
import librosa as lr
import pyreaper
")
    assign(".py_reaper_initialized", TRUE, envir = .GlobalEnv)
  }
  
  py <- reticulate::import_main()
  py$soundFile <- soundFile
  py$ws <- windowShift / 1000  # REAPER takes seconds
  py$fMax <- maxF
  py$fMin <- minF
  py$bt <- beginTime
  py$et <- endTime
  py$uc <- unvoiced_cost
  py$hp <- high.pass
  py$ht <- hilbert.transform
  
  reticulate::py_run_string("
if et > 0:
    x, fs = lr.load(soundFile, dtype=np.float64, offset=bt, duration=et - bt)
else:
    x, fs = lr.load(soundFile, dtype=np.float64, offset=bt)

# Convert to int16 for REAPER
raw_x = x * 2**15
int_x = raw_x.astype(np.int16)

pm_times, pm, f0_times, f0, corr = pyreaper.reaper(
    x=int_x, 
    fs=fs, 
    minf0=fMin, 
    maxf0=fMax, 
    do_high_pass=hp, 
    do_hilbert_transform=ht,
    frame_period=ws, 
    inter_pulse=ws, 
    unvoiced_cost=uc
)

del x, raw_x, int_x
gc.collect()
")
  
  inTable <- data.frame(
    f0 = py$f0,
    corr = py$corr
  )
  
  sampleRate <- 1 / windowShift * 1000
  startTime <- as.numeric(py$f0_times[[1]])
  
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
  
  # Add f0 track
  f0Table <- inTable %>%
    dplyr::select(f0) %>%
    replace(is.na(.), 0) %>%
    dplyr::mutate(dplyr::across(tidyselect::everything(), as.integer))
  
  names(f0Table) <- NULL
  outDataObj <- addTrack(outDataObj, "f0", as.matrix(f0Table[,1]), "INT16")
  
  # Add correlation track
  corrTable <- inTable %>%
    dplyr::select(corr) %>%
    replace(is.na(.), 0) %>%
    dplyr::mutate(dplyr::across(tidyselect::everything(), as.integer))
  
  names(corrTable) <- NULL
  outDataObj <- addTrack(outDataObj, "corr", as.matrix(corrTable[,1]), "INT16")
  
  # Fix missing values at start
  if(startTime > (1/sampleRate)) {
    nr_of_missing_samples <- as.integer(floor(startTime / (1/sampleRate)))
    
    missing_f0_vals <- matrix(0, nrow = nr_of_missing_samples, ncol = ncol(outDataObj$f0))
    missing_corr_vals <- matrix(0, nrow = nr_of_missing_samples, ncol = ncol(outDataObj$corr))
    
    outDataObj$f0 <- rbind(missing_f0_vals, outDataObj$f0)
    outDataObj$corr <- rbind(missing_corr_vals, outDataObj$corr)
    
    attr(outDataObj, "startTime") <- startTime - nr_of_missing_samples * (1/sampleRate)
  }
  
  assertthat::assert_that(is.AsspDataObj(outDataObj),
                          msg = "Invalid AsspDataObj created by reaper_opt")
  
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

attr(reaper_opt, "ext") <- "rp0"
attr(reaper_opt, "tracks") <- c("f0", "corr")
attr(reaper_opt, "outputType") <- "SSFF"
attr(reaper_opt, "nativeFiletypes") <- c("wav")
attr(reaper_opt, "suggestCaching") <- FALSE
