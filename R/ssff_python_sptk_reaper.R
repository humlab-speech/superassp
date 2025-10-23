

#' Extract f0 tracks using the REAPER algorithm (Non-optimized Python version)
#'
#' @description
#' \code{\link[lifecycle]{lifecycle-superseded}}
#'
#' This function has been superseded by \code{\link{reaper_cpp}}, which provides a native
#' C++ implementation that is 2-3x faster and does not require Python dependencies.
#' This Python-based version is retained for compatibility but is no longer exported.
#'
#' @details
#' Robust Epoch And Pitch EstimatoR (REAPER) algorithm from Google.
#' This implementation requires Python with the pyreaper package installed.
#'
#' @inheritParams trk_swipe
#' @param minF Minimum f0 in Hz
#' @param maxF Maximum f0 in Hz
#' @param unvoiced_cost Cost for unvoiced segments (0-1, higher = more f0)
#' @param high.pass Perform high-pass filtering?
#' @param hilbert.transform Use Hilbert transform for phase distortion?
#'
#' @return If \code{toFile} is \code{FALSE}, returns list of AsspDataObj objects.
#'   If \code{toFile} is \code{TRUE}, returns count of successfully processed files.
#'
#' @keywords internal
nonopt_reaper <- function(listOfFiles = NULL,
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

  # Load audio with av → convert to numpy (MEMORY-BASED, no disk I/O!)
  audio_result <- av_load_for_python(
    soundFile,
    start_time = beginTime,
    end_time = if(endTime == 0) NULL else endTime
  )

  # Initialize Python modules (reuse if already loaded)
  if(!exists(".py_reaper_initialized", envir = .GlobalEnv)) {
    reticulate::py_run_string("
import numpy as np
import gc
import pyreaper
")
    assign(".py_reaper_initialized", TRUE, envir = .GlobalEnv)
  }

  py <- reticulate::import_main()
  py$x <- audio_result$audio_np  # Audio already in memory as numpy array!
  py$fs <- audio_result$sample_rate
  py$ws <- windowShift / 1000  # REAPER takes seconds
  py$fMax <- maxF
  py$fMin <- minF
  py$uc <- unvoiced_cost
  py$hp <- high.pass
  py$ht <- hilbert.transform

  # Run Python computation (x is already loaded - no librosa.load()!)
  reticulate::py_run_string("
# Convert to int16 for REAPER
raw_x = x * 2**15
int_x = raw_x.astype(np.int16)

pm_times, pm, f0_times, f0, corr = pyreaper.trk_reaper(
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
                          msg = "Invalid AsspDataObj created by reaper")
  
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

attr(nonopt_reaper, "ext") <- "rp0"
attr(nonopt_reaper, "tracks") <- c("f0", "corr")
attr(nonopt_reaper, "outputType") <- "SSFF"
attr(nonopt_reaper, "nativeFiletypes") <- c("wav")
attr(nonopt_reaper, "suggestCaching") <- FALSE
