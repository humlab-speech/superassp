##' Calculate Linear Prediction smoothed spectrum 
##'
##' Short-term spectral analysis of the signal in `listOfFiles`
##' using the Fast Fourier Transform and linear predictive smoothing. 
##' 
##' The results will be will be written to an SSFF formated file with the base
##' name of the input file and extension *.lps* in a track *LPS[dB]* which contains amplitudes (on a dB scale) of
##' all frequencies in the computed spectrum.
##'
##' @details The function is a re-write of the [wrassp::lpsSpectrum] function, but
##' with media pre-conversion, better checking of preconditions such as the
##' input file existance, structured logging, and the use of a more modern
##' framework for user feedback.
##'
##' The native file type of this function is "wav" files (in "pcm_s16le"
##' format), SUNs "au", NIST, or CSL formats (kay or NSP extension). Input
##' signal conversion, when needed, is done by
##' [libavcodec](https://ffmpeg.org/libavcodec.html) and the excellent [av::av_audio_convert]
##' wrapper function
##'
##' @note
##' This function takes some time to apply but also result in data in a relatively large matrix.
##' It is therefore not usually efficient to store intermediate results in a cache. 
##' However, if the number of signals it will be applied to 
##' is *very* large, then caching of results may be warranted.
##' 
##' @inheritParams cssSpectrum
##' @param resolution = <freq>: set FFT length to the smallest value which
##' results in a frequency resolution of <freq> Hz or better (default: 40.0)
##' @param fftLength = <num>: set FFT length to <num> points (overrules default
##' and 'resolution' option)
##' @param order = <num>: set prediction order to <num> (default: sampling
##' rate in kHz + 3)
##' @param preemphasis = <val>: set pre-emphasis factor to <val> (default:
##' -0.95)
##' @param deemphasize (default: undo spectral tilt due to
##' pre-emphasis used in LP analysis, i.e. TRUE)
##' 
##' @return The number of successfully written files (if `toFile=TRUE`), or a vector of `AsspDataObj` objects (if `toFile=FALSE`).
##'
##' @author Raphael Winkelmann
##' @author Lasse Bombien
##' @author Fredrik Nylén
##'
##' @seealso \code{\link{dftSpectrum}}, \code{\link{cssSpectrum}}, \code{\link{cepstrum}}; 
##' all derived from *libassp* \insertCite{s5h}{superassp} spectrum function.
##' 
##' @useDynLib superassp, .registration = TRUE
##' @examples
##' # get path to audio file
##' path2wav <- list.files(system.file("samples","sustained", package = "superassp"), pattern = glob2rx("a1.wav"), full.names = TRUE)
##'
##' # calculate linear prediction smoothed spectrum 
##' res <- lpsSpectrum(path2wav, toFile=FALSE)
##' resolution <- attr(res,"origFreq") / ncol(res[[1]])
##'
##' # plot spectral values at midpoint of signal
##' plot(y=res[["CSS[dB]"]][400,],
##'     x=seq(1,ncol(res[[1]]),1)* resolution,
##'     type='l',
##'     xlab='Frequency (Hz)',
##'     ylab='Amplitude (dB)')
##' 
##' @export
##' 
lpsSpectrum <- function(listOfFiles = NULL, optLogFilePath = NULL,
                          beginTime = 0.0, centerTime = FALSE,
                          endTime = 0.0, resolution = 40.0,
                          fftLength = 0, windowSize = 20.0,
                          windowShift = 5.0, window = 'BLACKMAN',
                          order = 0, preemphasis = -0.95, 
                          deemphasize = TRUE, 
                        toFile = TRUE,
                        explicitExt = "lps",
                        outputDirectory = NULL,
                        assertLossless = NULL,
                        logToFile = FALSE,
                        keepConverted=FALSE,
                        verbose = TRUE){
  
  ## Initial constants -- specific to this function
  explicitExt <- ifelse(is.null(explicitExt),"lps",explicitExt)
  newTracknames <- "LPS[dB]"  ## Only used if SSFF tracks needs to be renamed from the called function (in C) before returning the SSFF track obj 
  nativeFiletypes <- c("wav","au","kay","nist","nsp")
  
  if(!isAsspWindowType(toupper(window))){
    cli::cli_abort(c("WindowFunction of type {.val {window}} is not supported!",
                     "i"="Accepted window types for routines implemented in *libassp* are {.field {AsspWindowTypes()}}.")
    )
  }
  
  ## Initial constants -- generics
  currCall <- rlang::current_call()
  funName <- rlang::call_name(currCall)
  preferedFiletype <- nativeFiletypes[[1]]
  
  knownLossless <- c(assertLossless,knownLossless()) #Use the user asserted information about lossless encoding, in addition to what is already known by superassp
  
  #Check begin and end times
  if(is.null(beginTime)) beginTime <- 0 # How the C function expects the argument
  if(is.null(endTime)) endTime <- 0 # How the C function expects the argument
  if(length(beginTime) > 1 && length(beginTime) != length(listOfFiles)) cli::cli_abort("The {.par beginTime} argument need to be a vector of the same length as the {.par listOfFiles} argument.")
  if(length(endTime) > 1 && length(endTime) != length(listOfFiles)) cli::cli_abort("The {.par endTime} argument need to be a vector of the same length as the {.par listOfFiles} argument.")
  
  toClear <- c() 
  
  
  #### Setup logging of the function call ####
  makeOutputDirectory(outputDirectory,logToFile, funName)
  
  
  
  #### [*] Input file conversion ####
  
  
  listOfFiles_toClear <- convertInputMediaFiles(listOfFiles,beginTime,endTime,windowShift,nativeFiletypes,preferedFiletype,knownLossless,funName,keepConverted,verbose)
  
  listOfFilesDF <- purrr::pluck(listOfFiles_toClear,1) |>
    dplyr::rename(x=dsp_input) |>
    dplyr::mutate(.beginTime = beginTime,
                  .endTime= endTime) |>
    dplyr::mutate(.beginTime = ifelse(convert_timewindow,0,.beginTime),
                  .endTime = ifelse(convert_timewindow,0,.endTime)) 
  
  toClear <- listOfFiles_toClear[[2]]
  
  
  assertthat::assert_that(all(tools::file_ext(listOfFilesDF$x) %in% nativeFiletypes )) #Make sure that we have a file that may now be handled
  
  #### Application of DSP C function  ####
  
  
  
  
  if(verbose) cli::cli_inform("Applying the {.fun {funName}} DSP function to {cli::no(length(listOfFiles))} speech recording{?s}")
  
  applyC_DSPfunction <- function(x,.beginTime=0,.endTime=0,...){
    assertthat::assert_that(file.exists(x))
    assertthat::assert_that(tools::file_ext(x) %in% nativeFiletypes)
    
    
    
    externalRes = invisible(.External("performAssp", x, 
                                      fname = "spectrum", beginTime = .beginTime, 
                                      centerTime = centerTime, endTime = .endTime, 
                                      spectrumType = 'LPS',
                                      resolution = resolution, 
                                      fftLength = as.integer(fftLength), windowSize = windowSize, 
                                      windowShift = windowShift, window = window, 
                                      effectiveLength = TRUE, 
                                      order = as.integer(order), preemphasis = preemphasis, 
                                      deemphasize = deemphasize, 
                                      toFile = FALSE, explicitExt = explicitExt, 
                                      progressBar = NULL, outputDirectory = outputDirectory,
                                      PACKAGE = "superassp"))
  
    return(externalRes)
  }
  
  ## Prepare for processing: progress bar
  
  
  process_pb <- FALSE
  if(verbose && FALSE){
    process_pb <- list(name="Applying DSP function",
                       format="{cli::pb_extra$currFunName} {cli::pb_bar} {cli::pb_current}/{cli::pb_total}",
                       show_after=1,
                       clear=FALSE,
                       extra=list(currFunName=funName)
    )
  }
  
  ## Process files
  
  externalRes <- purrr::pmap(.l=listOfFilesDF,.f=applyC_DSPfunction)
  #Rename SSFF track names if needed
  if(! is.null(newTracknames)){
    if(length(names(externalRes[[1]])) != length(newTracknames)) cli::cli_abort(c("Wrong number of track names supplied:",
                                                                                  "i"="The track{?s} in the {.cls SSFF} object {?is/are} named {.field {names(externalRes)}}")
    )
    for(i in 1:length(externalRes)){
      names(externalRes[[i]]) <- newTracknames
    }
    
  }
  
  if(toFile){
    toWriteDF <- tibble::tibble(ssffobj=externalRes,filename=listOfFilesDF[["audio"]],ext=explicitExt,outputDirectory=outputDirectory, verbose=verbose)
    filesCreated <- purrr::pwalk(.l=toWriteDF,.f=writeSSFFOutputFile)
    
  }
  
  #Simplify output if just one file is processed 
  if(length(listOfFiles) == 1) externalRes <- purrr::pluck(externalRes,1)
  
  
  #### [*] Cleanup of possibly converted files  ####
  
  cleanupConvertedInputMediaFiles(toClear[["output"]], keepConverted,verbose)
  
  return(externalRes)
}
attr(lpsSpectrum,"ext") <-  "lps" 
attr(lpsSpectrum,"tracks") <-  c("LPS[dB]")
attr(lpsSpectrum,"outputType") <-  "SSFF"
attr(lpsSpectrum,"nativeFiletypes") <-  c("wav","au","kay","nist","nsp")
attr(lpsSpectrum,"suggestCaching") <-  FALSE

### INTERACTIVE TESTING
#
#f <- normalizePath(list.files(file.path("..","inst","samples"),recursive = TRUE,full.names = TRUE))
#f <- f[grepl("*.aiff",f)]

#lpsSpectrum(f,beginTime=1.2, endTime=2.2, toFile=FALSE,keepConverted = FALSE,verbose = TRUE) -> a
#lpsSpectrum(f, toFile=FALSE,keepConverted = FALSE,verbose = TRUE) -> a

#r <- normalizePath(list.files(file.path("..","inst","samples"),recursive = TRUE,full.names = TRUE,pattern = attr(lpsSpectrum,"ext")))
#unlink(r)
#outputDirectory = "/Users/frkkan96/Desktop/output/"
#Error in cssSpectrum(f, toFile = FALSE, keepConverted = FALSE, outputDirectory = "/Users/frkkan96/Desktop/output/",  : 
#object '"/Users/frkkan96/Desktop/output/"' not found
