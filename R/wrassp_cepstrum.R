##' Performs a short-term cepstral analysis of the signal
##'
##' @description 
##' This function performs a cepstral analysis of the signal in  `listOfFiles`
##' using the Fast Fourier Transform. The number of
##' coefficients per output record will also equal the
##' FFT length / 2 + 1 (which means that it will not be mirrored).
##' 
##' 
##' The results will be will be written to an SSFF formated file with the base
##' name of the input file and extension *.cep* in a track *Cepstrum* which contains amplitudes of
##' all Quefrencies (in ms).
##'
##' @details The function is a re-write of the [wrassp::cepstrum] function, but
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
##' @inheritParams acfana
##' @param resolution = <freq>: set FFT length to the smallest value which
##' results in a frequency resolution of <freq> Hz or better (default: 40.0)
##' @param fftLength = <num>: set FFT length to <num> points (overrules default
##' and 'resolution' option)
##' 
##' @author Raphael Winkelmann
##' @author Lasse Bombien
##' @author Fredrik Nylén
##' 
##' @seealso \code{\link{dftSpectrum}}, \code{\link{cssSpectrum}}, \code{\link{lpsSpectrum}}; 
##' all derived from libassp's spectrum function
##' 
##' @return The number of successfully written files (if `toFile=TRUE`), or a vector of `AsspDataObj` objects (if `toFile=FALSE`).
##'
##' @seealso [wrassp::cepstrum]
##' @seealso [superassp::AsspWindowTypes]
##' @seealso [av::av_audio_convert]
##' 
##' @useDynLib superassp, .registration = TRUE
##' @examples
##' # get path to audio file
##' path2wav <- list.files(system.file("samples","sustained", package = "superassp"), pattern = glob2rx("a1.wav"), full.names = TRUE)
##'
##' # calulate cepstrum
##' res <- cepstrum(path2wav, toFile=FALSE)
##'
##' # plot cepstral values at midpoint of signal
##' plot(res$Cepstrum[dim(res$Cepstrum)[1]/2,],
##'     type='l',
##'     xlab='Cepstral value index',
##'     ylab='Cepstral value')
##'
##' @export
cepstrum<- function(listOfFiles = NULL, optLogFilePath = NULL,
                       beginTime = 0.0, centerTime = FALSE,
                       endTime = 0.0, resolution = 40.0,
                       fftLength = 0, windowShift = 5.0,
                       window = 'BLACKMAN', 
                    toFile = TRUE,
                    explicitExt = "cep",
                    outputDirectory = NULL,
                    assertLossless = NULL,
                    logToFile = FALSE,
                    keepConverted=FALSE,
                    verbose = TRUE){
  
  ## Initial constants -- specific to this function
  explicitExt <- ifelse(is.null(explicitExt),"cep",explicitExt)
  newTracknames <- "Cepstrum"  ## Only used if SSFF tracks needs to be renamed from the called function (in C) before returning the SSFF track obj 
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
                                      spectrumType = 'CEP',
                                      resolution = resolution, 
                                      fftLength = as.integer(fftLength),
                                      windowShift = windowShift, window = window, 
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
attr(cepstrum,"ext") <-  "cep" 
attr(cepstrum,"tracks") <-  c("Cepstrum")
attr(cepstrum,"outputType") <-  "SSFF"
attr(cepstrum,"nativeFiletypes") <-  c("wav","au","kay","nist","nsp")
attr(cepstrum,"suggestCaching") <-  FALSE

### INTERACTIVE TESTING
#
#f <- normalizePath(list.files(file.path("..","inst","samples"),recursive = TRUE,full.names = TRUE))
#f <- f[grepl("*.aiff",f)]

#acfana(f,beginTime=1.2, endTime=2.2, toFile=FALSE,keepConverted = FALSE,verbose = TRUE) -> a
#cepstrum(f, toFile=FALSE,keepConverted = FALSE,verbose = TRUE) -> a

#r <- normalizePath(list.files(file.path("..","inst","samples"),recursive = TRUE,full.names = TRUE,pattern = attr(acfana,"ext")))
#unlink(r)
#outputDirectory = "/Users/frkkan96/Desktop/output/"
#Error in acfana(f, toFile = FALSE, keepConverted = FALSE, outputDirectory = "/Users/frkkan96/Desktop/output/",  : 
#object '"/Users/frkkan96/Desktop/output/"' not found
