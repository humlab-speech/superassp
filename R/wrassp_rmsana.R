##' Short-term Root Mean Square amplitude of signals
##'
##' @description The RMS amplitude is computed for each window of `windowSize`
##'   length in the input signals files listed in `listOfFiles`. Per default, the
##'   RMS values are expressed in decibel (dB) so that they correspond to the
##'   short-term power of the signal. Input signals not in a file format natively
##'   supported will be converted before the autocorrelation functions are
##'   computed. The conversion process will display warnings about input files
##'   that are not in known losslessly encoded formats.
##'
##'   The results will be will be written to an SSFF formated file with the base
##'   name of the input file and extension *.rms* in a track *rms*.
##'
##' @details The function is a re-write of the [wrassp::rmsana] function, but
##'   with media pre-conversion, better checking of preconditions such as the
##'   input file existence, structured logging, and the use of a more modern
##'   framework for user feedback.
##'
##'   The native file type of this function is "wav" files (in "pcm_s16le"
##'   format), SUNs "au", NIST, or CSL formats (kay or NSP extension). Input
##'   signal conversion, when needed, is done by
##'   [libavcodec](https://ffmpeg.org/libavcodec.html) and the excellent [av::av_audio_convert]
##'   wrapper function
##'   
##' @note
##' This function is not considered computationally expensive enough to require caching of 
##' results if applied to many signals. However, if the number of signals it will be applied to 
##' is *very* long, then caching of results may be warranted.
##' 
##'
##' @inheritParams acfana
##' @param linear Should linear RMS values be computed? The default (`FALSE`)
##'   means that the output will be on a logarithmic decibel scale (dB).
##' @useDynLib superassp, .registration = TRUE
##' @seealso [wrassp::rmsana]
##' @seealso [superassp::AsspWindowTypes]
##' @seealso [av::av_audio_convert]
##'
##' @examples
##' # get path to audio file
##'path2wav <- list.files(system.file("samples","sustained", package = "superassp"), pattern = glob2rx("a.wav"), full.names = TRUE)
##'
##'# calculate short-term autocorrelation
##'res <- rmsana(path2wav, toFile=FALSE)
##'
##'
##'# plot rms values
##'plot(seq(0,numRecs.AsspDataObj(res) - 1) / rate.AsspDataObj(res) +
##'       attr(res, 'startTime'),
##'     res$rms,
##'     type='l',
##'     xlab='time (s)',
##'     ylab='RMS energy (dB)')
##'
##' @export
##' 
rmsana <- function(listOfFiles = NULL,
                   beginTime = 0.0,
                   centerTime = FALSE,
                   endTime = 0.0,
                   windowShift = 5.0,
                   windowSize = 20.0,
                   effectiveLength = TRUE,
                   linear = FALSE,
                   window = 'HAMMING',
                   toFile = TRUE,
                   explicitExt = "rms",
                   outputDirectory = NULL,
                   assertLossless = NULL,
                   logToFile = FALSE,
                   convertOverwrites=FALSE,
                   keepConverted=FALSE,
                   verbose = TRUE) {
  
  ## Initial constants
  currCall <- rlang::current_call()
  funName <- rlang::call_name(currCall)
  nativeFiletypes <- c("wav","au","kay","nist","nsp")
  preferedFiletype <- nativeFiletypes[[1]]
  knownLossless <- c(assertLossless,knownLossless()) #Use the user asserted information about lossless encoding, in addition to what is already known by superassp
  
  if(is.null(beginTime)) beginTime <- 0 # How the C function expects the argument
  if(is.null(endTime)) endTime <- 0 # How the C function expects the argument
  
  
  if(!isAsspWindowType(window)){
    cli::cli_abort("WindowFunction of type {.val window} is not supported!")
  }
  
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
    assertthat::assert_that(length(x) == 1)
    
    ret <- invisible(.External("performAssp", x, 
                                      fname = "rmsana", beginTime = .beginTime, 
                                      centerTime = centerTime, endTime = .endTime, 
                                      windowShift = windowShift, windowSize = windowSize, 
                                      effectiveLength = effectiveLength, linear = linear, 
                                      window = window, toFile = toFile, 
                                      explicitExt = explicitExt, 
                                      progressBar = NULL, outputDirectory = outputDirectory,
                                      PACKAGE = "superassp"))
    
    return(ret)
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
  if(toFile){
    externalRes <- purrr::pwalk(.l=listOfFilesDF,.f=applyC_DSPfunction)
    externalRes <- nrow(externalRes) ## TODO: Please add validation that the file actually was created
  }else{
    externalRes <- purrr::pmap(.l=listOfFilesDF,.f=applyC_DSPfunction)
  }
  #Simplify output if just one file is processed 
  if(length(listOfFiles) == 1) externalRes <- purrr::pluck(externalRes,1)
  
  
  #### [*] Cleanup of possibly converted files  ####
  
  cleanupConvertedInputMediaFiles(toClear[["output"]], keepConverted,verbose)
  
  return(externalRes)
}
attr(rmsana,"ext") <-  "rms" 
attr(rmsana,"tracks") <-  c("rms")
attr(rmsana,"outputType") <-  "SSFF"
attr(rmsana,"nativeFiletypes") <-  c("wav","au","kay","nist","nsp")
attr(rmsana,"suggestCaching") <-  FALSE


### INTERACTIVE TESTING
#
#f <- normalizePath(list.files(file.path("..","inst","samples"),recursive = TRUE,full.names = TRUE))
#  f <- f[grepl("*.aiff",f)]
#  
#rmsana(f,toFile=FALSE,keepConverted = FALSE,verbose = TRUE,convertOverwrites=TRUE) -> a 

#r <- normalizePath(list.files(file.path("..","inst","samples"),recursive = TRUE,full.names = TRUE,pattern = attr(rmsana,"ext")))
#unlink(r)
