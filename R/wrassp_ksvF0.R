##' Finds the f0 using the K.Schaefer-Vincent periodicity detection algorithm
##'
##' Applies Schäefer-Vincent periodicity analysis to find f~o~
##' \insertCite{Schäfer-Vincent.1983.10.1159/000261691}{superassp} on 
##' signals listed in `listOfFiles`. Input signals not in a file format natively
##' supported will be converted before the autocorrelation functions are
##' computed. The conversion process will display warnings about input files
##' that are not in known losslessly encoded formats.
##'
##' The results will be will be written to an SSFF formated file with the base
##' name of the input file and extension *.fo* in a track *fo*. 
##'
##' @details The function is a re-write of the [wrassp::acfana] function, but
##' with media pre-conversion, better checking of preconditions such as the
##' input file existance, structured logging, and the use of a more modern
##' framework for user feedback.
##' 
##' Optionally, location and type of the signal extrema on
##' which the F0 data are based, may be stored in a label
##' file. The name of this file will consist of the base
##' name of the f~o~ file and the extension '.prd'. 
##'
##' The native file type of this function is "wav" files (in "pcm_s16le"
##' format), SUNs "au", NIST, or CSL formats (kay or NSP extension). Input
##' signal conversion, when needed, is done by
##' [libavcodec](https://ffmpeg.org/libavcodec.html) and the excellent [av::av_audio_convert]
##' wrapper function
##' 
##' @note
##' This function is not considered computationally expensive enough to require caching of 
##' results if applied to many signals. However, if the number of signals it will be applied to 
##' is *very* long, then caching of results may be warranted.
##' 
##' @inheritParams acfana
##' @param gender = <code>  set gender-specific f~o~ ranges; <code> may be:
##' "f[emale]" (80.0 - 640.0 Hz)
##' "m[ale]" (50.0 - 400.0 Hz)
##' "u[nknown]" (default; 50.0 - 600.0 Hz)
##' @param maxF = <freq>: set maximum F0 value to <freq> Hz (default: 500.0)
##' @param minF = <freq>: set minimum F0 value to <freq> Hz (default: 50.0)
##' @param minAmp = <amp>: set amplitude threshold for voiced samples to <amp> (default: 100)
##' @param maxZCR maximum zero crossing rate in Hz (for voicing detection)
##' 
##' @return nrOfProcessedFiles or if only one file to process return AsspDataObj of that file
##' 
##' @author Raphael Winkelmann
##' @author Lasse Bombien
##' @author Fredrik Nylén
##' 
##' @aliases fo_ana ksvfo
##' 
##' @useDynLib superassp, .registration = TRUE
##' @examples
##' # get path to audio file
##' path2wav <- list.files(system.file("extdata", package = "wrassp"), 
##'                        pattern = glob2rx("*.wav"), 
##'                        full.names = TRUE)[1]
##' 
##' # calculate fundamental frequency contour
##' res <- ksv_fo(path2wav, toFile=FALSE)
##' 
##' # plot the fundamental frequency contour
##' plot(seq(0,numRecs.AsspDataObj(res) - 1) / rate.AsspDataObj(res) +
##'        attr(res, 'startTime'),
##'      res$fo, 
##'      type='l', 
##'      xlab='time (s)', 
##'      ylab='F0 frequency (Hz)')
##'      
##' @export
##' @references
#'  \insertAllCited{}
#'  
fo_ana <- ksv_fo <- ksvfo <- function(listOfFiles = NULL, 
                                           beginTime = 0.0, 
                                           endTime = 0.0, 
                                           windowShift = 5.0, 
                                           gender = 'u',
                                           maxF = 600, 
                                           minF = 50, 
                                           minAmp = 50, 
                                           maxZCR = 3000.0, 
                                           toFile = TRUE,
                                           explicitExt = "fo",
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
  
    externalRes = invisible(.External("performAssp", listOfFiles, 
                                      fname = "f0ana", beginTime = .beginTime, 
                                      endTime = .endTime, windowShift = windowShift, 
                                      gender = gender, maxF = maxF, 
                                      minF = minF, minAmp = minAmp, 
                                      maxZCR = maxZCR, explicitExt = explicitExt, 
                                      toFile = FALSE, progressBar = NULL, 
                                      outputDirectory = outputDirectory, PACKAGE = "superassp"))
    
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
attr(ksv_fo,"ext") <-  "fo" 
attr(ksv_fo,"tracks") <-  c("fo")
attr(ksv_fo,"outputType") <-  "SSFF"
attr(ksv_fo,"nativeFiletypes") <-  c("wav","au","kay","nist","nsp")
attr(ksv_fo,"suggestCaching") <-  FALSE


### INTERACTIVE TESTING
#
#f <- normalizePath(list.files(file.path("..","inst","samples"),recursive = TRUE,full.names = TRUE))
#f <- f[grepl("*.aiff",f)]

#acfana(f,beginTime=1.2, endTime=2.2, toFile=FALSE,keepConverted = FALSE,verbose = TRUE) -> a
#acfana(f, toFile=FALSE,keepConverted = FALSE,verbose = TRUE) -> a

