##' Finds the f0 using the K.Schaefer-Vincent periodicity detection algorithm
##'
##' Applies Schäefer-Vincent periodicity analysis \insertCite{Schäfer-Vincent.1983.10.1159/000261691}{superassp} 
##' to find \ifelse{html}{\out{f<sub>o</sub>}}{\eqn{f_o}} (in Hz)
##' along signals listed in `listOfFiles`. Input signals not in a file format natively
##' supported will be converted before the autocorrelation functions are
##' computed. The conversion process will display warnings about input files
##' that are not in known losslessly encoded formats.
##'
##' The results will be will be written to an SSFF formated file with the base
##' name of the input file and extension *.fo* in a track *fo[Hz]*. 
##'
##' @details The function is a re-write of the [wrassp::ksvF0] function, but
##' with media pre-conversion, better checking of preconditions such as the
##' input file existance, structured logging, and the use of a more modern
##' framework for user feedback.
##' 
##' Optionally, location and type of the signal extrema on
##' which the \ifelse{html}{\out{f<sub>o</sub>}}{\eqn{f_o}} data are based, may be stored in a label
##' file. The name of this file will consist of the basename of the `.fo` file and the extension '.prd'. 
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
##' @param gender = <code>  set gender-specific \ifelse{html}{\out{f<sub>o</sub>}}{\eqn{f_o}} ranges; <code> may be:
##' "f[emale]" (80.0 - 640.0 Hz)
##' "m[ale]" (50.0 - 400.0 Hz)
##' "u[nknown]" (default; 50.0 - 600.0 Hz)
##' @param maxF = <freq>: set maximum \ifelse{html}{\out{f<sub>o</sub>}}{\eqn{f_o}} value to <freq> Hz (default: 500.0)
##' @param minF = <freq>: set minimum \ifelse{html}{\out{f<sub>o</sub>}}{\eqn{f_o}} value to <freq> Hz (default: 50.0)
##' @param minAmp = <amp>: set amplitude threshold for voiced samples to <amp> (default: 100)
##' @param maxZCR maximum zero crossing rate in Hz (for voicing detection)
##' 
##' @return The number of successfully written files (if `toFile=TRUE`), or a vector of `AsspDataObj` objects (if `toFile=FALSE`).
##' 
##' @author Raphael Winkelmann
##' @author Lasse Bombien
##' @author Fredrik Nylén
##' 
##' @aliases foana fo_ksv
##' 
##' @references 
##'   \insertAllCited{}
##'
##' @aliases foana fo_ksv
##' 
##' @seealso \code{\link{fo_mhs}} for an alternative pitch tracker
##' @useDynLib superassp, .registration = TRUE
##' @examples
##' # get path to audio file
##'path2wav <- list.files(system.file("samples","sustained", package = "superassp"), pattern = glob2rx("a1.wav"), full.names = TRUE)
##'
##'# calculate fundamental frequency contour
##'res <- ksvfo(path2wav, toFile=FALSE)
##'
##'# plot the fundamental frequency contour
##'plot(seq(0,numRecs.AsspDataObj(res) - 1) / rate.AsspDataObj(res) +
##'       attr(res, 'startTime'),
##'     res[["fo[Hz]"]],
##'     type='l',
##'     xlab='time (s)',
##'     ylab=expression(paste(f[o]," frequency (Hz)")))
##'      
##' @export
ksvfo <- foana <- fo_ksv <- function(listOfFiles = NULL, 
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
  ## Initial constants -- specific to this function
  explicitExt <- ifelse(is.null(explicitExt),"fo",explicitExt)
  newTracknames <- c("fo[Hz]") ## Only used if SSFF tracks needs to be renamed from the called function (in C) before returning the SSFF track obj 
  nativeFiletypes <- c("wav","au","kay","nist","nsp")
  
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
    
  
    externalRes = invisible(.External("performAssp", x, 
                                      fname = "f0ana", beginTime = .beginTime, 
                                      endTime = .endTime, windowShift = windowShift, 
                                      gender = gender, maxF = maxF, 
                                      minF = minF, minAmp = minAmp, 
                                      maxZCR = maxZCR, explicitExt = explicitExt, 
                                      toFile = FALSE, progressBar = FALSE, 
                                      outputDirectory = outputDirectory, PACKAGE = "superassp"))
  
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
attr(ksvfo,"ext") <-  "fo" 
attr(ksvfo,"tracks") <-  c("fo[Hz]")
attr(ksvfo,"outputType") <-  "SSFF"
attr(ksvfo,"nativeFiletypes") <-  c("wav","au","kay","nist","nsp")
attr(ksvfo,"suggestCaching") <-  FALSE

attr(foana,"ext") <-  "fo" 
attr(foana,"tracks") <-  c("fo[Hz]")
attr(foana,"outputType") <-  "SSFF"
attr(foana,"nativeFiletypes") <-  c("wav","au","kay","nist","nsp")
attr(foana,"suggestCaching") <-  FALSE

attr(fo_ksv,"ext") <-  "fo" 
attr(fo_ksv,"tracks") <-  c("fo[Hz]")
attr(fo_ksv,"outputType") <-  "SSFF"
attr(fo_ksv,"nativeFiletypes") <-  c("wav","au","kay","nist","nsp")
attr(fo_ksv,"suggestCaching") <-  FALSE

### INTERACTIVE TESTING
#
#f <- normalizePath(list.files(file.path("..","inst","samples"),recursive = TRUE,full.names = TRUE))
#f <- f[grepl("*.aiff",f)]
#f <- f[1:2]

#fo_ksv(f,beginTime=1.2, endTime=2.2, toFile=FALSE,keepConverted = FALSE,verbose = TRUE) -> a
#fo_ksv(f, toFile=TRUE,keepConverted = FALSE,verbose = TRUE) -> a

