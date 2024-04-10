##' Pitch analysis using a Modified Harmonic Sieve
##'
##' This function finds the pitch (in Hz) along files in <listOfFile> using
##' Michel Scheffers' Modified Harmonic Sieve algorithm implmented in *libassp* \insertCite{s5h}{superassp}. 
##' Input signals not in a file format natively
##' supported will be converted before the autocorrelation functions are
##' computed. The conversion process will display warnings about input files
##' that are not in known losslessly encoded formats.
##'
##' The results will be will be written to an SSFF formated file with the base
##' name of the input file and extension *.pit* in a track *pitch[Hz]*. 
##'
##' @details The function is a re-write of the [wrassp::mhsF0] function, but
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
##' This function is not considered computationally expensive enough to require caching of 
##' results if applied to many signals. However, if the number of signals it will be applied to 
##' is *very* large, then caching of results may be warranted.
##' 
##' @inheritParams ksvfo
##' @param minAmp = <amp>:  minimum signal amplitude (default: 50)
##' @param minAC1 = <freq>: minimum 1st correlation coefficient (default: 0.250)
##' @param minRMS = <num>:  minimum RMS amplitude in dB (default: 18.0)
##' @param maxZCR = <freq>: maximum zero crossing rate in Hz (default: 3000)
##' @param minProb = <num>: minimum quality value of \ifelse{html}{\out{f<sub>o</sub>}}{\eqn{f_o}} fit (default: 0.520)
##' 
##' @return The number of successfully written files (if `toFile=TRUE`), or a vector of `AsspDataObj` objects (if `toFile=FALSE`).
##' 
##' @author Raphael Winkelmann
##' @author Lasse Bombien
##' @author Fredrik Nylén
##' 
##' @aliases mhspitch pitch_mhs
##' 
##' @seealso \code{\link{ksv_fo}} for an alternative algorithm for tracking the fundamental frequency \ifelse{html}{\out{f<sub>o</sub>}}{\eqn{f_o}}.
##' 
##' @useDynLib superassp, .registration = TRUE
##' @references 
##'   \insertAllCited{}
##'
##' @examples
##' # get path to audio file
##' path2wav <- list.files(system.file("samples","sustained", package = "superassp"), pattern = glob2rx("a1.wav"), full.names = TRUE)
##' 
##' # calculate short-term autocorrelation
##' res <- mhspitch(path2wav, toFile=FALSE)
##' 
##' # plot fundamental frequency contour
##' plot(seq(0,numRecs.AsspDataObj(res) - 1) / rate.AsspDataObj(res) +
##'        attr(res, 'startTime'),
##'      res[["pitch[Hz]"]],
##'      type='l',
##'      xlab='time (s)',
##'      ylab="Pitch (Hz)")
##' 
##' @export
'mhspitch' <- 'pitch_mhs' <-function(listOfFiles = NULL,
                                             beginTime = 0.0, 
                                             centerTime = FALSE, 
                                             endTime = 0.0, 
                                             windowShift = 5.0, 
                                             gender = 'u', 
                                             maxF = 600.0, 
                                             minF = 50.0, 
                                             minAmp = 50.0, 
                                             minAC1 = 0.25, 
                                             minRMS = 18.0, 
                                             maxZCR = 3000.0, 
                                             minProb = 0.52, 
                                             plainSpectrum = FALSE, 
                                             toFile = TRUE,
                                             explicitExt = "pit",
                                             outputDirectory = NULL,
                                             assertLossless = NULL,
                                             logToFile = FALSE,
                                             convertOverwrites=FALSE,
                                             keepConverted=FALSE,
                                             verbose = TRUE){
  
  ## Initial constants -- specific to this function
  explicitExt <- ifelse(is.null(explicitExt),"pit",explicitExt)
  newTracknames <- c("pitch[Hz]") ## Only used if SSFF tracks needs to be renamed from the called function (in C) before returning the SSFF track obj 
  nativeFiletypes <- c("wav","au","kay","nist","nsp")
  
  if(! gender %in% c('u','f','m')) cli::cli_abort(c("Incorrect specification of gender: {.val {gender}}",
                                                  "i"="Valid  specifications are {.val {c('u','f','m')}} for 'unspecified', and cis 'male', and 'female' genders, respectively."))
  
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
                                      fname = "mhspitch", beginTime = .beginTime, 
                                      centerTime = centerTime, endTime = .endTime, 
                                      windowShift = windowShift, gender = gender, 
                                      maxF = maxF, minF = minF, 
                                      minAmp = minAmp, minAC1 = minAC1, 
                                      minRMS = minRMS, maxZCR = maxZCR, 
                                      minProb = minProb, plainSpectrum = plainSpectrum, 
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

attr(mhspitch,"ext") <-  "pit" 
attr(mhspitch,"tracks") <-  c("pitch[Hz]")
attr(mhspitch,"outputType") <-  "SSFF"
attr(mhspitch,"nativeFiletypes") <-  c("wav","au","kay","nist","nsp")
attr(mhspitch,"suggestCaching") <-  FALSE

attr(pitch_mhs,"ext") <-  "pit" 
attr(pitch_mhs,"tracks") <-  c("pitch[Hz]")
attr(pitch_mhs,"outputType") <-  "SSFF"
attr(pitch_mhs,"nativeFiletypes") <-  c("wav","au","kay","nist","nsp")
attr(pitch_mhs,"suggestCaching") <-  FALSE

### INTERACTIVE TESTING
#
#f <- normalizePath(list.files(file.path("..","inst","samples"),recursive = TRUE,full.names = TRUE))
#f <- f[grepl("*.aiff",f)]
#f <- f[!grepl("*.pit$",f,)]
#f <- f[1:2]

#mhspitch(f,beginTime=1.2, endTime=2.2, toFile=FALSE,keepConverted = FALSE,verbose = TRUE) -> a1
#mhspitch(f, toFile=TRUE,keepConverted = FALSE,verbose = TRUE) -> a2

#r <- normalizePath(list.files(file.path("..","inst","samples"),recursive = TRUE,full.names = TRUE,pattern = attr(fo_mhs,"ext")))
#unlink(r)
