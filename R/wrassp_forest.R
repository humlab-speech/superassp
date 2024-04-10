##' Estimate formant frequencies and bandwidths
##'
##' Formant estimation of the signal(s) in `listOfFiles`. Raw
##' resonance frequency and bandwidth values are obtained by root-solving of the
##' Linear Prediction polynomial from the autocorrelation method and the
##' Split-Levinson-Algorithm (SLA). Resonances are then classified as formants
##' using the so-called Pisarenko frequencies (by-product of the SLA) and a
##' formant frequency range table derived from the nominal \ifelse{html}{\out{F<sub>1</sub>}}{\eqn{F_1}}. The
##' latter may have to be increased by about 12% for female voices (see
##' `nominalF1` and `gender` parameters). This function uses the *libassp* C library 
##'  \insertCite{s5h}{superassp} for the DSP work.
##'
##' Input signals not in a natively supported file format will be converted
##' before the autocorrelation functions are computed. The conversion process
##' will display warnings about input files that are not in known losslessly
##' encoded formats.
##'
##' Default output is in SSFF binary format, with tracks containing the
##' estimated mid formant frequency of each formant (track *F[Hz]*, one column per
##' formant) and the associated formant bandwidth  (track *B[Hz]*, one column per
##' formant). If `toFile` is `TRUE`, the results will be written to a file with the
##' same name as the input file, but with an extension *.fms*.
##'
##' @details The function is a re-write of the [wrassp::forest] function, but
##' with media pre-conversion, better checking of preconditions such as the
##' input file existance, structured logging, and the use of a more modern
##' framework for user feedback.
##'
##' The native file type of this function is "wav" files (in "pcm_s16le"
##' format), SUNs "au", NIST, or CSL formats (kay or NSP extension). Input
##' signal conversion, when needed, is done by
##' [libavcodec](https://ffmpeg.org/libavcodec.html) and the excellent [av]
##' wrapper package.
##' 
##' @note
##' This function is not considered computationally expensive enough to require caching of 
##' results if applied to many signals. However, if the number of signals it will be applied to 
##' is *very* long, then caching of results may be warranted.
##' 
##' 
##' @inheritParams ksvfo
##' @param effectiveLength make window size effective rather than exact
##' @param nominalF1 = The nominal (assumed)  \ifelse{html}{\out{F<sub>1</sub>}}{\eqn{F_1}} frequency (default: 500.0 Hz)
##' @param gender = Use gender specific parameters? Permitted codes are  "f"[emale], "m"[ale] or "u"[nknown]. When "f", the effective window length is set to 12.5 ms and the nominal \ifelse{html}{\out{F<sub>1</sub>}}{\eqn{F_1}} to 560 Hz.
##' @param estimate insert rough frequency estimates of missing formants? By default, the frequency is set to zero.
##' @param order decrease default LPC filter order by 2 (one resonance less)
##' @param incrOrder increase default LPC filter order by 2 (one resonance more)
##' @param numFormants = The number of formants to identify. Defaults to 4, and the maximum value is 8 or half the LPC filter order)
##' @param window = <type>: set analysis window function to <type> (default: BLACKMAN)
##' @param preemphasis = <val>: set pre-emphasis factor to <val> (-1 <= val <= 0) 
##' (default: dependent on sample rate and nominal \ifelse{html}{\out{F<sub>1</sub>}}{\eqn{F_1}})
##' @author Raphael Winkelmann
##' @author Lasse Bombien
##' @author Fredrik Nylén 
##' 
##' @return If `toFile` is `FALSE`, the function returns a list of [AsspDataObj]
##'   objects. If `toFile` is `TRUE`, the number (integer) of successfully
##'   processed and stored output files is returned.
##'
##' @seealso [wrassp::acfana]
##' @seealso [superassp::AsspWindowTypes]
##' @seealso [av::av_audio_convert]
##' 
##' @useDynLib superassp, .registration = TRUE
##' @examples
##' # get path to audio file
##' path2wav <- list.files(system.file("samples","sustained", package = "superassp"), pattern = glob2rx("a1.wav"), full.names = TRUE)
##'
##' 
##' # calculate formant values
##' res <- forest(path2wav, toFile=FALSE)
##' 
##' # plot formant values
##' matplot(seq(0,numRecs.AsspDataObj(res) - 1) / rate.AsspDataObj(res) + 
##'           attr(res, 'startTime'), 
##'         res[["F[Hz]"]], 
##'         type='l', 
##'         xlab='time (s)', 
##'         ylab='Formant frequency (Hz)')
##' 
##' @export
forest <- function(listOfFiles = NULL,
                   beginTime = 0.0,
                   endTime = 0.0,
                   windowShift = 5.0,
                   windowSize = 20.0,
                   effectiveLength = TRUE,
                   nominalF1 = 500,
                   gender = 'm',
                   estimate = FALSE,
                   order = 0,
                   incrOrder = 0,
                   numFormants = 4,
                   window = 'BLACKMAN',
                   preemphasis = -0.8,
                   toFile = TRUE,
                   explicitExt = "fms",
                   outputDirectory = NULL,
                   assertLossless = NULL,
                   logToFile = FALSE,
                   convertOverwrites=FALSE,
                   keepConverted=FALSE,
                   verbose = TRUE) {
  ## Initial constants -- specific to this function
  explicitExt <- ifelse(is.null(explicitExt),"fms",explicitExt)
  newTracknames <- c("F[Hz]","B[Hz]") ## Only used if SSFF tracks needs to be renamed from the called function (in C) before returning the SSFF track obj 
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
                                      fname = "forest", beginTime =  .beginTime, 
                                      endTime = .endTime, windowShift = windowShift, 
                                      windowSize = windowSize, effectiveLength = effectiveLength, 
                                      nominalF1 = nominalF1, gender = gender, 
                                      estimate = estimate, order = as.integer(order), 
                                      incrOrder = as.integer(incrOrder), numFormants = as.integer(numFormants), 
                                      window = window, preemphasis = preemphasis, 
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

attr(forest,"ext") <-  "fms" 
attr(forest,"tracks") <-  c("F[Hz]","B[Hz]")
attr(forest,"outputType") <-  "SSFF"
attr(forest,"nativeFiletypes") <-  c("wav","au","kay","nist","nsp")
attr(forest,"suggestCaching") <-  FALSE


### INTERACTIVE TESTING
#
#f <- normalizePath(list.files(file.path("..","inst","samples"),recursive = TRUE,full.names = TRUE))
#f <- f[grepl("*.aiff",f)]
#f <- f[!grepl("*.fms$",f,)]
#f <- f[1:2]

#forest(f,beginTime=1.2, endTime=2.2, toFile=FALSE,keepConverted = FALSE,verbose = TRUE) -> a1
#forest(f, toFile=FALSE,keepConverted = FALSE,verbose = TRUE) -> a2

#r <- normalizePath(list.files(file.path("..","inst","samples"),recursive = TRUE,full.names = TRUE,pattern = attr(forest,"ext")))

#unlink(r)

