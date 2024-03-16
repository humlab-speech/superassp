##' Analysis of short-term autocorrelation function of signals
##'
##' @description Applies the autocorrelation function to windows of the input
##' signals listed in `listOfFiles`. Input signals not in a file format natively
##' supported will be converted before the autocorrelation functions are
##' computed. The conversion process will display warnings about input files
##' that are not in known losslessly encoded formats.
##'
##' The results will be will be written to an SSFF formated file with the base
##' name of the input file and extension *.acf* in a track *acf*.
##'
##' @details The function is a re-write of the [wrassp::acfana] function, but
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
##' is *very* long, then caching of results may be warranted.
##' 
##' @param listOfFiles vector of file paths to be processed by function
##' @param beginTime the time point (in seconds) of the start of the analysed
##'   interval. A NULL or 0 is interpreted as the start of the signal file. 
##'   If a vector of time points is supplied, the length of that vector needs 
##'   to correspond with the length of `listOfFiles`.
##' @param centerTime sets a single-frame analysis time point (in seconds).
##'   Overrides `beginTime`, `endTime` and `windowShift` parameters.
##' @param endTime the time point (in seconds) of the end of the analysed
##'   interval. A NULL or 0 is interpreted as the end of the signal file. 
##'   If a vector of time points is supplied, the length of that vector needs 
##'   to correspond with the length of `listOfFiles`.
##' @param windowShift the amount of time (in ms) that the analysis window will
##'   be shifted between analysis frames
##' @param windowSize the analysis window size (in ms); overrides the effect of
##'   the `effectiveLength` parameter
##' @param effectiveLength make window size effective rather than exact
##' @param window = the analysis window function type ("BLACKMAN" by default).
##'   See [superassp::AsspWindowTypes] for a list of supported window types.
##' @param analysisOrder the analysis order. The `NULL` or `0` sets the analysis
##'   order to the sample rate (in kHz) + 3, so that a signal with a 16000 Hz
##'   sampling rate will be analysed using an `analysisOrder` of 19.
##' @param energyNormalization calculate energy-normalized autocorrelation
##' @param lengthNormalization calculate length-normalized autocorrelation
##' @param toFile Should the function write the results to a file, with the
##'   (default) file extension *.acf* (`TRUE`) or returned as a list of
##'   [AsspDataObj] objects (`FALSE`)?
##' @param explicitExt by default an *.acf* file extension will be used when
##'   result files are written (`toFile=TRUE`), but the file extension can be
##'   set to something else using this function argument.
##' @param outputDirectory directory in which output files are stored. Defaults
##'   to NULL which means that the result file will be stored in the same
##'   directory as the input file.
##' @param verbose display verbose information about processing steps taken, as
##'   well as progress bars.
##' @param assertLossless an optional list of file extensions that the user wants to assert 
##'   contains losslessly encoded signals data.
##' @param logToFile whether to log commands to a separate logfile in the
##'   `outputDirectory`. Logging will otherwise be in the `acfana` logging
##'   namespace of [logger] and will be put wherever this namespace is defined to place its output.
##'   See [logger::log_appender] for details.
##'
##' @return The number of successfully written files (if `toFile=TRUE`), or a vector of `AsspDataObj` objects (if `toFile=FALSE`).
##'
##' @seealso [wrassp::acfana]
##' @seealso [superassp::AsspWindowTypes]
##' @seealso [av::av_audio_convert]
##'
##' @useDynLib superassp, .registration = TRUE
##' @examples
##' # get path to audio file
#' path2wav <- list.files(system.file("samples","sustained", package = "superassp"), pattern = glob2rx("a1.wav"), full.names = TRUE)
#'
#' # calculate short-term autocorrelation
#' res <- acfana(path2wav, toFile=FALSE)
#'
#' # plot short-term autocorrelation values
#' matplot(seq(0,numRecs.AsspDataObj(res) - 1) / rate.AsspDataObj(res) +
#'         attr(res, 'startTime'),
#'         res$acf,
#'         type='l',
#'         xlab='time (s)',
#'         ylab='short-term autocorrelation values')
##'
##' @export
##' @references 
##'   \insertAllCited{}
##'
   
acfana <- function(listOfFiles = NULL,
                   beginTime = 0.0,
                   centerTime = FALSE,
                   endTime = 0.0,
                   windowShift = 5.0,
                   windowSize = 20.0,
                   effectiveLength = TRUE,
                   window = "BLACKMAN",
                   analysisOrder = NULL,
                   energyNormalization = FALSE,
                   lengthNormalization = FALSE,
                   toFile = TRUE,
                   explicitExt = "acf",
                   outputDirectory = NULL,
                   assertLossless = NULL,
                   logToFile = FALSE,
                   keepConverted=FALSE,
                   verbose = TRUE) {
  

  
  ## Initial constants -- specific to this function
  explicitExt <- ifelse(is.null(explicitExt),"fo",explicitExt)
  newTracknames <- NULL ## Only used if SSFF tracks needs to be renamed from the called function (in C) before returning the SSFF track obj 
  nativeFiletypes <- c("wav","au","kay","nist","nsp")
  
  if(is.null(analysisOrder)) analysisOrder <- 0 # How the C function expects the argument
  if(!isAsspWindowType(window)){
    cli::cli_abort("WindowFunction of type {.val window} is not supported!")
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
    
    ret <- invisible(.External("performAssp", x, 
                                      fname = "acfana", beginTime = .beginTime, 
                                      centerTime = centerTime, endTime = .endTime, 
                                      windowShift = windowShift, windowSize = windowSize, 
                                      effectiveLength = effectiveLength, window = window, 
                                      analysisOrder = as.integer(analysisOrder), energyNormalization = energyNormalization, 
                                      lengthNormalization = lengthNormalization, toFile = toFile, 
                                      explicitExt = explicitExt, progressBar = NULL,
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
attr(acfana,"ext") <-  "acf" 
attr(acfana,"tracks") <-  c("acf")
attr(acfana,"outputType") <-  "SSFF"
attr(acfana,"nativeFiletypes") <-  c("wav","au","kay","nist","nsp")
attr(acfana,"suggestCaching") <-  FALSE



### INTERACTIVE TESTING
#
#f <- normalizePath(list.files(file.path("..","inst","samples"),recursive = TRUE,full.names = TRUE))
 #f <- f[grepl("*.aiff",f)]
 
#acfana(f,beginTime=1.2, endTime=2.2, toFile=FALSE,keepConverted = FALSE,verbose = TRUE) -> a
#acfana(f, toFile=FALSE,keepConverted = FALSE,verbose = TRUE) -> a

#r <- normalizePath(list.files(file.path("..","inst","samples"),recursive = TRUE,full.names = TRUE,pattern = attr(acfana,"ext")))
#unlink(r)
#outputDirectory = "/Users/frkkan96/Desktop/output/"
#Error in acfana(f, toFile = FALSE, keepConverted = FALSE, outputDirectory = "/Users/frkkan96/Desktop/output/",  : 
#object '"/Users/frkkan96/Desktop/output/"' not found

