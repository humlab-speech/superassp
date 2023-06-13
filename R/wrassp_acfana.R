##' Analysis of short-term autocorrelation function of signals
##' 
##' @description
##' Applies the autocorrelation function to windows of the input signals listen
##' in `listOfFiles`. Input signals not in the native "wav" file format will be converted before the autocorrelation functions are computed.
##' The conversion process will display warnings about input files that are not in known losslessly encoded formats.
##'
##' The results will be will be written to an SSFF formated file with the base
##' name of the input file and extension *.acf* in a track *acf*.
##' 
##' @details
##' The function is a re-write of the [wrassp::acfana] function, but with media
##' pre-conversion, better checking of preconditions such as the input file
##' existance, structured logging, and the use of a more modern framework for
##' user feedback. 
##' 
##' The native file type of this function is "wav" files (in "pcm_s16le" format). Input signal conversion, when needed, is done by [libavcodec](https://ffmpeg.org/libavcodec.html) and the excellent [av] wrapper package.
##' 
##' 
##' 

##' @param listOfFiles vector of file paths to be processed by function

##' @param beginTime the time point (in seconds) of the start of the analysed
##'   interval. A NULL or 0 is interpreted as the start of the signal file.
##' @param centerTime sets a single-frame analysis time point (in seconds).
##'   Overrides `beginTime`, `endTime` and `windowShift` parameters.
##' @param endTime the time point (in seconds) of the end of the analysed
##'   interval. A NULL or 0 is interpreted as the end of the signal file.
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
##' @param knownLossless a list of file extensions associated with known
##'   lossless file encodings.
##' @param logToFile whether to log commands to a separate logfile in the
##'   `outputDirectory`. Logging will otherwise be in the `acfana` logging
##'   namespace of [logger] and will be put wherever this namespace is defined to place its output.
##'   See [logger::log_appender] for details.
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
##' path2wav <- list.files(system.file("extdata", package = "wrassp"),
##'                        pattern = glob2rx("*.wav"),
##'                        full.names = TRUE)[1]
##'
##' # calculate short-term autocorrelation
##' res <- acfana(path2wav, toFile=FALSE)
##'
##' # plot short-term autocorrelation values
##' matplot(seq(0,numRecs.AsspDataObj(res) - 1) / rate.AsspDataObj(res) +
##'         attr(res, 'startTime'),
##'         res$acf,
##'         type='l',
##'         xlab='time (s)',
##'         ylab='short-term autocorrelation values')
##'
##' @export
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
                   knownLossless = c("wav","flac","aiff","wv","tta","caf"),
                   logToFile = FALSE,
                   convertOverwrites=FALSE,
                   keepConverted=FALSE,
                   verbose = TRUE) {
  
  ## Initial constants
  funName <- "acfana"
  nativeFiletypes <- c("wav")
  preferedFiletype <- nativeFiletypes[[1]]
  currCall <- rlang::current_call()
  
  if(is.null(analysisOrder)) analysisOrder <- 0 # How the C function expects the argument
  if(is.null(beginTime)) beginTime <- 0 # How the C function expects the argument
  if(is.null(endTime)) endTime <- 0 # How the C function expects the argument
  
  #### Setup logging of the function call ####
  
  
  if (!is.null(outputDirectory) && is.character(outputDirectory)) {
    outputDirectory = normalizePath(path.expand(outputDirectory),mustWork = FALSE)
    finfo  <- file.info(outputDirectory)
    if (is.na(finfo$isdir))
      if (!dir.create(outputDirectory, recursive=TRUE))
        cli::cli_abort("Unable to create the output directory {.path outputDirectory}.")
    else if (!finfo$isdir)
      cli::cli_abort("The path {.path outputDirectory} exists but is not a directory.")
    
    if(logToFile){
        logger::log_appender()
        cli::cli_inform("Storing the processing log in {.path {outputDirectory}}.")
        logger::log_threshold(logger::TRACE,namespace=funName)
        logger::log_layout(logger::layout_glue_generator(
          format = "{level} [{format(time, \"%Y-%m-%d %H:%M:%S\")}] {msg}"
        ),namespace=funName)
        logger::log_appender(logger::appender_file(
          file=file.path(outputDirectory,paste(funName,"log",sep="."))
        ),namespace=funName)
        logger::log_trace("{currCall}",namespace=funName)
    }
    
  }
  
  #### Progress bars #####
  
  if(verbose){
    process_pb <- list(name="Applying DSP function",
                       format="{cli::pb_extra$currFunName} {cli::pb_bar} {cli::pb_current}/{cli::pb_total}",
                       show_after=1,
                       clear=FALSE,
                       extra=list(currFunName=funName)
                       )
    convert_pb <- list(name="Converting media files",
                       format="Converting non-{.field {preferedFiletype}} media files to {.field {preferedFiletype}} format {cli::pb_bar} {cli::pb_current}/{cli::pb_total}",
                       show_after=1,
                       clear=FALSE)
  }else{
    process_pb <- FALSE
    convert_pb <- FALSE
  }

  #### [*] Input file conversion ####
  
  ## Check and fix input file paths
  listOfFiles <- prepareFiles(listOfFiles)
  
  if (is.null(listOfFiles) || length(listOfFiles) == 0 || ! all(file.exists(listOfFiles)) ) {
    cli::cli_abort(c("!"="The {.arg listOfFiles} has to contain a vector of working full paths to speech recordings."))
  }
  
  if(verbose) cli::cli_inform("Applying the {.fun {funName}} DSP function to {cli::no(length(listOfFiles))} speech recordings")

  # Not used now
  getaudiotype <- function(x){
    out <- rep(NA,length(x))
    for(c in seq_along(x)){
      out[c] <- av::av_media_info(x[c])$audio$codec
    }
    return(out)
  }
  
  notLossless <- listOfFiles[! tools::file_ext(listOfFiles) %in% knownLossless]

  if(length(notLossless) > 0){
    cli::cli_warn(c("w"="Found {.val {length(notLossless)}} recording{?s} that may have lossy compression",
                    "i"="If lossy compression was used when storing the signal, the result {.fun acfana} may not be accurate",
                    "x"="Please use known lossless formats (file extensions {.or {.val { knownLossless}}}) for acoustic analysis of speech recordings."))
  }
  
  # Convertion code for non-"native" files
  isNotNative <- ! (tools::file_ext(listOfFiles) %in% nativeFiletypes)
  
  listOfFilesDF <- data.frame(audio=listOfFiles,
                              isNative = ! isNotNative,
                              output = paste(tools::file_path_sans_ext(listOfFiles),preferedFiletype,sep=".")
                              )


  toConvert <- subset(listOfFilesDF, ! isNative)
  
  if(!convertOverwrites){ 
    toConvert <- subset( toConvert,!file.exists(output))
  }
  toConvert$isNative <- NULL
  

  purrr::pwalk(.l=toConvert,.f=av::av_audio_convert,verbose = FALSE,channels=1,.progress = convert_pb)

  
  # Here we explicitly make a new listOfFiles that points to wav files
  # so that legacy code can be used
  listOfFiles <- listOfFilesDF$output
  
  #### Application of DSP C function  ####
  
  
  if(!isAsspWindowType(window)){
    cli::cli_abort("WindowFunction of type {.val window} is not supported!")
  }
  
  


  

  
  insideFunction <- function(x){
    
      
      
    externalRes = invisible(.External("performAssp", x, 
                                      fname = "acfana", beginTime = beginTime, 
                                      centerTime = centerTime, endTime = endTime, 
                                      windowShift = windowShift, windowSize = windowSize, 
                                      effectiveLength = effectiveLength, window = window, 
                                      analysisOrder = as.integer(analysisOrder), energyNormalization = energyNormalization, 
                                      lengthNormalization = lengthNormalization, toFile = toFile, 
                                      explicitExt = explicitExt, progressBar = NULL,
                                      outputDirectory = outputDirectory, PACKAGE = "superassp"))
    
    return(externalRes)
  }



  if(toFile){
    externalRes <- purrr::walk(.x=listOfFiles,.f=insideFunction,.progress = process_pb)
  }else{
    externalRes <- purrr::map(.x=listOfFiles,.f=insideFunction,.progress = process_pb)
  }
  
  #### [*] Cleanup of possibly converted files  ####
  
  
  #Clear the wavs created in the conversion step
  # It is assumed that files that had not been created in the conversion should also not be cleared
  if(! keepConverted && nrow(toConvert) > 0){
    if(verbose) cli::cli_inform("Cleaning up temporary (converted) versions of media files")
    purrr::pwalk(.l=toConvert,.f= \(audio, output) unlink(output,recursive = FALSE, force = FALSE, expand = FALSE))
  }

  #Simplify output if just one AsspDataObj is returned
  
  if(length(listOfFiles) == 1) externalRes <- purrr::pluck(externalRes,1)
  
  return(externalRes)
}
attr(acfana,"ext") <-  "acf" 
attr(acfana,"tracks") <-  c("acf")
attr(acfana,"outputType") <-  "SSFF"
attr(acfana,"nativeFiletypes") <-  c("wav")



### INTERACTIVE TESTING
#f <- list.files("~/Desktop/input/",full.names = TRUE)
#acfana(f,toFile=FALSE,keepConverted = FALSE,outputDirectory = "/Users/frkkan96/Desktop/output/",verbose = FALSE) -> a


