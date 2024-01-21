##' Short-term Root Mean Square amplitude of signals
##'
##' @description The RMS amplitude is computed for each window of `windowSize`
##' length in the input signals files liste in `listOfFiles`. Per default, the
##' RMS values are expressed in decibel (dB) so that they correspond to the
##' short-term power of the signal. Input signals not in the native "wav" file
##' format will be converted before the function is applied.
##' The conversion process will display warnings about input files that are not
##' in known losslessly encoded formats.
##'
##' The results will be will be written to an SSFF formated file with the base
##' name of the input file and extension *.rms* in a track *rms*.
##'
##' @details The function is a re-write of the [wrassp::rmsana] function, but
##' with media pre-conversion, better checking of preconditions such as the
##' input file existence, structured logging, and the use of a more modern
##' framework for user feedback.
##'
##' The native file type of this function is "wav" files (in "pcm_s16le"
##' format). Input signal conversion, when needed, is done by
##' [libavcodec](https://ffmpeg.org/libavcodec.html) and the [av]
##' wrapper package.
##'
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
##' path2wav <- list.files(system.file("samples","sustained", package = "superassp"), pattern = glob2rx("a1.wav"), full.names = TRUE)
##'
##' # calculate short-term autocorrelation
##' res <- rmsana(path2wav, toFile=FALSE)
##'
##'
##' # plot rms values
##' plot(seq(0,numRecs.AsspDataObj(res) - 1) / rate.AsspDataObj(res) +
##'        attr(res, 'startTime'),
##'      res$rms,
##'      type='l',
##'      xlab='time (s)',
##'      ylab='RMS energy (dB)')
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
                   knownLossless = c("wav","flac","aiff","wv","tta","caf"),
                   logToFile = FALSE,
                   convertOverwrites=FALSE,
                   keepConverted=FALSE,
                   verbose = TRUE) {
  
  ## Initial constants
  funName <- "rmsana"
  nativeFiletypes <- c("wav")
  preferedFiletype <- nativeFiletypes[[1]]
  currCall <- rlang::current_call()
  
  if(is.null(beginTime)) beginTime <- 0 # How the C function expects the argument
  if(is.null(endTime)) endTime <- 0 # How the C function expects the argument
  
  #### Setup logging of the function call ####
  
  
  if (!is.null(outputDirectory) && is.character(outputDirectory)) {
    outputDirectory = normalizePath(path.expand(outputDirectory),mustWork = FALSE)
    
    if( file.exists(outputDirectory) ){
      finfo  <- file.info(outputDirectory)
      if (! finfo$isdir)
        cli::cli_abort("The path {.path outputDirectory} exists but is not a directory.")
    }else{
      if (!dir.create(outputDirectory, recursive=TRUE))
        cli::cli_abort("Unable to create the output directory {.path outputDirectory}.")
    }
    
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
                       format="Preparation: converting to {.field {cli::pb_extra$preferedFT}} format {cli::pb_bar} {cli::pb_current}/{cli::pb_total}",
                       show_after=1,
                       clear=FALSE,
                       extra=list(preferedFT=preferedFiletype)
    )
  }else{
    process_pb <- FALSE
    convert_pb <- FALSE
  }
  
  #### [*] Input file conversion ####
  
  ## Check and fix input file paths
  listOfFiles = gsub("^file://","", listOfFiles)
  listOfFiles = normalizePath(path.expand(listOfFiles))
  
  if (is.null(listOfFiles) || length(listOfFiles) == 0 || ! all(file.exists(listOfFiles)) ) {
    cli::cli_abort(c("!"="The {.arg listOfFiles} contains no working full paths to speech recordings."))
  }
  
  if(verbose) cli::cli_inform("Applying the {.fun {funName}} DSP function to {cli::no(length(listOfFiles))} speech recording{?s}")
  
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
    cli::cli_warn(c("w"="Found {.val {length(notLossless)}} recording{?s} stored in lossy compression formats",
                    "i"="If lossy compression was used when storing the signal, the result {.fun {funName}} may not be accurate",
                    "x"="Please use a lossless file format (file extensions {.or {.val { knownLossless}}}) for acoustic analysis of speech recordings if possible."))
    if(verbose)
      cli::cli_inform(c("i"="Lossy files: {.file {basename(notLossless)}}"))
  }
  
  # Convertion code for non-"native" files
  
  listOfFilesDF <- data.frame(audio=listOfFiles,
                              output = paste(tools::file_path_sans_ext(listOfFiles),preferedFiletype,sep=".")
  )
  
  
  toConvert <- subset(listOfFilesDF, ! (tools::file_ext(listOfFiles) %in% nativeFiletypes) )
  
  if(!convertOverwrites){ 
    toConvert <- subset( toConvert,!file.exists(output) && audio != output)
  }
  
  if(nrow(toConvert) > 0 ){
    cli::cli_inform(c("Found {.val {nrow(toConvert)}} recording{?s} that require conversion",
                      "x"="If a file format is not nativelly suppored by {.fun {funName}} it will have to be converted before application of the function",
                      "i"="Please use {.or {.val {nativeFiletypes}}} which are nativelly supported by {.fun {funName}} to eliminate this conversion."))
    if(verbose){
      
      cli::cli_inform(c("i"="Converted files: {.file {basename(toConvert$audio)}}"))
    }
    
    
    purrr::pwalk(.l=toConvert,.f=av::av_audio_convert,verbose = FALSE,channels=1,format=NULL,.progress = convert_pb)
    
    # Here we explicitly make a new listOfFiles that points to wav files
    # so that upcoming legacy code can be used
    listOfFiles <- listOfFilesDF$output
    
  }
  
  
  assertthat::assert_that(all(tools::file_ext(listOfFiles) %in% nativeFiletypes )) #Make sure that we have a file that may now be handled
  
  #### Application of DSP C function  ####
  
  
  if(!isAsspWindowType(window)){
    cli::cli_abort("WindowFunction of type {.val window} is not supported!")
  }
  
  insideFunction <- function(x){
    
    externalRes = invisible(.External("performAssp", listOfFiles, 
                                      fname = "rmsana", beginTime = beginTime, 
                                      centerTime = centerTime, endTime = endTime, 
                                      windowShift = windowShift, windowSize = windowSize, 
                                      effectiveLength = effectiveLength, linear = linear, 
                                      window = window, toFile = toFile, 
                                      explicitExt = explicitExt, 
                                      progressBar = NULL, outputDirectory = outputDirectory,
                                      PACKAGE = "superassp"))
    
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

attr(rmsana,"ext") <-  "rms" 
attr(rmsana,"tracks") <-  c("acf")
attr(rmsana,"outputType") <-  "SSFF"
attr(rmsana,"nativeFiletypes") <-  c("wav","au","kay","nist")


### INTERACTIVE TESTING
#
#f <- normalizePath(list.files(file.path("..","inst","samples","sustained"),full.names = TRUE))
#f <- f[grepl("*.aiff",f)]

#rmsana(f,toFile=FALSE,keepConverted = FALSE,outputDirectory = "/Users/frkkan96/Desktop/output/",verbose = TRUE,convertOverwrites=TRUE) -> a 
