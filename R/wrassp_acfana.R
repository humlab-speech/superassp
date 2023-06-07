##' Analysis of short-term autocorrelation function of signals
##'
##' The function applies the autocorrelation function to windows of the input signals listen in `listOfFiles`.
##' The results will be will be written to an SSFF formated file with the
##' base name of the input file and extension *.acf* in a track *acf*.
##' 
##' @title acfana
##' @param listOfFiles vector of file paths to be processed by function
##' @param optLogFilePath path to option log file
##' @param beginTime = <time>: set begin of analysis interval to <time> seconds (default: 0 = beginning of file)
##' @param centerTime = <time>: set single-frame analysis with the analysis window centred at <time> seconds; 
##' overrules BeginTime, EndTime and WindowShift options
##' @param endTime = <time>: set end of analysis interval to <time> seconds (default: 0 = end of file)
##' @param windowShift = <dur>: set analysis window shift to <dur> ms (default: 5.0)
##' @param windowSize = <dur>: set analysis window size to <dur> ms; overrules EffectiveLength parameter
##' @param effectiveLength make window size effective rather than exact
##' @param window = <type>: set analysis window function to <type> (default: BLACKMAN)
##' @param analysisOrder = <num>: set analysis order to <num> (default: 0 = sample rate in kHz + 3)
##' @param energyNormalization calculate energy-normalized autocorrelation
##' @param lengthNormalization calculate length-normalized autocorrelation
##' @param toFile write results to file (default extension is .acf)
##' @param explicitExt set if you wish to override the default extension
##' @param outputDirectory directory in which output files are stored. Defaults to NULL, i.e. 
##' the directory of the input files
##' @param verbose display infos & show progress bar
##' @return A list of objects of class [AsspDataObj]
##' @author Raphael Winkelmann
##' @author Lasse Bombien
##' @author Fredrik Karlsson
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
                   optLogFilePath = NULL,
                   beginTime = 0.0,
                   centerTime = FALSE,
                   endTime = 0.0,
                   windowShift = 5.0,
                   windowSize = 20.0,
                   effectiveLength = TRUE,
                   window = "BLACKMAN",
                   analysisOrder = 0,
                   energyNormalization = FALSE,
                   lengthNormalization = FALSE,
                   toFile = TRUE,
                   explicitExt = "acf",
                   outputDirectory = NULL,
                   knownLossless = c("pcm_s16le", "flac"),
                   forceToLog = FALSE,
                   convertOverwrites=FALSE,
                   keepConverted=FALSE,
                   verbose = TRUE) {
  
  ## Initial constants
  funName <- "acfana"
  
  if(verbose){
    process_pb <- list(name="Applying DSP function",
                       format="{cli::pb_name} {cli::pb_extra$currFunName} {cli::pb_bar} {cli::pb_current}/{cli::pb_total}",
                       show_after=1,
                       clear=FALSE,
                       extra=list(currFunName=funName)
                       )
    convert_pb <- list(name="Converting media files",
                       format="Converting non-{.field wav} media files to {.field wav} format {cli::pb_bar} {cli::pb_current}/{cli::pb_total}",
                       show_after=1,
                       clear=FALSE)
  }else{
    process_pb <- NULL
    convert_pb <- NULL
  }

  ## Check and fix input file paths
  listOfFiles <- prepareFiles(listOfFiles)
  
  if (is.null(listOfFiles) || length(listOfFiles) == 0 || ! all(file.exists(listOfFiles)) ) {
    cli::cli_abort(c("!"="The {.arg listOfFiles} has to contain a vector of working full paths to speech recordings."))
  }
  
  if(verbose) cli::cli_inform("Processing {cli::no(length(listOfFiles))} speech recordings using the {.fun {funName}} DSP function")

  getaudiotype <- function(x){
    out <- rep(NA,length(x))
    for(c in seq_along(x)){
      out[c] <- av::av_media_info(x[c])$audio$codec
    }
    return(out)
  }
  notLossless <- listOfFiles[! getaudiotype(listOfFiles) %in% knownLossless]

  if(length(notLossless) > 0){
    cli::cli_warn(c("w"="Found {.val {length(notLossless)}} recording{?s} with lossy compression",
                    "i"="If the signal has been stored with lossy compression the result {.fun acfana} may not be accurate",
                    "x"="Please use known lossless formats ({.or {.val { knownLossless}}}) for speech recordings"))
  }
  
  # Convertion code for non-wav files
  isNotWavs <- ! (tools::file_ext(listOfFiles) == "wav")
  
  listOfFilesDF <- data.frame(audio=listOfFiles,
                              isWav = (tools::file_ext(listOfFiles) == "wav"),
                              output = paste(tools::file_path_sans_ext(listOfFiles),"wav",sep=".")
                              )


  toConvert <- subset(listOfFilesDF, ! isWav)
  
  if(!convertOverwrites){ 
    toConvert <- subset( toConvert,!file.exists(output))
  }
  toConvert$isWav <- NULL
  

  purrr::pwalk(.l=toConvert,.f=av::av_audio_convert,verbose = FALSE,channels=1,.progress = convert_pb)

  
  # Here we explicitly make a new listOfFiles that points to wav files
  # so that legacy code can be used
  listOfFiles <- listOfFilesDF$output
  
  ## END OF CONVERSION CODE
  
  
  if(!isAsspWindowType(window)){
    cli::cli_abort("WindowFunction of type {.val window} is not supported!")
  }
  
  
  if (!is.null(outputDirectory) ) {
    if(toFile){
      outputDirectory = normalizePath(path.expand(outputDirectory),mustWork = FALSE)
      finfo  <- file.info(outputDirectory)
      if (is.na(finfo$isdir))
        if (!dir.create(outputDirectory, recursive=TRUE))
          cli::cli_abort("Unable to create the output directory {.path outputDirectory}.")
      else if (!finfo$isdir)
        cli::cli_abort("The path {.path outputDirectory} exists but is not a directory.")
    }else{
      if(forceToLog && verbose){
        cli::cli_inform("Storing the processing log in {.path {outputDirectory}}.")
      }else{
        cli::cli_abort("You have specified an output directory but neither result or log files should be stored on disk.")
      }
    }

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
  
  #Clear the wavs created in the conversion step
  # It is assumed that files that had not been created in the conversion should also not be cleared
  if(! keepConverted){
    purrr::pwalk(.l=toConvert,.f= \(audio, output) unlink(output,recursive = FALSE, force = FALSE, expand = FALSE))
  }

  
  return(externalRes)
}
attr(acfana,"ext") <-  "acf" 
attr(acfana,"tracks") <-  c("acf")
attr(acfana,"outputType") <-  "SSFF"


### INTERACTIVE TESTING
#f <- list.files("~/Desktop/input/",full.names = TRUE)
#acfana(f,toFile=FALSE,keepConverted = FALSE,outputDirectory = "/Users/frkkan96/Desktop/output/") -> a


