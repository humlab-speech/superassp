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
                   explicitExt = wrasspOutputInfos[["acfana"]]$ext,
                   outputDirectory = NULL,
                   knownLossless = c("wav", "flac"),
                   forceToLog = NULL,
                   verbose = TRUE) {
  
  
  
  ###########################
  # a few parameter checks and expand paths
  
  ###########################
  # Pre-process file list
  listOfFiles <- prepareFiles(listOfFiles)
  
  if (is.null(listOfFiles) ||! all(file.exists(listOfFiles))) {
    cli::cli_abort(c("!"="The {.arg listOfFiles} has to contain a vector of working full paths to speech recordings."))
  }
  
  notLossless <- listOfFiles[! tools::file_ext(listOfFiles) %in% knownLossless]

  if(length(notLossless) > 0){
    cli::cli_warn(c("w"="Found {no(notLossless)}) recordings stored likelly stored with lossy compression",
                    "i"="If the signal has been stored with lossy compression the result {.fun acfana} may not be accurate",
                    "x"="Please use known lossless formats ({.or {.val { knownLossless}}}) for speech recordings"))
  }
  
  # Convertion code for non-wav files
  isNotWavs <- ! (tools::file_ext(listOfFiles) == "wav")
  
  listOfFilesDF <- data.frame(originalListOfFiles=listOfFiles,
                              isWav = ! (tools::file_ext(listOfFiles) == "wav"),
                              listOfFiles = paste(tools::file_path_sans_ext(listOfFiles),"wav",sep=".")
                              )


  notWavs <- nrow(listOfFilesDF[! listOfFilesDF$isWav,])
  
  if(nrow(notWavs) > 0){
    
    for(t in 1:length(notWavs)){
      fromFile <- notWavs[[r,"originalListOfFiles"]]
      outputFile <- notWavs[[r,"listOfFiles"]]
      
      if(! file.exists(outputFile)){ #Perhaps left from a previous attempt
        cli::cli_inform(c("i"="Converting {basename(fromFile)}."))
        av::av_audio_convert(audio = fromFile,output = outputFile,verbose = FALSE,channels=1)
        if(! file.exists(outputFile)){
          cli::cli_abort("Could not convert file {.path {basename(fromFile)}.")
        }
      }else{
        cli::cli_inform("Conversion of {.path {basename(fromFile)}} due to an existing {.field wav} version.")
      }
      
      
    }
  }
  # Here we explicitly make a new listOfFiles that points to wav files
  # so that legacy code can be used
  listOfFiles <- listOfFilesDF[[ ,listOfFiles]]
  
  ## END OF CONVERSION CODE
  
  
  if(!isAsspWindowType(window)){
    cli::cli_abort("WindowFunction of type {.val window} is not supported!")
  }
  
  if (!is.null(outputDirectory)) {
    outputDirectory = normalizePath(path.expand(outputDirectory))
    finfo  <- file.info(outputDirectory)
    if (is.na(finfo$isdir))
      if (!dir.create(outputDirectory, recursive=TRUE))
        cli::cli_abort("Unable to create the output directory {.path outputDirectory}.")
    else if (!finfo$isdir)
      cli::cli_abort("The path {.path outputDirectory} exists but is not a directory.")
  }

  

  
  insideFunction <- function(x){
    invisible(.External("performAssp", x, 
                                      fname = "acfana", beginTime = beginTime, 
                                      centerTime = centerTime, endTime = endTime, 
                                      windowShift = windowShift, windowSize = windowSize, 
                                      effectiveLength = effectiveLength, window = window, 
                                      analysisOrder = as.integer(analysisOrder), energyNormalization = energyNormalization, 
                                      lengthNormalization = lengthNormalization, toFile = toFile, 
                                      explicitExt = explicitExt, progressBar = NULL, # To be removed later
                                      outputDirectory = outputDirectory, PACKAGE = "superassp"))
  }

  pb <- 
  if(toFile){
    externalRes <- purrr::walk(.x=listOfFiles,.f=insideFunction,.progress = verbose)
  }else{
    externalRes <- purrr::map(.x=listOfFiles,.f=insideFunction,.progress = verbose)
  }
  

  ############################
  # write options to options log file
  
  # 
  # if (forceToLog){
  #   optionsGivenAsArgs = as.list(match.call(expand.dots = TRUE))
  #   wrassp.logger(optionsGivenAsArgs[[1]], optionsGivenAsArgs[-1],
  #                 optLogFilePath, listOfFiles)
  #   
  # }
  
  return(externalRes)
}
attr(acfana,"ext") <-  wrasspOutputInfos[["acfana"]]$ext 
attr(acfana,"tracks") <-  wrasspOutputInfos[["acfana"]]$tracks
attr(acfana,"outputType") <-  wrasspOutputInfos[["acfana"]]$outputType


### INTERACTIVE TESTING
#acfana(c("~/Desktop/test.wav","~/Desktop/ChatGPT/Skärminspelning 2023-01-19 kl. 18.41.25.mov")) -> a
