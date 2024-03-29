##' Calculate Linear Prediction smoothed spectrum using libassp
##'
##' Short-term spectral analysis of the signal in <listOfFiles>
##' using the Fast Fourier Transform and linear predictive smoothing.
##' Analysis results will be written to a file with the
##' base name of the input file and the spectrum type in
##' lower case as extension (i.e. '.lps').
##' Default output is in SSFF format with the
##' spectrum type in lower case as track name.
##' @title lpsSpectrum
##' @param listOfFiles vector of file paths to be processed by function 
##' @param optLogFilePath path to option log file
##' @param beginTime = <time>: set begin of analysis interval to <time> seconds
##' (default: begin of data)
##' @param centerTime = <time>: set single-frame analysis with the analysis
##' window centred at <time> seconds; overrules beginTime, endTime and
##' windowShift options
##' @param endTime = <time>: set end of analysis interval to <time> seconds
##' (default: end of data)
##' @param resolution = <freq>: set FFT length to the smallest value which
##' results in a frequency resolution of <freq> Hz or better (default: 40.0)
##' @param fftLength = <num>: set FFT length to <num> points (overrules default
##' and 'resolution' option)
##' @param windowShift = <dur>: set analysis window shift to <dur> ms
##' (default: 5.0)
##' @param window = <type>: set analysis window function to <type> (default:
##' BLACKMAN)
##' @param windowSize = <dur>: set effective analysis window size to <dur> ms 
##' @param order = <num>: set prediction order to <num> (default: sampling
##' rate in kHz + 3)
##' @param preemphasis = <val>: set pre-emphasis factor to <val> (default:
##' -0.95)
##' @param deemphasize (default: undo spectral tilt due to
##' pre-emphasis used in LP analysis, i.e. TRUE)
##' @param toFile write results to file (default extension depends on )
##' @param explicitExt set if you wish to override the default extension
##' @param outputDirectory directory in which output files are stored. Defaults to NULL, i.e.
##' the directory of the input files
##' @param forceToLog is set by the global package variable useWrasspLogger. This is set
##' to FALSE by default and should be set to TRUE is logging is desired.
##' @param verbose display infos & show progress bar
##' @return nrOfProcessedFiles or if only one file to process return
##' AsspDataObj of that file
##' @author Raphael Winkelmann
##' @author Lasse Bombien
##' @seealso \code{\link{dftSpectrum}}, \code{\link{cssSpectrum}}, \code{\link{cepstrum}}; 
##' all derived from libassp's spectrum function.
##' @useDynLib superassp, .registration = TRUE
##' @examples
##' # get path to audio file
##' path2wav <- list.files(system.file("extdata", package = "wrassp"), 
##'                        pattern = glob2rx("*.wav"), 
##'                        full.names = TRUE)[1]
##' 
##' # calculate lps spectrum
##' res <- lpsSpectrum(path2wav, toFile=FALSE)
##' 
##' # plot spectral values at midpoint of signal
##' plot(res$lps[dim(res$lps)[1]/2,], 
##'      type='l', 
##'      xlab='spectral value index', 
##'      ylab='spectral value')
##' 
##' @export
'lpsSpectrum' <- function(listOfFiles = NULL, optLogFilePath = NULL,
                          beginTime = 0.0, centerTime = FALSE,
                          endTime = 0.0, resolution = 40.0,
                          fftLength = 0, windowSize = 20.0,
                          windowShift = 5.0, window = 'BLACKMAN',
                          order = 0, preemphasis = -0.95, 
                          deemphasize = TRUE, toFile = TRUE,
                          explicitExt = NULL, outputDirectory = NULL,
                          forceToLog = useWrasspLogger, verbose = TRUE){
  
  ## ########################
  ## a few parameter checks and expand paths
  
  if (is.null(listOfFiles)) {
    stop(paste("listOfFiles is NULL! It has to be a string or vector of file",
               "paths (min length = 1) pointing to valid file(s) to perform",
               "the given analysis function."))
  }
  
  if (is.null(optLogFilePath) && forceToLog){
    stop("optLogFilePath is NULL! -> not logging!")
  }else{
    if(forceToLog){
      optLogFilePath = path.expand(optLogFilePath)  
    }
  }
  
  if(!isAsspWindowType(window)){
    stop("WindowFunction of type '", window,"' is not supported!")
  }
  
  if (!is.null(outputDirectory)) {
    outputDirectory = normalizePath(path.expand(outputDirectory))
    finfo  <- file.info(outputDirectory)
    if (is.na(finfo$isdir))
      if (!dir.create(outputDirectory, recursive=TRUE))
        stop('Unable to create output directory.')
    else if (!finfo$isdir)
      stop(paste(outputDirectory, 'exists but is not a directory.'))
  }
  
  ###########################
  # Pre-process file list
  listOfFiles <- prepareFiles(listOfFiles)
  
  ## #######################
  ## perform analysis
  
  if(length(listOfFiles) == 1 | !verbose){
    pb <- NULL
  }else{
    if(toFile==FALSE){
      stop("length(listOfFiles) is > 1 and toFile=FALSE! toFile=FALSE only permitted for single files.")
    }
    cat('\n  INFO: applying lpsSpectrum to', length(listOfFiles), 'files\n')
    pb <- utils::txtProgressBar(min = 0, max = length(listOfFiles), style = 3)
  }	
  
  externalRes = invisible(.External("performAssp", listOfFiles, 
                                    fname = "spectrum", beginTime = beginTime, 
                                    centerTime = centerTime, endTime = endTime, 
                                    spectrumType = 'LPS',
                                    resolution = resolution, 
                                    fftLength = as.integer(fftLength), windowSize = windowSize, 
                                    windowShift = windowShift, window = window, 
                                    effectiveLength = TRUE, 
                                    order = as.integer(order), preemphasis = preemphasis, 
                                    deemphasize = deemphasize, 
                                    toFile = toFile, explicitExt = explicitExt, 
                                    progressBar = pb, outputDirectory = outputDirectory,
                                    PACKAGE = "superassp"))
  
  
  ## #########################
  ## write options to options log file
  if (forceToLog){
    optionsGivenAsArgs = as.list(match.call(expand.dots = TRUE))
    wrassp.logger(optionsGivenAsArgs[[1]], optionsGivenAsArgs[-1],
                  optLogFilePath, listOfFiles)
  }
  
  ## #########################
  ## return dataObj if length only one file
  
  if(!is.null(pb)){
    close(pb)
  }else{
    return(externalRes)
  }
}
