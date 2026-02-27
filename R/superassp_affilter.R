##' affilter function adapted from libassp (Rcpp-optimized)
##'
##' Filters the audio signal in `listOfFiles`.
##' By specifying the high-pass and/or low-pass cut-off
##' frequency one of four filter characteristics may be
##' selected as shown in the table below. Per default a high-pass filter from 0 to 4000 Hz is applied.
##' The Kaiser-window design method is used to compute the
##' coefficients of a linear-phase FIR filter with unity gain
##' in the pass-band. The cut-off frequencies (-6 dB points)
##'
##' of the filters are in the middle of the transition band.
##' The filtered signal will be written to a file with the
##' base name of the input file and an extension corresponding
##' to the filter characteristic (see table). The format of
##' the output file will be the same as that of the input file.
##'
##' \tabular{ccll}{
##'   \strong{hp}   \tab \strong{lp} \tab \strong{filter characteristic} \tab \strong{extension}\cr
##'    '>0' \tab 0    \tab high-pass from hp \tab '.hpf'\cr
##'     0 \tab > 0  \tab low-pass up to lp \tab '.lpf'\cr
##'    '>0' \tab '>hp' \tab band-pass from hp to lp \tab '.bpf'\cr
##'    '>lp' \tab '>0' \tab band-stop between lp and hp \tab '.bsf'\cr
##' }
##'
##' @details This version includes Rcpp optimizations for improved
##' performance on large batches of files.
##'
##' @param listOfFiles vector of file paths to be processed by function
##' @param highPass = <num>: set the high-pass cut-off frequency to <num> Hz (default: 4000, high-pass filtering is applied)
##' @param lowPass = <num>: set the low-pass cut-off frequency to <num> Hz (default: 0, no low-pass filtering)
##' @param stopBand = <num>: set the stop-band attenuation to <num> dB (default: 93.0 dB, minimum: 21.0 dB)
##' @param transition = <num>: set the width of the transition band to <num> Hz (default: 250.0 Hz)
##' @param useIIR switch from the default FIR to IIR filter
##' @param numIIRsections = <num>: set the number of 2nd order sections to <num> (default: 4) where each section
##' adds 12dB/oct to the slope of the filter
##' @param toFile write results to file (for default extension see details section))
##' @param explicitExt set if you wish to override the default extension
##' @param outputDirectory directory in which output files are stored. Defaults to NULL, i.e.
##' the directory of the input files
##' @param verbose display infos & show progress bar
##'
##' @return nrOfProcessedFiles or if only one file to process return AsspDataObj of that file
##'
##' @author Raphael Winkelmann
##' @author Lasse Bombien
##' @useDynLib superassp, .registration = TRUE
##' @importFrom Rcpp sourceCpp
##' @examples
##' # get path to audio file
##' path2wav <- list.files(system.file("extdata", package = "wrassp"),
##'                        pattern = glob2rx("*.wav"),
##'                        full.names = TRUE)[1]
##'
##' # band-pass filter signal between 4000 and 5000 Hz
##' res <- affilter(path2wav, highPass=4000, lowPass=5000, toFile=FALSE)
##'
##' # plot samples
##' # (only plot every 10th element to accelerate plotting)
##' plot(seq(0,numRecs.AsspDataObj(res) - 1, 10) / rate.AsspDataObj(res),
##'      res$audio[c(TRUE, rep(FALSE,9))],
##'      type='l',
##'      xlab='time (s)',
##'      ylab='Audio samples')
##'
'affilter' <- function(listOfFiles = NULL,
                       highPass = 4000,
                       lowPass = 0,
                       stopBand = 96,
                       transition = 250,
                       useIIR = FALSE,
                       numIIRsections = 4,
                       toFile = TRUE,
                       explicitExt = NULL,
                       outputDirectory = NULL,
                       verbose = TRUE){

  ###########################
  ### a few parameter checks and expand paths

  if (is.null(listOfFiles)) {
    cli::cli_abort("listOfFiles is NULL! It has to be a string or vector of file paths.")
  }

  if (!is.null(outputDirectory)) {
    outputDirectory = normalizePath(path.expand(outputDirectory))
    finfo  <- file.info(outputDirectory)
    if (is.na(finfo$isdir)) {
      if (!dir.create(outputDirectory, recursive=TRUE)) {
        cli::cli_abort('Unable to create output directory.')
      }
    } else if (!finfo$isdir) {
      cli::cli_abort(paste(outputDirectory, 'exists but is not a directory.'))
    }
  }

  ###########################
  # Use unified memory-based processing

  n_files <- length(listOfFiles)

  if(verbose && n_files > 1){
    if(toFile==FALSE){
      cli::cli_abort("length(listOfFiles) is > 1 and ToFile=FALSE! ToFile=FALSE only permitted for single files.")
    }
    cli::cli_inform("Applying {.fun affilter} to {cli::no(n_files)} file{?s}")
  }

  # Use unified load-and-process helper (works for all file formats)
  nativeFiletypes <- c("wav", "au", "kay", "nist", "nsp")
  result <- processMediaFiles_LoadAndProcess(
    listOfFiles = listOfFiles,
    beginTime = rep(0.0, n_files),
    endTime = rep(0.0, n_files),
    nativeFiletypes = nativeFiletypes,
    fname = "affilter",
    toFile = toFile,
    verbose = FALSE,  # Already handled above
    highPass = highPass,
    lowPass = lowPass,
    stopBand = stopBand,
    transition = transition,
    useIIR = useIIR,
    numIIRsections = as.integer(numIIRsections),
    explicitExt = explicitExt,
    outputDirectory = outputDirectory
  )

  externalRes <- result$externalRes

  #############################
  # return dataObj if length only one file

  if(n_files == 1 && !toFile){
    return(externalRes[[1]])
  }

  invisible(NULL)
}
