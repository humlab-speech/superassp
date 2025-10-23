##' afdiff function adapted from libassp (Rcpp-optimized)
##'
##' Computes the first difference of the signal in the audio-
##' formatted file(s) <listOfFiles>. The differentiated signal will
##' be written to a file with the base name of the input file
##' and an extension consisting of '.d', followed by the
##' extension of the input file. The format of the output file
##' will be the same as that of the input file.
##' Differentiation can improve results on F0 analysis of e.g.
##' EGG signals because it removes a DC offset, attenuates
##' very low frequency components - and in the case of central
##' differentiation also very high ones - and enhances the
##' moment of glottal closure.
##'
##' @details This version includes Rcpp optimizations for improved
##' performance on large batches of files.
##'
##' @title afdiff
##' @param listOfFiles vector of file paths to be processed by function
##' @param computeBackwardDifference compute backward difference (s'\[n\] = s\[n\] - s\[n-1\]) (default: forward difference s'\[n\] = s\[n+1\] - s\[n\])
##' @param computeCentralDifference compute central/interpolated/3-point difference
##' @param channel = <num>: for multi-channel input files: extract and differentiate channel <num> (1 <= <num> <= 8  default: channel 1)
##' @param toFile write results to file (default extension is .d+(extensionsOfAudioFile))
##' @param explicitExt set if you wish to override the default extension
##' @param outputDirectory directory in which output files are stored. Defaults to NULL, i.e.
##' the directory of the input files
##' @param verbose display infos & show progress bar
##' @return nrOfProcessedFiles or if only one file to process return AsspDataObj of that file
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
##' # compute the first forward difference of the signal
##' res <- afdiff(path2wav, toFile=FALSE)
##'
##' # plot samples
##' # (only plot every 10th element to accelerate plotting)
##' plot(seq(0,numRecs.AsspDataObj(res) - 1, 10) / rate.AsspDataObj(res),
##'      res$audio[c(TRUE, rep(FALSE,9))],
##'      type='l',
##'      xlab='time (s)',
##'      ylab='Audio samples')
##'
##' @export
'afdiff' <- function(listOfFiles = NULL,
                     computeBackwardDifference = FALSE,
                     computeCentralDifference = FALSE,
                     channel = 1,
                     toFile = TRUE,
                     explicitExt=NULL,
                     outputDirectory = NULL,
                     verbose = TRUE){

  ###########################
  # a few parameter checks and expand paths

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
      cli::cli_abort("length(listOfFiles) is > 1 and toFile=FALSE! toFile=FALSE only permitted for single files.")
    }
    cli::cli_inform("Applying {.fun afdiff} to {cli::no(n_files)} file{?s}")
  }

  # Use unified load-and-process helper (works for all file formats)
  nativeFiletypes <- c("wav", "au", "kay", "nist", "nsp")
  result <- processMediaFiles_LoadAndProcess(
    listOfFiles = listOfFiles,
    beginTime = rep(0.0, n_files),
    endTime = rep(0.0, n_files),
    nativeFiletypes = nativeFiletypes,
    fname = "afdiff",
    toFile = toFile,
    verbose = FALSE,  # Already handled above
    computeBackwardDifference = computeBackwardDifference,
    channel = as.integer(channel),
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
