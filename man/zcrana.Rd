% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrassp_zcrana.R
\name{zcrana}
\alias{zcrana}
\title{Analysis of the averages of the short-term positive and negative zero-crossing rates}
\usage{
zcrana(
  listOfFiles = NULL,
  beginTime = 0,
  centerTime = FALSE,
  endTime = 0,
  windowShift = 5,
  windowSize = 25,
  toFile = TRUE,
  explicitExt = "zcr",
  outputDirectory = NULL,
  assertLossless = NULL,
  logToFile = FALSE,
  convertOverwrites = FALSE,
  keepConverted = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{listOfFiles}{vector of file paths to be processed by function}

\item{beginTime}{the time point (in seconds) of the start of the analysed
interval. A NULL or 0 is interpreted as the start of the signal file.
If a vector of time points is supplied, the length of that vector needs
to correspond with the length of \code{listOfFiles}.}

\item{centerTime}{sets a single-frame analysis time point (in seconds).
Overrides \code{beginTime}, \code{endTime} and \code{windowShift} parameters.}

\item{endTime}{the time point (in seconds) of the end of the analysed
interval. A NULL or 0 is interpreted as the end of the signal file.
If a vector of time points is supplied, the length of that vector needs
to correspond with the length of \code{listOfFiles}.}

\item{windowShift}{the amount of time (in ms) that the analysis window will
be shifted between analysis frames}

\item{windowSize}{the analysis window size (in ms); overrides the effect of
the \code{effectiveLength} parameter}

\item{toFile}{Should the function write the results to a file, with the
(default) file extension (\code{TRUE}) or returned as a list of
\link{AsspDataObj} objects (\code{FALSE})?}

\item{explicitExt}{the file extension will be used when
result files are written (\code{toFile=TRUE}), but the file extension can be
set to something else using this function argument.}

\item{outputDirectory}{directory in which output files are stored. Defaults
to NULL which means that the result file will be stored in the same
directory as the input file.}

\item{assertLossless}{an optional list of file extensions that the user wants to assert
contains losslessly encoded signals data.}

\item{logToFile}{whether to log commands to a separate logfile in the
\code{outputDirectory}. Logging will otherwise be in the function-specific logging
namespace of \link{logger} and will be put wherever this namespace is defined to place its output.
See \link[logger:log_appender]{logger::log_appender} for details.}

\item{verbose}{display verbose information about processing steps taken, as
well as progress bars.}
}
\value{
If \code{toFile} is \code{FALSE}, the function returns a list of \link{AsspDataObj}
objects. If \code{toFile} is \code{TRUE}, the number (integer) of successfully
processed and stored output files is returned.
}
\description{
Analysis of the averages of the short-term positive and
negative zero-crossing rates of the signal in \code{listOfFiles} using the \emph{libassp} C library
\insertCite{s5h}{superassp} function. If \code{toFile} is \code{TRUE}, the results will be written to an output fil in the SSFF binary format, with the
same name as the input file, but with an extension \emph{.zcr} and with a track named 'ZCR\link{Hz}'.
}
\details{
Input signals not in a natively supported file format will be converted
before the autocorrelation functions are computed. The conversion process
will display warnings about input files that are not in known losslessly
encoded formats.

The function is a re-write of the \link[wrassp:zcrana]{wrassp::zcrana} function, but
with media pre-conversion, better checking of preconditions such as the
input file existance, structured logging, and the use of a more modern
framework for user feedback.

The native file type of this function is "wav" files (in "pcm_s16le"
format), SUNs "au", NIST, or CSL formats (kay or NSP extension). Input
signal conversion, when needed, is done by
\href{https://ffmpeg.org/libavcodec.html}{libavcodec} and the excellent \link{av}
wrapper package.
}
\note{
This function is not considered computationally expensive enough to require caching of
results if applied to many signals. However, if the number of signals it will be applied to
is \emph{very} long, then caching of results may be warranted.
}
\examples{
# get path to audio file
path2wav <- list.files(system.file("samples","sustained", package = "superassp"), pattern = glob2rx("a1.wav"), full.names = TRUE)

# calculate zcr values
res <- zcrana(path2wav, toFile=FALSE)

# plot zcr values
plot(seq(0,numRecs.AsspDataObj(res) - 1) / rate.AsspDataObj(res) +
      attr(res, 'startTime'),
    res[["ZCR[Hz]"]],
    type='l',
    xlab='time (s)',
    ylab='Zero Crossing Rates (Hz)')

}
\author{
Raphael Winkelmann

Lasse Bombien

Fredrik Nylén
}
