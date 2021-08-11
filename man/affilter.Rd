% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrassp_re-export.R
\name{affilter}
\alias{affilter}
\title{Filters the audio signal (From the 'wrassp' package)}
\usage{
affilter(
  listOfFiles = NULL,
  optLogFilePath = NULL,
  highPass = 4000,
  lowPass = 0,
  stopBand = 96,
  transition = 250,
  useIIR = FALSE,
  numIIRsections = 4,
  toFile = TRUE,
  explicitExt = NULL,
  outputDirectory = NULL,
  forceToLog = useWrasspLogger,
  verbose = TRUE
)
}
\arguments{
\item{listOfFiles}{vector of file paths to be processed by function}

\item{optLogFilePath}{path to option log file}

\item{highPass}{= <num>: set the high-pass cut-off frequency to <num> Hz (default: 4000, high-pass filtering is applied)}

\item{lowPass}{= <num>: set the low-pass cut-off frequency to <num> Hz (default: 0, no low-pass filtering)}

\item{stopBand}{= <num>: set the stop-band attenuation to <num> dB (default: 93.0 dB, minimum: 21.0 dB)}

\item{transition}{= <num>: set the width of the transition band to <num> Hz (default: 250.0 Hz)}

\item{useIIR}{switch from the default FIR to IIR filter}

\item{numIIRsections}{= <num>: set the number of 2nd order sections to <num> (default: 4) where each section 
adds 12dB/oct to the slope of the filter}

\item{toFile}{write results to file (for default extension see details section))}

\item{explicitExt}{set if you wish to override the default extension}

\item{outputDirectory}{directory in which output files are stored. Defaults to NULL, i.e. 
the directory of the input files}

\item{forceToLog}{is set by the global package variable useWrasspLogger. This is set
to FALSE by default and should be set to TRUE is logging is desired.}

\item{verbose}{display infos & show progress bar}
}
\value{
nrOfProcessedFiles or if only one file to process return AsspDataObj of that file
}
\description{
affilter function adapted from libassp
}
\details{
Filters the audio signal in <listOfFiles>.
By specifying the high-pass and/or low-pass cut-off
frequency one of four filter characteristics may be
selected as shown in the table below.

\tabular{ccll}{
\strong{hp} \tab \strong{lp} \tab \strong{filter characteristic} \tab \strong{extension}\cr
> 0 \tab 0 \tab high-pass from hp \tab '.hpf'\cr
 0 \tab > 0 \tab low-pass up to lp \tab '.lpf'\cr
> 0 \tab > hp \tab band-pass from hp to lp \tab '.bpf'\cr
> lp \tab > 0 \tab band-stop between lp and hp \tab '.bsf'\cr
}

Please note: per default a high-pass filter from 0 to 4000 Hz is applied.

The Kaiser-window design method is used to compute the
coefficients of a linear-phase FIR filter with unity gain
in the pass-band. The cut-off frequencies (-6 dB points)
of the filters are in the middle of the transition band.
The filtered signal will be written to a file with the
base name of the input file and an extension corresponding
to the filter characteristic (see table). The format of
the output file will be the same as that of the input file.
}
\seealso{
[wrassp::affilter]
}