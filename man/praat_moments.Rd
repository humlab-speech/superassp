% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/praat_ssff.R
\name{praat_moments}
\alias{praat_moments}
\title{Compute spectral moments using Praat}
\usage{
praat_moments(
  listOfFiles,
  beginTime = NULL,
  endTime = NULL,
  windowShift = 5,
  windowSize = 25,
  freqBinSize = 20,
  power = 2,
  maximumFrequency = NULL,
  windowShape = "Gaussian1",
  relativeWidth = 1,
  toFile = TRUE,
  explicitExt = "pmo",
  outputDirectory = NULL,
  verbose = FALSE,
  praat_path = NULL
)
}
\arguments{
\item{listOfFiles}{The full paths of the files to be processed.}

\item{beginTime}{The start time of the portion of a wave that should be included. This argument needs to be either be a single value or a vector of the same length as \code{listOfFiles}.}

\item{endTime}{The end time of the portion of a wave that should be included. Like \code{beginTime} this argument needs to be either be a single value or a vector of the same length as \code{listOfFiles}.}

\item{windowShift}{The time step between (time) analysis windows (in ms).}

\item{windowSize}{The size of the time aligned analysis window (in ms).}

\item{freqBinSize}{The spectral resolution.}

\item{power}{The power to be used when computing the spectral moments. If \code{power=1} the spectral moments will be computed from the absolute spectrum, and if \code{power=2} they will be computed based on the power spectrum.}

\item{maximumFrequency}{The cutoff frequency (in Hz) used when computing the spectrum. Frequencies above this cutoff will not be included when computing spectral moments.}

\item{windowShape}{The window type used for extracting a section of the wave file for analysis. Permitted values are "rectangular", "triangular", "parabolic", "Hanning", "Hamming", "Gaussian1", "Gaussian2", "Gaussian3", "Gaussian4", "Gaussian5", "Kaiser1", and "Kaiser2". See the Praat manual for a descriptio of these window shapes.}

\item{relativeWidth}{The relative width of the windowing function used for extracting part of the sound file.}

\item{toFile}{Should the SSFF signal tracks be store on disk or returned?}

\item{explicitExt}{The signal file extension.}

\item{outputDirectory}{Where should the signal file be stored. If \code{NULL}, the signal file will be stored in the same directory as the wave file.}

\item{verbose}{Produce verbose output?}

\item{praat_path}{The location where the Praat executable is stored.}
}
\value{
An SSFF file (if toFile=FALSE), or nothing.
}
\description{
This function takes a sound file, or an indicated part of a sound file, and
computes the first (Center of gravity), second (Standard deviation), third
(skewness), and fourth (kurtosis) spectral moments at regular intervals. The
moments are stored as separate tracks in an SSFF signal object on disk, or
optionally returned as an object (toFile=FALSE).
}
\details{
By default, the moments are based on the entire spectrum of the sound file, which means that they will depend on the sampling frequency of the file. The user may however choose to include just a portion of the spectrum up to a cutoff frequency if required. This way, the user may opt to limit the analysis to half of the smallest sampling rate for a collection of sound files and thereby produce comparable results.
}
