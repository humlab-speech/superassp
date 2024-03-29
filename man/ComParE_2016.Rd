% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/python_slicefunctions.R
\name{ComParE_2016}
\alias{ComParE_2016}
\title{Compute the ComParE 2016 openSMILE feature set}
\usage{
ComParE_2016(listOfFiles, beginTime = 0, endTime = 0, explicitExt = "ocp")
}
\arguments{
\item{listOfFiles}{The full path to the sound file.}

\item{beginTime}{The starting time of the section of the sound files that should be analysed.}

\item{endTime}{The end time of the section of the sound files that should be analysed.}

\item{explicitExt}{The file extension of the slice file where the results should be stored.}
}
\value{
A list of 6 373 acoustic values, with the names as reported by
openSMILE. Please consult the
\insertCite{10.3389/fpsyg.2013.00292}{superassp} for a description of the
features.
}
\description{
This function applies the "The INTERSPEECH 2016 Computational Paralinguistics
Challenge: Deception, Sincerity & Native Language" ComParE
\insertCite{10.3389/fpsyg.2013.00292}{superassp} to a portion of a recording.
}
\details{
The ComParE feature set consists of of 6 373 static acoustic
features resulting from the computation of various functionals over low-level
descriptor features, and is applied by this function using the openSMILE
\insertCite{Eyben:2010fq,Jaimes.2013.10.1145/2502081.2502224}{superassp}
acoustic feature extraction library.
}
\references{
\insertAllCited{}
}
