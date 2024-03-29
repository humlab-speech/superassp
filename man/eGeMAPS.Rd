% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/python_slicefunctions.R
\name{eGeMAPS}
\alias{eGeMAPS}
\title{Compute the eGeMAPS openSMILE feature set}
\usage{
eGeMAPS(listOfFiles, beginTime = 0, endTime = 0, explicitExt = "ocp")
}
\arguments{
\item{listOfFiles}{The full path to the sound file.}

\item{beginTime}{The starting time of the section of the sound files that
should be analysed.}

\item{endTime}{The end time of the section of the sound files that should be
analysed.}

\item{explicitExt}{The file extension of the slice file where the results
should be stored.}
}
\value{
A list of 88 acoustic values, with the names as reported by
openSMILE. The extendedacoustic parameter set contains the following compact set of 18 low-level descriptors (LLD), sorted by parameter groups:

Frequency related parameters:
\itemize{
\item Pitch, logarithmic f0 on a semitone frequency scale,starting at 27.5 Hz (semitone 0).
\item Jitter, deviations in individual consecutive f0 period lengths.
\item Formant 1, 2, and 3 frequency, centre frequency of first, second, and third formant
\item Formant 1, bandwidth of first formant.Energy/Amplitude related parameters:
\item Shimmer, difference of the peak amplitudes of consecutive f0 periods.
\item Loudness, estimate of perceived signal intensity from an auditory spectrum.
\item Harmonics-to-noise ratio (HNR), relation of energy in harmonic components to energy in noise-like components.
\item Formant 2-3 bandwidth
}

Spectral (balance) parameters:
\itemize{
\item Alpha Ratio, ratio of the summed energy from50-1000 Hz and 1-5 kHz
\item Hammarberg Index, ratio of the strongest energy peak in the 0-2 kHz region to the strongest peak in the 2–5 kHz region.
\item Spectral Slope 0-500 Hz and 500-1500 Hz, linear regression slope of the logarithmic power spectrum within the two given bands.
\item Formant 1, 2, and 3 relative energy, as well as the ratio of the energy of the spectral harmonic peak at the first, second, third formant’s centre frequency to the energy of the spectral peak atF0.
\item Harmonic difference H1-H2, ratio of energy of the first f0 harmonic (H1) to the energy of the second f0 harmonic (H2).
\item Harmonic difference H1-A3, ratio of energy of the first f0harmonic (H1) to the energy of the highest harmonic in the third formant range (A3).
\item MFCC 1-4 Mel-Frequency Cepstral Coefficients 1-4.
\item Spectral flux difference of the spectra of two consecutive frames.
}

which are analysed in terms of mean and coefficient of variation, as well as 20th, median (50th), and 80th percentile (pitch and loudness),
the arithmetic mean of the Alpha Ratio, the Hammarberg Index, and the spectral slopes from 0-500 Hz and 500-1500 Hz over all unvoiced segments, and
the equivalent sound level.

Temporal features:
\itemize{
\item the rate of loudness peaks, i.e., the number of loudness peaks per second,
\item the mean length and the standard deviation of continuously voiced regions(f0>0),
\item the mean length and the standard deviation of unvoiced regions (f0 == 0; approximating pauses),
\item the number of continuous voiced regions per second(pseudo syllable rate).
}

Please consult the \insertCite{Eyben.2015.10.1109/taffc.2015.2457417}{superassp} for a description of the
features.
}
\description{
This function applies the extended version of the "Geneva Minimalistic Acoustic Parameter Set
(eGeMAPS) for Voice Research and Affective Computing" (the \emph{Extended Geneva Minimalistic Standard Parameter Set}, eGeMAPS v02)
\insertCite{Eyben.2015.10.1109/taffc.2015.2457417}{superassp} to a portion of a recording.
}
\details{
The GeMAPS feature set consists of of 88 static acoustic
features resulting from the computation of various functionals over
low-level descriptor features, and is applied by this function using the
openSMILE.
\insertCite{Eyben:2010fq,Jaimes.2013.10.1145/2502081.2502224}{superassp}
acoustic feature extraction library.
}
\references{
\insertAllCited{}
}
