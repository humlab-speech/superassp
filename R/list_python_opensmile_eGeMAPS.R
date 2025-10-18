#' Compute the eGeMAPS openSMILE feature set
#'
#' This function applies the extended version of the "Geneva Minimalistic Acoustic Parameter Set
#' (eGeMAPS) for Voice Research and Affective Computing" (the *Extended Geneva Minimalistic Standard Parameter Set*, eGeMAPS v02)
#' \insertCite{Eyben.2015.10.1109/taffc.2015.2457417}{superassp} to a portion of a recording.
#' @details The GeMAPS feature set consists of of 88 static acoustic
#'   features resulting from the computation of various functionals over
#'   low-level descriptor features, and is applied by this function using the
#'   openSMILE.
#'   \insertCite{Eyben:2010fq,Jaimes.2013.10.1145/2502081.2502224}{superassp}
#'   acoustic feature extraction library.
#'
#'
#' @param listOfFiles The full path to the sound file.
#' @param beginTime The starting time of the section of the sound files that
#'   should be analysed.
#' @param endTime The end time of the section of the sound files that should be
#'   analysed.
#' @param explicitExt The file extension of the slice file where the results
#'   should be stored.
#'
#' @return A list of 88 acoustic values, with the names as reported by
#'   openSMILE. The extendedacoustic parameter set contains the following compact set of 18 low-level descriptors (LLD), sorted by parameter groups:
#'   
#'   
#' Frequency related parameters:
#' * Pitch, logarithmic f0 on a semitone frequency scale,starting at 27.5 Hz (semitone 0).
#' * Jitter, deviations in individual consecutive f0 period lengths.
#' * Formant 1, 2, and 3 frequency, centre frequency of first, second, and third formant
#' * Formant 1, bandwidth of first formant.Energy/Amplitude related parameters:
#' * Shimmer, difference of the peak amplitudes of consecutive f0 periods.
#' * Loudness, estimate of perceived signal intensity from an auditory spectrum.
#' * Harmonics-to-noise ratio (HNR), relation of energy in harmonic components to energy in noise-like components.
#' * Formant 2-3 bandwidth
#' 
#' Spectral (balance) parameters:
#' * Alpha Ratio, ratio of the summed energy from50-1000 Hz and 1-5 kHz
#' * Hammarberg Index, ratio of the strongest energy peak in the 0-2 kHz region to the strongest peak in the 2–5 kHz region.
#' * Spectral Slope 0-500 Hz and 500-1500 Hz, linear regression slope of the logarithmic power spectrum within the two given bands.
#' * Formant 1, 2, and 3 relative energy, as well as the ratio of the energy of the spectral harmonic peak at the first, second, third formant’s centre frequency to the energy of the spectral peak atF0.
#' * Harmonic difference H1-H2, ratio of energy of the first f0 harmonic (H1) to the energy of the second f0 harmonic (H2).
#' * Harmonic difference H1-A3, ratio of energy of the first f0harmonic (H1) to the energy of the highest harmonic in the third formant range (A3).
#' * MFCC 1-4 Mel-Frequency Cepstral Coefficients 1-4.
#' * Spectral flux difference of the spectra of two consecutive frames.
#' 
#' which are analysed in terms of mean and coefficient of variation, as well as 20th, median (50th), and 80th percentile (pitch and loudness),  
#' the arithmetic mean of the Alpha Ratio, the Hammarberg Index, and the spectral slopes from 0-500 Hz and 500-1500 Hz over all unvoiced segments, and
#'  the equivalent sound level.
#' 
#' Temporal features:
#' * the rate of loudness peaks, i.e., the number of loudness peaks per second,
#' * the mean length and the standard deviation of continuously voiced regions(f0>0),
#' * the mean length and the standard deviation of unvoiced regions (f0 == 0; approximating pauses),
#' * the number of continuous voiced regions per second(pseudo syllable rate).
#' 
#'  Please consult the \insertCite{Eyben.2015.10.1109/taffc.2015.2457417}{superassp} for a description of the
#'   features.
#' @export
#'
#' @references \insertAllCited{}

eGeMAPS<- function(listOfFiles,
                   beginTime=0,
                   endTime=0,
                   explicitExt="ocp"){
  
  origSoundFile <- normalizePath(listOfFiles,mustWork = TRUE)
  if(! file.exists(origSoundFile)){
    stop("Unable to open sound file '",listOfFiles,"'.")
  }
  if(endTime == 0){
    endTime <- NULL
  } 
  if(beginTime == 0){
    beginTime <- NULL
  } 
  
  
  py$soundFile <- reticulate::r_to_py(origSoundFile)
  py$beginTime <- reticulate::r_to_py(beginTime)
  py$endTime <- reticulate::r_to_py(endTime)
  
  
  
  reticulate::py_run_string("import opensmile\
import numpy as np\
import gc\
\
smile = opensmile.Smile(\
    feature_set=opensmile.FeatureSet.eGeMAPSv02,\
    feature_level=opensmile.FeatureLevel.Functionals,\
)\
\
smile_results = smile.process_file(file=soundFile,start=beginTime,end=endTime)\
del soundFile\
gc.collect()")
  
  out <- py$smile_results
  
  return(as.list(out))
}

attr(eGeMAPS,"ext") <-  c("ogs") 
attr(eGeMAPS,"outputType") <-  c("list")
attr(eGeMAPS,"tracks") <- c("F0semitoneFrom27.5Hz_sma3nz_amean", "F0semitoneFrom27.5Hz_sma3nz_stddevNorm", 
                            "F0semitoneFrom27.5Hz_sma3nz_percentile20.0", "F0semitoneFrom27.5Hz_sma3nz_percentile50.0", 
                            "F0semitoneFrom27.5Hz_sma3nz_percentile80.0", "F0semitoneFrom27.5Hz_sma3nz_pctlrange0-2", 
                            "F0semitoneFrom27.5Hz_sma3nz_meanRisingSlope", "F0semitoneFrom27.5Hz_sma3nz_stddevRisingSlope", 
                            "F0semitoneFrom27.5Hz_sma3nz_meanFallingSlope", "F0semitoneFrom27.5Hz_sma3nz_stddevFallingSlope", 
                            "loudness_sma3_amean", "loudness_sma3_stddevNorm", "loudness_sma3_percentile20.0", 
                            "loudness_sma3_percentile50.0", "loudness_sma3_percentile80.0", 
                            "loudness_sma3_pctlrange0-2", "loudness_sma3_meanRisingSlope", 
                            "loudness_sma3_stddevRisingSlope", "loudness_sma3_meanFallingSlope", 
                            "loudness_sma3_stddevFallingSlope", "spectralFlux_sma3_amean", 
                            "spectralFlux_sma3_stddevNorm", "mfcc1_sma3_amean", "mfcc1_sma3_stddevNorm", 
                            "mfcc2_sma3_amean", "mfcc2_sma3_stddevNorm", "mfcc3_sma3_amean", 
                            "mfcc3_sma3_stddevNorm", "mfcc4_sma3_amean", "mfcc4_sma3_stddevNorm", 
                            "jitterLocal_sma3nz_amean", "jitterLocal_sma3nz_stddevNorm", 
                            "shimmerLocaldB_sma3nz_amean", "shimmerLocaldB_sma3nz_stddevNorm", 
                            "HNRdBACF_sma3nz_amean", "HNRdBACF_sma3nz_stddevNorm", "logRelF0-H1-H2_sma3nz_amean", 
                            "logRelF0-H1-H2_sma3nz_stddevNorm", "logRelF0-H1-A3_sma3nz_amean", 
                            "logRelF0-H1-A3_sma3nz_stddevNorm", "F1frequency_sma3nz_amean", 
                            "F1frequency_sma3nz_stddevNorm", "F1bandwidth_sma3nz_amean", 
                            "F1bandwidth_sma3nz_stddevNorm", "F1amplitudeLogRelF0_sma3nz_amean", 
                            "F1amplitudeLogRelF0_sma3nz_stddevNorm", "F2frequency_sma3nz_amean", 
                            "F2frequency_sma3nz_stddevNorm", "F2bandwidth_sma3nz_amean", 
                            "F2bandwidth_sma3nz_stddevNorm", "F2amplitudeLogRelF0_sma3nz_amean", 
                            "F2amplitudeLogRelF0_sma3nz_stddevNorm", "F3frequency_sma3nz_amean", 
                            "F3frequency_sma3nz_stddevNorm", "F3bandwidth_sma3nz_amean", 
                            "F3bandwidth_sma3nz_stddevNorm", "F3amplitudeLogRelF0_sma3nz_amean", 
                            "F3amplitudeLogRelF0_sma3nz_stddevNorm", "alphaRatioV_sma3nz_amean", 
                            "alphaRatioV_sma3nz_stddevNorm", "hammarbergIndexV_sma3nz_amean", 
                            "hammarbergIndexV_sma3nz_stddevNorm", "slopeV0-500_sma3nz_amean", 
                            "slopeV0-500_sma3nz_stddevNorm", "slopeV500-1500_sma3nz_amean", 
                            "slopeV500-1500_sma3nz_stddevNorm", "spectralFluxV_sma3nz_amean", 
                            "spectralFluxV_sma3nz_stddevNorm", "mfcc1V_sma3nz_amean", "mfcc1V_sma3nz_stddevNorm", 
                            "mfcc2V_sma3nz_amean", "mfcc2V_sma3nz_stddevNorm", "mfcc3V_sma3nz_amean", 
                            "mfcc3V_sma3nz_stddevNorm", "mfcc4V_sma3nz_amean", "mfcc4V_sma3nz_stddevNorm", 
                            "alphaRatioUV_sma3nz_amean", "hammarbergIndexUV_sma3nz_amean", 
                            "slopeUV0-500_sma3nz_amean", "slopeUV500-1500_sma3nz_amean", 
                            "spectralFluxUV_sma3nz_amean", "loudnessPeaksPerSec", "VoicedSegmentsPerSec", 
                            "MeanVoicedSegmentLengthSec", "StddevVoicedSegmentLengthSec", 
                            "MeanUnvoicedSegmentLength", "StddevUnvoicedSegmentLength", "equivalentSoundLevel_dBp"
)

attr(eGeMAPS,"nativeFiletypes") <-  NA

