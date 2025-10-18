#' Compute the GeMAPS openSMILE feature set
#'
#' This function applies the "The Geneva Minimalistic Acoustic Parameter Set
#' (GeMAPS) for Voice Research and Affective Computing" (the *Geneva Minimalistic Standard Parameter Set*, GeMAPS v0.1b)
#' \insertCite{Eyben.2015.10.1109/taffc.2015.2457417}{superassp} to a portion of a recording.
#' @details The GeMAPS feature set consists of of 62 static acoustic
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
#' @return A list of 62 acoustic values, with the names as reported by
#'   openSMILE. The minimalistic acoustic parameter set contains the following compact set of 18 low-level descriptors (LLD), sorted by parameter groups:
#'   
#' Frequency related parameters:
#' * Pitch, logarithmic f0 on a semitone frequency scale,starting at 27.5 Hz (semitone 0).
#' * Jitter, deviations in individual consecutive f0 period lengths.
#' * Formant 1, 2, and 3 frequency, centre frequency of first, second, and third formant
#' * Formant 1, bandwidth of first formant.Energy/Amplitude related parameters:
#' * Shimmer, difference of the peak amplitudes of consecutive f0 periods.
#' * Loudness, estimate of perceived signal intensity from an auditory spectrum.
#' * Harmonics-to-noise ratio (HNR), relation of energy in harmonic components to energy in noise-like components.
#' 
#' Spectral (balance) parameters:
#' * Alpha Ratio, ratio of the summed energy from50-1000 Hz and 1-5 kHz
#' * Hammarberg Index, ratio of the strongest energy peak in the 0-2 kHz region to the strongest peak in the 2–5 kHz region.
#' * Spectral Slope 0-500 Hz and 500-1500 Hz, linear regression slope of the logarithmic power spectrum within the two given bands.
#' * Formant 1, 2, and 3 relative energy, as well as the ratio of the energy of the spectral harmonic peak at the first, second, third formant’s centre frequency to the energy of the spectral peak atF0.
#' * Harmonic difference H1-H2, ratio of energy of the first f0 harmonic (H1) to the energy of the second f0 harmonic (H2).
#' * Harmonic difference H1-A3, ratio of energy of the first f0harmonic (H1) to the energy of the highest harmonic in the third formant range (A3).
#' 
#' which are analysed in terms of mean and coefficient of variation, as well as 20th, median (50th), and 80th percentile (pitch and loudness), and 
#' the arithmetic mean of the Alpha Ratio, the Hammarberg Index, and the spectral slopes from 0-500 Hz and 500-1500 Hz over all unvoiced segments.
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

GeMAPS <- function(listOfFiles,
                         beginTime=0,
                         endTime=0,
                         explicitExt="ocp"){

  origSoundFile <- normalizePath(listOfFiles,mustWork = TRUE)
  if(! file.exists(origSoundFile)){
    stop("Unable to open sound file '",listOfFiles,"'.")
  }

  # Convert time parameters (openSMILE uses seconds, av uses seconds)
  bt <- if(beginTime == 0) 0 else beginTime
  et <- if(endTime == 0) NULL else endTime

  # Load audio with av → convert to numpy (MEMORY-BASED, no disk I/O!)
  audio_result <- av_load_for_python(
    origSoundFile,
    start_time = bt,
    end_time = et
  )

  # Pass numpy array and sample rate to Python
  py <- reticulate::import_main()
  py$audio_np <- audio_result$audio_np
  py$fs <- audio_result$sample_rate

  # Process with openSMILE (using audio signal instead of file)
  reticulate::py_run_string("import opensmile
import numpy as np
import gc

smile = opensmile.Smile(
    feature_set=opensmile.FeatureSet.GeMAPSv01b,
    feature_level=opensmile.FeatureLevel.Functionals,
)
# openSMILE can process signal directly (no file I/O!)
smile_results = smile.process_signal(signal=audio_np, sampling_rate=fs)
del audio_np
gc.collect()")

  out <- py$smile_results

  return(as.list(out))

}

attr(GeMAPS,"ext") <-  c("oge") 
attr(GeMAPS,"outputType") <-  c("list")
attr(GeMAPS,"tracks") <- c("F0semitoneFrom27.5Hz_sma3nz_amean", "F0semitoneFrom27.5Hz_sma3nz_stddevNorm", 
                           "F0semitoneFrom27.5Hz_sma3nz_percentile20.0", "F0semitoneFrom27.5Hz_sma3nz_percentile50.0", 
                           "F0semitoneFrom27.5Hz_sma3nz_percentile80.0", "F0semitoneFrom27.5Hz_sma3nz_pctlrange0-2", 
                           "F0semitoneFrom27.5Hz_sma3nz_meanRisingSlope", "F0semitoneFrom27.5Hz_sma3nz_stddevRisingSlope", 
                           "F0semitoneFrom27.5Hz_sma3nz_meanFallingSlope", "F0semitoneFrom27.5Hz_sma3nz_stddevFallingSlope", 
                           "loudness_sma3_amean", "loudness_sma3_stddevNorm", "loudness_sma3_percentile20.0", 
                           "loudness_sma3_percentile50.0", "loudness_sma3_percentile80.0", 
                           "loudness_sma3_pctlrange0-2", "loudness_sma3_meanRisingSlope", 
                           "loudness_sma3_stddevRisingSlope", "loudness_sma3_meanFallingSlope", 
                           "loudness_sma3_stddevFallingSlope", "jitterLocal_sma3nz_amean", 
                           "jitterLocal_sma3nz_stddevNorm", "shimmerLocaldB_sma3nz_amean", 
                           "shimmerLocaldB_sma3nz_stddevNorm", "HNRdBACF_sma3nz_amean", 
                           "HNRdBACF_sma3nz_stddevNorm", "logRelF0-H1-H2_sma3nz_amean", 
                           "logRelF0-H1-H2_sma3nz_stddevNorm", "logRelF0-H1-A3_sma3nz_amean", 
                           "logRelF0-H1-A3_sma3nz_stddevNorm", "F1frequency_sma3nz_amean", 
                           "F1frequency_sma3nz_stddevNorm", "F1bandwidth_sma3nz_amean", 
                           "F1bandwidth_sma3nz_stddevNorm", "F1amplitudeLogRelF0_sma3nz_amean", 
                           "F1amplitudeLogRelF0_sma3nz_stddevNorm", "F2frequency_sma3nz_amean", 
                           "F2frequency_sma3nz_stddevNorm", "F2amplitudeLogRelF0_sma3nz_amean", 
                           "F2amplitudeLogRelF0_sma3nz_stddevNorm", "F3frequency_sma3nz_amean", 
                           "F3frequency_sma3nz_stddevNorm", "F3amplitudeLogRelF0_sma3nz_amean", 
                           "F3amplitudeLogRelF0_sma3nz_stddevNorm", "alphaRatioV_sma3nz_amean", 
                           "alphaRatioV_sma3nz_stddevNorm", "hammarbergIndexV_sma3nz_amean", 
                           "hammarbergIndexV_sma3nz_stddevNorm", "slopeV0-500_sma3nz_amean", 
                           "slopeV0-500_sma3nz_stddevNorm", "slopeV500-1500_sma3nz_amean", 
                           "slopeV500-1500_sma3nz_stddevNorm", "alphaRatioUV_sma3nz_amean", 
                           "hammarbergIndexUV_sma3nz_amean", "slopeUV0-500_sma3nz_amean", 
                           "slopeUV500-1500_sma3nz_amean", "loudnessPeaksPerSec", "VoicedSegmentsPerSec", 
                           "MeanVoicedSegmentLengthSec", "StddevVoicedSegmentLengthSec", 
                           "MeanUnvoicedSegmentLength", "StddevUnvoicedSegmentLength")

