###############TREMOR#################
# tremor.praat is a Praat[6.0.14] script (http://www.praat.org/)
# that extracts 14 measures of vocal tremor from a wav-file that
# captures a (sustained) phonation.
# Input: wav-file(s)
# Output: 
# (in run mode: a tab-separated csv-file containing)
# (1) frequency contour magnitude (FCoM),
# (2) frequency tremor frequency (FTrF),
# (3) frequency tremor cyclicality (FTrC),
# (4) frequency tremor intensity index (FTrI),
# (5) frequency tremor power index (FTrP),
# (6) frequency tremor cyclicality intensity product (FTrCIP),
# (7) frequency tremor product sum (FTrPS),
# (8) amplitude contour magnitude (ACoM),
# (9) amplitude tremor frequency (ATrF),
# (10) amplitude tremor cyclicality (ATrC),
# (11) amplitude tremor intensity index (ATrI),
# (12) amplitude tremor power index (ATrP),
# (13) amplitude tremor cyclicality intensity product (ATrCIP), and
# (14) amplitude tremor product sum (ATrPS).
######################################
# Documentation, help:
# Markus Brückl (2012): Vocal Tremor Measurement Based on Autocorrelation 
# of Contours. In: Proceedings of the ISCA Conference Interspeech '12, 
# Portland (OR).
######################################
# Author: Markus Brückl (markus.brueckl@tu-berlin.de)
# License: GNU GPL > v3 (http://www.gnu.org/licenses/gpl.html)
######################################
# 1.00: the version used in Brückl (2011)
# 2.01: the version presented in 2012 at Interspeech
# new in...
# 2.02: modularization of the script into separate files
# 2.03: second way to extract an AmplitudeTier from Sound & PointProcess: RMS per 
#       pitch period
# 2.04: removed a bug in amptrem.praat causing a (the more) raised "zero"-level in 
#       amplitude contour (the lower the mean sound intensity was) and therefore 
#       raised ATrI and ATrP values
# 2.05: removed a bug causing the script to stop if the sound is considered voiceless
#       at the beginnig or at the end -- e.g. because of wrong pitch range
# 2.06: improved the resampling at a constant rate (time step) of the amplitude 
#       contour (per period)
# 2.07: output of measures of the contours' (mean) magnitudes (technically speaking: 
#       the contours' mean (relative) distance to zero or as Praat names it in Pitch 
#       objects: (a frame's) intensity)
# 2.07: output of measures for tremor cyclicality (measures of tremor 
#       periodicity, technically speaking: the contours' auto-correlation coefficients, 
#       or as Praat names it in Pitch objects: the frequency's "strength")
# 3.01: invention of new measures that combine cyclicality and intensity
#       A version that takes just a single file and allow easy calling from the console 
#       created by Fredrik Karlsson
######################################



######################################
# Global Settings
######################################
form Tremor 3.01
   positive Analysis_time_step_(s) 0.015
comment Arguments for mandatory pitch extraction
   positive Minimal_pitch_(Hz) 60
   positive Maximal_pitch_(Hz) 350
   positive Silence_threshold 0.03
   positive Voicing_threshold 0.3
   positive Octave_cost 0.01
   positive Octave-jump_cost 0.35
   positive Voiced_/_unvoiced_cost 0.14
comment Arguments for tremor extraction from contours
   optionmenu Amplitude_extraction_method 2
      option Integral [RMS per pitch period]
      option Envelope [To AmplitudeTier (period)]
   positive Minimal_tremor_frequency_(Hz) 1.5
   positive Maximal_tremor_frequency_(Hz) 15
   positive Contour_magnitude_threshold 0.01
   positive Tremor_cyclicality_threshold 0.15
   positive Frequency_tremor_octave_cost 0.01
   positive Amplitude_tremor_octave_cost 0.01
	sentence Path_of_sound_to_be_analyzed /Users/frkkan96/Desktop/aaa.wav
	sentence Path_and_name_of_result_csv /Users/frkkan96/Desktop/aaa.csv
endform
#Run mode always selected
mode = 2

ts = analysis_time_step; [s]

minPi = minimal_pitch; [Hz]
maxPi = maximal_pitch; [Hz]
silThresh = silence_threshold
voiThresh = voicing_threshold
ocCo = octave_cost
ocjCo = 'octave-jump_cost'
vuvCo = 'voiced_/_unvoiced_cost'

minTr = minimal_tremor_frequency; [Hz]
maxTr = maximal_tremor_frequency; [Hz]
tremthresh = tremor_cyclicality_threshold
# equals "the strength of the unvoiced candidate, 
# relative to the maximum possible autocorrelation." (cf. Praat manual)
# The max. possible autocorrelation varies (e.g. between different sounds)
# and could be around 0.9.
# tremthresh = tremthresh * 0.9
# same for the 3 other thresholds
tremMagThresh = contour_magnitude_threshold

ocFtrem = frequency_tremor_octave_cost
ocAtrem = amplitude_tremor_octave_cost


include ./procedures/amptrem.praat
include ./procedures/freqtrem.praat
include ./procedures/analysisinout.praat
include ./procedures/singleruninout.praat
include ./procedures/getCyclicality.praat
include ./procedures/tremIntIndex.praat
include ./procedures/tremProdSum.praat


@singleruninout: path_of_sound_to_be_analyzed$, path_and_name_of_result_csv$
