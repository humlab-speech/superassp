###############TREMOR#################
# This is a console version of the tremor.praat 
# script for Praat[6.1.47] (http://www.praat.org/)
# that extracts 18 measures of vocal tremor from a sound that
# captures a (sustained) phonation.
# Input: wav-file 
# Output: text in csv
######################################
# Author: Markus Brückl (markus.brueckl@tu-berlin.de)
# Copyright: Markus Brückl (2011-2021)
# Adapted for console use by Fredrik Karlsson (fredrik.k.karlsson@umu.se)
# License: GNU GPL > v3 (http://www.gnu.org/licenses/gpl.html)
######################################


######################################
# Global Settings
######################################
form Tremor 3.05
	real StartTime_(s) 0.0
	real EndTime_(s) 0.0
	real SelectionOffset 0.0
	real SelectionLength 2.0
	word WindowType Gaussian1
	real WindowWidth 1.0
   positive Analysis_time_step_(s) 0.015
comment Arguments for initial pitch extraction
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
   optionmenu Output_of_indeterminate_values 2
      option indeterminate values are replaced by zeros
      option indeterminate values are -- undefined --
	sentence Path_of_sound_to_be_analyzed /Users/frkkan96/Desktop/aaa.wav
	sentence Path_and_name_of_result_csv /Users/frkkan96/Desktop/aaa.csv
endform

#Run mode always selected
program_mode = 4

conPath$ = path_of_sound_to_be_analyzed$
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
tremMagThresh = contour_magnitude_threshold

ocFtrem = frequency_tremor_octave_cost
ocAtrem = amplitude_tremor_octave_cost

nan_out = output_of_indeterminate_values


include ./procedures/freqtrem.praat
include ./procedures/amptrem.praat
include ./procedures/readPitchOb.praat
include ./procedures/tremIntIndex.praat
include ./procedures/tremProdSum.praat
include ./procedures/definout.praat
include ./procedures/signlerunconinout.praat
include ./procedures/analysisinout.praat
include ./procedures/runinout.praat
include ./procedures/tremHNR.praat

if program_mode = 1
   call anainout
elsif program_mode = 2
   call dinout 
elsif program_mode = 3
   call rinout
elsif program_mode = 4
 	@cinout: startTime, endTime, selectionOffset, selectionLength, windowType$, windowWidth
endif

