
# TITLE OF THE SCRIPT: ACOUSTIC VOICE QUALITY INDEX (AVQI) v.02.03
# Originally created by Youri Maryn (PhD) and Paul Corthals (PhD)
#
# Modified 2021-07-13 by Fredrik Karlsson to allow application of the procedure to 
# a concatentaion of all sv*.wav and cs*.wav files in a directory.
# The output will also be stored in a PDF in an output directory,
# along with a .csv file containing the results of the procedure in tabulated form.
#
# Form for introduction and/or parameterization

form Acoustic Voice Quality Index v.02.03
comment >>> It is advocated to estimate someone's dysphonia severity in both 
comment continuous speech (i.e., 'cs') and sustained vowel (i.e., 'sv') (Maryn et al.,
comment 2010). This script therefore runs on these two types of recordings, and it
comment important to name these recordings 'cs' and 'sv', respectively.
comment >>> This script automatically (a) searches, extracts and then concatenates
comment the voiced segments of the continuous speech recording to a new sound; (b)
comment concatenates the sustained vowel recording to the new sound, (c) determines
comment the Smoothed Cepstral Peak Prominence, the Shimmer Local, the Shimmer
comment Local dB, the LTAS-slope, the LTAS-tilt and the Harmonics-to-Noise Ratio of
comment the concatenated sound signal, (d) calculates the AVQI-score mostly based
comment on the method of Maryn et al. (2010), and draws the oscillogram, the narrow-
comment band spectrogram with LTAS and the power-cepstrogram with power-
comment cepstrum of the concatenated sound signal to allow further interpretation.
comment >>> For the AVQI to be reliable, it is imperative that the sound recordings
comment are made in an optimal data acquisition conditions.
comment >>> There are two versions in this script: (1) a simple version (only AVQI with
comment data of acoustic measures), and (2) an illustrated version (AVQI with data of 
comment acoustic measures and above-mentioned graphs). 
choice version: 2
button simple
button illustrated
comment >>> Additional information (optional):
sentence name_patient
sentence left_dates_(birth_-_assessment) 
sentence right_dates_(birth_-_assessment) 
comment 
comment Script credits: Youri Maryn (PhD) and Paul Corthals (PhD)
# BEGIN first addition of Fredrik Karlsson 2021-07-13 from original script
comment Modified for batch application by Fredrik Karlsson (PhD)
sentence Input_directory /Users/frkkan96/Documents/src/superassp/tests/signalfiles/AVQI/input
boolean Generate_PDF_files 1
sentence Speaker_ID 1
sentence Output_directory /Users/frkkan96/Documents/src/superassp/tests/signalfiles/AVQI/output
sentence Output_file /Users/frkkan96/Documents/src/superassp/tests/signalfiles/AVQI/output/avqi.csv
endform

# Load all sustained vowels and concatenate them

svLst = Create Strings as file list: "svList", "'input_directory$'/sv*.wav"
noSv = Get number of strings

for sv from 1 to noSv
	select svLst
	currSv$ = Get string: sv
	currSv = Read from file: "'input_directory$'/'currSv$'"

endfor
removeObject: svLst
select all
sv = Concatenate
Rename: "sv"
select all
minus sv
Remove

# Load all continous speech files and concatenate them

csLst = Create Strings as file list: "csList", "'input_directory$'/cs*.wav"
noCS = Get number of strings

for cs from 1 to noCS
	select csLst
	currCS$ = Get string: cs
	currCS = Read from file: "'input_directory$'/'currCS$'"

endfor
removeObject: csLst
select all
minus sv
cs = Concatenate
Rename: "cs"
select all
minus sv
minus cs
Remove
select cs


# END first addition of Fredrik Karlsson 2021-07-13 from original script

Erase all
Select inner viewport... 0.5 7.5 0.5 4.5
Axes... 0 1 0 1
Black
Text special... 0.5 centre 0.6 half Helvetica 12 0 Please wait an instant. Depending on the duration and/or the sample rate of the recorded
Text special... 0.5 centre 0.4 half Helvetica 12 0 sound files, this script takes more or less time to process the sound and search for the AVQI.

# --------------------------------------------------------------------------------------------
# PART 0:
# HIGH-PASS FILTERING OF THE SOUND FILES. 
# --------------------------------------------------------------------------------------------

select Sound cs
Filter (stop Hann band)... 0 34 0.1
Rename... cs
select Sound sv
Filter (stop Hann band)... 0 34 0.1
Rename... sv

# --------------------------------------------------------------------------------------------
# PART 1:
# DETECTION, EXTRACTION AND CONCATENATION OF
# THE VOICED SEGMENTS IN THE RECORDING
# OF CONTINUOUS SPEECH.
# --------------------------------------------------------------------------------------------

select Sound cs
Copy... original
samplingRate = Get sampling frequency
intermediateSamples = Get sampling period
Create Sound... onlyVoice 0 0.001 'samplingRate' 0 
select Sound original
To TextGrid (silences)... 50 0.003 -25 0.1 0.1 silence sounding
select Sound original
plus TextGrid original
Extract intervals where... 1 no "does not contain" silence
Concatenate
select Sound chain
Rename... onlyLoud
globalPower = Get power in air
select TextGrid original
Remove

select Sound onlyLoud
signalEnd = Get end time
windowBorderLeft = Get start time
windowWidth = 0.03
windowBorderRight = windowBorderLeft + windowWidth
globalPower = Get power in air
voicelessThreshold = globalPower*(30/100)

select Sound onlyLoud
extremeRight = signalEnd - windowWidth
while windowBorderRight < extremeRight
	Extract part... 'windowBorderLeft' 'windowBorderRight' Rectangular 1.0 no
	select Sound onlyLoud_part
	partialPower = Get power in air
	if partialPower > voicelessThreshold
		call checkZeros 0
		if (zeroCrossingRate <> undefined) and (zeroCrossingRate < 3000)
			select Sound onlyVoice
			plus Sound onlyLoud_part
			Concatenate
			Rename... onlyVoiceNew
			select Sound onlyVoice
			Remove
			select Sound onlyVoiceNew
			Rename... onlyVoice
		endif
	endif
	select Sound onlyLoud_part
	Remove
	windowBorderLeft = windowBorderLeft + 0.03
	windowBorderRight = windowBorderLeft + 0.03
	select Sound onlyLoud
endwhile
select Sound onlyVoice

procedure checkZeros zeroCrossingRate

	start = 0.0025
	startZero = Get nearest zero crossing... 'start'
	findStart = startZero
	findStartZeroPlusOne = startZero + intermediateSamples
	startZeroPlusOne = Get nearest zero crossing... 'findStartZeroPlusOne'
	zeroCrossings = 0
	strips = 0

	while (findStart < 0.0275) and (findStart <> undefined)
		while startZeroPlusOne = findStart
			findStartZeroPlusOne = findStartZeroPlusOne + intermediateSamples
			startZeroPlusOne = Get nearest zero crossing... 'findStartZeroPlusOne'
		endwhile
		afstand = startZeroPlusOne - startZero
		strips = strips +1
		zeroCrossings = zeroCrossings +1
		findStart = startZeroPlusOne
	endwhile
	zeroCrossingRate = zeroCrossings/afstand
endproc

# --------------------------------------------------------------------------------------------
# PART 2:
# DETERMINATION OF THE SIX ACOUSTIC MEASURES
# AND CALCULATION OF THE ACOUSTIC VOICE QUALITY INDEX.
# --------------------------------------------------------------------------------------------

select Sound sv
durationVowel = Get total duration
durationStart=durationVowel-3
if durationVowel>3
Extract part... durationStart durationVowel rectangular 1 no
Rename... sv2
elsif durationVowel<=3
Copy... sv2
endif

select Sound onlyVoice
durationOnlyVoice = Get total duration
plus Sound sv2
Concatenate
Rename... avqi
durationAll = Get total duration
minimumSPL = Get minimum... 0 0 None
maximumSPL = Get maximum... 0 0 None

# Narrow-band spectrogram and LTAS

To Spectrogram... 0.03 4000 0.002 20 Gaussian
select Sound avqi
To Ltas... 1
minimumSpectrum = Get minimum... 0 4000 None
maximumSpectrum = Get maximum... 0 4000 None

# Power-cepstrogram, Cepstral peak prominence and Smoothed cepstral peak prominence

select Sound avqi
To PowerCepstrogram... 60 0.002 5000 50
cpps = Get CPPS... no 0.01 0.001 60 330 0.05 Parabolic 0.001 0 Straight Robust
To PowerCepstrum (slice)... 0.1
maximumCepstrum = Get peak... 60 330 None

# Slope of the long-term average spectrum

select Sound avqi
To Ltas... 1
slope = Get slope... 0 1000 1000 10000 energy

# Tilt of trendline through the long-term average spectrum

select Ltas avqi
Compute trend line... 1 10000
tilt = Get slope... 0 1000 1000 10000 energy

# Amplitude perturbation measures

select Sound avqi
To PointProcess (periodic, cc)... 50 400
Rename... avqi1
select Sound avqi
plus PointProcess avqi1
percentShimmer = Get shimmer (local)... 0 0 0.0001 0.02 1.3 1.6
shim = percentShimmer*100
shdb = Get shimmer (local_dB)... 0 0 0.0001 0.02 1.3 1.6

# Harmonic-to-noise ratio

select Sound avqi
To Pitch (cc)... 0 75 15 no 0.03 0.45 0.01 0.35 0.14 600
select Sound avqi
plus Pitch avqi
To PointProcess (cc)
Rename... avqi2
select Sound avqi
plus Pitch avqi
plus PointProcess avqi2
voiceReport$ = Voice report... 0 0 75 600 1.3 1.6 0.03 0.45
hnr = extractNumber (voiceReport$, "Mean harmonics-to-noise ratio: ")

# Calculation of the AVQI

avqi = ((3.295-(0.111*cpps)-(0.073*hnr)-(0.213*shim)+(2.789*shdb)-(0.032*slope)+(0.077*tilt))*2.208)+1.797

# --------------------------------------------------------------------------------------------
# PART 3:
# DRAWINGS ALL THE INFORMATION AND THE GRAPHS.
# --------------------------------------------------------------------------------------------

# Title and patient information

Erase all
Solid line
Line width... 1
Black
Helvetica
Select inner viewport... 0 8 0 0.5
Font size... 1
Select inner viewport... 0.5 7.5 0.1 0.15
Axes... 0 1 0 1
Text... 0 Left 0.5 Half Script: Youri Maryn (PhD) and Paul Corthals (PhD)
Font size... 12
Select inner viewport... 0.5 7.5 0 0.5
Axes... 0 1 0 1
Text... 0 Left 0.5 Half ##ACOUSTIC VOICE QUALITY INDEX (AVQI) v.02.03#
Font size... 8
Select inner viewport... 0.5 7.5 0 0.5
Axes... 0 1 0 3
Text... 1 Right 2.3 Half %%'name_patient$'%
Text... 1 Right 1.5 Half %%°'left_dates$'%
Text... 1 Right 0.7 Half %%'right_dates$'%

# Simple version

if version = 1

	# Data

Font size... 10
Select inner viewport... 0.5 7.5 0.5 2
Axes... 0 7 6 0
Text... 0.05 Left 0.5 Half Smoothed cepstral peak prominence (CPPS): ##'cpps:2'#
Text... 0.05 Left 1.5 Half Harmonics-to-noise ratio: ##'hnr:2' dB#
Text... 0.05 Left 2.5 Half Shimmer local: ##'shim:2' \% #
Text... 0.05 Left 3.5 Half Shimmer local dB: ##'shdb:2' dB#
Text... 0.05 Left 4.5 Half Slope of LTAS: ##'slope:2' dB#
Text... 0.05 Left 5.5 Half Tilt of trendline through LTAS: ##'tilt:2' dB#
Select inner viewport... 0.5 3.8 0.5 2
Draw inner box
Font size... 7
Arrow size... 1
Select inner viewport... 4 7.5 1.25 2
Axes... 0 10 1 0
Paint rectangle... green 0 2.91 0 1
Paint rectangle... red 2.91 10 0 1
Draw arrow... avqi 1 avqi 0 
Draw inner box
Marks top every... 1 1 yes yes no
Font size... 16
Select inner viewport... 4 7.5 0.5 1.15
Axes... 0 1 0 1
Text... 0.5 Centre 0.5 Half AVQI: ##'avqi:2'#

	# Copy Praat picture

Select inner viewport... 0.5 7.5 0 2
Copy to clipboard

# Illustrated version

elsif version = 2

	# Oscillogram

Font size... 7
Select inner viewport... 0.5 5 0.5 2.0
select Sound avqi
Draw... 0 0 0 0 no Curve
Draw inner box
One mark left... minimumSPL no yes no 'minimumSPL:2'
One mark left... maximumSPL no yes no 'maximumSPL:2'
Text left... no Sound pressure level (Pa)
One mark bottom... 0 no yes no 0.00
One mark bottom... durationOnlyVoice no no yes
One mark bottom... durationAll no yes no 'durationAll:2'
Text bottom... no Time (s)

	# Narrow-band spectrogram

Select inner viewport... 0.5 5 2.3 3.8
select Spectrogram avqi
Paint... 0 0 0 4000 100 yes 50 6 0 no
Draw inner box
One mark left... 0 no yes no 0
One mark left... 4000 no yes no 4000
Text left... no Frequency (Hz)
One mark bottom... 0 no yes no 0.00
One mark bottom... durationOnlyVoice no no yes
One mark bottom... durationAll no yes no 'durationAll:2'
Text bottom... no Time (s)

	# LTAS

Select inner viewport... 5.4 7.5 2.3 3.8
select Ltas avqi
Draw... 0 4000 minimumSpectrum maximumSpectrum no Curve
Draw inner box
One mark left... minimumSpectrum no yes no 'minimumSpectrum:2'
One mark left... maximumSpectrum no yes no 'maximumSpectrum:2'
Text left... no Sound pressure level (dB/Hz)
One mark bottom... 0 no yes no 0
One mark bottom... 4000 no yes no 4000
Text bottom... no Frequency (Hz)

	# Power-cepstrogram

Select inner viewport... 0.5 5 4.1 5.6
select PowerCepstrogram avqi
Paint: 0, 0, 0, 0, 80, "no", 30, 0, "yes"
Draw inner box
One mark left... 0.00303 no yes no 0.003
One mark left... 0.01667 no yes no 0.017
Text left... no Quefrency (s)
One mark bottom... 0 no yes no 0.00
One mark bottom... durationOnlyVoice no no yes
One mark bottom... durationAll no yes no 'durationAll:2'
Text bottom... no Time (s)

	# Power-cepstrum

Select inner viewport... 5.4 7.5 4.1 5.6
select PowerCepstrum avqi_0_100
Draw... 0.00303 0.01667 0 0 no
Draw tilt line... 0.00303 0.01667 0 0 0.00303 0.01667 Straight Robust
Draw inner box
One mark left... maximumCepstrum no yes no 'maximumCepstrum:2'
Text left... no Amplitude (dB)
One mark bottom... 0.00303 no yes no 0.003
One mark bottom... 0.01667 no yes no 0.017
Text bottom... no Quefrency (s)

	# Data

Font size... 10
Select inner viewport... 0.5 7.5 5.9 7.4
Axes... 0 7 6 0
Text... 0.05 Left 0.5 Half Smoothed cepstral peak prominence (CPPS): ##'cpps:2'#
Text... 0.05 Left 1.5 Half Harmonics-to-noise ratio: ##'hnr:2' dB#
Text... 0.05 Left 2.5 Half Shimmer local: ##'shim:2' \% #
Text... 0.05 Left 3.5 Half Shimmer local dB: ##'shdb:2' dB#
Text... 0.05 Left 4.5 Half Slope of LTAS: ##'slope:2' dB#
Text... 0.05 Left 5.5 Half Tilt of trendline through LTAS: ##'tilt:2' dB#
Select inner viewport... 0.5 3.8 5.9 7.4
Draw inner box
Font size... 7
Arrow size... 1
Select inner viewport... 4 7.5 6.75 7.4
Axes... 0 10 1 0
Paint rectangle... green 0 2.91 0 1
Paint rectangle... red 2.91 10 0 1
Draw arrow... avqi 1 avqi 0 
Draw inner box
Marks top every... 1 1 yes yes no
Font size... 16
Select inner viewport... 4 7.5 5.9 6.65
Axes... 0 1 0 1
Text... 0.5 Centre 0.5 Half AVQI: ##'avqi:2'#

	# Copy Praat picture

Select inner viewport... 0.5 7.5 0 7.4
Copy to clipboard

endif


# BEGIN second addition of Fredrik Karlsson 2021-07-15 from original script

#Now store the results
if generate_PDF_files == 1
	Save as PDF file: "'output_directory$'/'speaker_ID$'.pdf"
endif
outTab = Create Table with column names: "outTab", 1, "SpeakerID CPPS HNR Shim_local Shim_local_DB LTAS_Slope LTAS_Tilt AVQI"
Set string value: 1, "SpeakerID", speaker_ID$
Set numeric value: 1, "CPPS", 'cpps:2'
Set numeric value: 1, "HNR", 'hnr:2'
Set numeric value: 1, "Shim_local", 'shim:2'
Set numeric value: 1, "Shim_local_DB", 'shdb:2'
Set numeric value: 1, "LTAS_Slope", 'slope:2'
Set numeric value: 1, "LTAS_Tilt", 'tilt:2'
Set numeric value: 1, "AVQI", 'avqi:2'
Save as comma-separated file: output_file$

# END second addition of Fredrik Karlsson 2021-07-15 from original script


# Remove intermediate objects

select all
minus Sound cs
minus Sound sv
Remove
