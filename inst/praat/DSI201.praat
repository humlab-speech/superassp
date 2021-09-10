form Dysphonia Severity Index in Praat v2.01
	comment >>> The original Dysphonia Severity Index (i.e., DSI) is developed 
	comment by Wuyts et al. (2000). It is a multivariate measure of the degree of 
	comment of dysphonia and consists of a set of four weighted measures: 
	comment maximum phonation time, highest fundamental frequency (from the comment Voice Range Profile of the KayPentax CSL), softest intensity (also from the 
	comment Voice Range Profile of the KayPentax CSL) and percent jitter (from the 
	comment Multi-Dimensional Voice Program of the KayPentax CSL). 
	comment >>> The output and the measures of this script are not equal to but highly 
	comment resemble the output and the measures of the original DSI: maximum 
	comment phonation time (i.e., 'mpt'), highest fundamental frequency or f0-high (i.e. 
	comment 'fh'), softest intensity of voiced speech or I-low (i.e., 'il') and jitter 
	comment ppq5 (i.e., 'ppq'). Their weighted combination is strongly correlated to 
	comment the original DSI-data. 
	comment >>> Make sure that the intensity measurements are calibrated. Does the 
	comment calibration method necessitate the implementation of a calibration factor?
	choice choose: 2
		button no
		button yes
	sentence Calibration_factor 1*self+10
	positive Maximum_phonation_time_(s) 2
	comment >>> Additional information (optional):
	sentence name_patient
	sentence left_dates_(birth_-_assessment)
	sentence right_dates_(birth_-_assessment)
endform


mpt = maximum_phonation_time
Erase all
Select inner viewport... 0.5 7.5 0.5 4.5
Axes... 0 1 0 1 
Black
Text special... 0.5 centre 0.6 half Helvetica 12 0 Please wait an instant. Depending on the duration and/or the sample rate of the recorded
Text special... 0.5 centre 0.4 half Helvetica 120 sound files, this script takes more or less time to process the sounds and calculate the DSI-2.

# Part 2 of Praat-script: determination of IMIN. The sound recording
# with the softest phonations should be named "Sound im".

select Sound im
To Pitch (cc)... 0 70 15 no 0.03 0.8 0.01 0.35 0.14 600

select Sound im
plus Pitch im
To PointProcess (cc)

select PointProcess im_im
To TextGrid (vuv)... 0.02 0.01
select Sound im
plus TextGrid im_im

Extract intervals where... 1 no "is equal to" V
Concatenate

if choose = 1 
	To Intensity... 60 0.0 yes 
elsif choose = 2 
	To Intensity... 60 0.0 yes 
	Formula... 'calibration_factor$'
endif

minimumIntensity = Get minimum... 0 0 none

# Part 3 of Praat-script: determination of FOHIGH. The sound recording 
# with the highest phonations should be named "Sound fh".

select Sound fh
To Pitch (cc)... 0 70 15 no 0.03 0.8 0.01 0.35 0.14 1300
maximumF0 = Get maximum... 0 0 Hertz none

# Part 4 of Praat-script: determination of JITTER PP05. The sound recording
# with sustained (a:] phonation should be named "Sound ppq".

select Sound ppq
duration_Vowel = Get total duration
durationStart = duration_Vowel - 3 

if duration_Vowel > 3
	Extract part... duration_Start duration_Vowel rectangular 1 no
	Rename... ppq2
elsif duration_Vowel <=3
	Copy... ppq2
endif

To Pitch... 0 70 600

select Sound ppq2
plus Pitch ppq2
To PointProcess (cc)

select Sound ppq2
plus Pitch ppq2
plus PointProcess ppq2_ppq2

voiceReport$ = Voice report... 0 0 70 600 1.3 1.6 0.03 0.45

jitterPpq5Pre = extractNumber (voiceReport$, "Jitter (ppq5): ") 
jitterPpq5 = jitterPpq5Pre*100

# Part 5 of Praat-script: calculation of DSIbeta.
dsi2 = 1.127+ 0.164*mpt - 0.038*minimumIntensity + 0.0053*maximumF0 - 5.30*jitterPpq5

# Part 6 of Praat-script: drawing all the information and relevant graphs. # To insert this output in another program, just use to 'paste'-function in e.g. # the text editor after this script is terminated.
Erase all
Solid line
Line width... 1 
Black
Helvetica
Font size... 1 
Select inner viewport... 0.5 7.5 0.1 0.15
Axes... 0 1 0 1

Text... 0 left 0.5 half Script: Youri Maryn, PhD 12
Select inner viewport... 0.5 7.5 0 0.5 
Axes... 0 1 0 1 
Text... 0 left 0.5 half 

##DYSPHONIA SEVERITY INDEX (DSI) IN PRAAT, v.02.01# 

Font size... 8 
Select inner viewport... 0.5 7.5 0 0.5 
Axes... 0 1 0 3
Text... 1 right 2.3 half 'name_patients'
Text... 1 right 1.5 half 880 'left_dates$' 
Text... 1 right 0.7 half is right_dates $'s
Font size... 10 
Select inner viewport... 0.5 7.5 0.5 2 
Axes... 0 7 4 0
Text... 0.05 left 0.5 half Maximum phonation time: ## 'mpt:2' s#
Text... 0.05 left 1.5 half Softest intensity of voiced speech: ## 'minimumIntensity:2' dB# 
Text... 0.05 left 2.5 half Maximum fundamental frequency: ## 'maximumF0:2' Hz# 
Text... 0.05 left 3.5 half Jitter ppq5: ## 'jitterPpq5:2' \% #
Select inner viewport... 0.5 3.8 0.5 2 
Draw inner box
Font size... 8
Arrow size... 1 
Select inner viewport... 4 7.5 1.25 2 
Axes... -10 10 1 0
Paint rectangle... green 1.6 5 0 1 
Paint rectangle... red -5 1.6 0 1 
Draw arrow... dsi2 1 dsi2 0 
One mark top... -5 yes yes no 
One mark top... -4 yes yes no 
One mark top... -3 yes yes no 
One mark top... -2 yes yes no 
One mark top... -1 yes yes no 

One mark top... 1 yes yes no 
One mark top... 2 yes yes no 
One mark top... 3 yes yes no 
One mark top... 4 yes yes no 
One mark top... 5 yes yes no 
Font size... 16

Select inner viewport... 4 7.5 0.5 1.15 
Axes... 0 1 0 1 
Text... 0.5 centre 0.5 half DSI: ##'dsi2:2' 
Select inner viewport... 4.875 6.625 1.25 2
Draw inner box
Select inner viewport... 0.5 7.5 0 2 
Copy to clipboard

select all
minus Sound im
minus Sound fh
minus Sound ppq
Remove
