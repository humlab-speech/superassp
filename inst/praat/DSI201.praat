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
	boolean Apply_calibration 1
	sentence Calibration_factor 1*self+10
	#positive Maximum_phonation_time_(s) 2 #This are now computed from extracted performances
	comment >>> Additional information (optional):
	sentence Name_of_patient Fredrik Karlsson
	sentence Date_of_birth 1975-12-31
	sentence Assessment_date 2021-12-31
	# BEGIN first addition of Fredrik Karlsson 2021-09-10 from original script
	comment Modified for batch application by Fredrik Karlsson (PhD)
	sentence Input_directory ../../tests/signalfiles/DSI/input
	#/Users/frkkan96/Documents/src/superassp/tests/signalfiles/DSI/input
	boolean Generate_PDF_files 1
	sentence Speaker_ID 1
	sentence Output_directory /Users/frkkan96/Documents/src/superassp/tests/signalfiles/DSI/output
	sentence Output_file /Users/frkkan96/Documents/src/superassp/tests/signalfiles/DSI/output/dsi.csv
endform


# Make a clean workspace

select all
nOSelected = numberOfSelected ()

if nOSelected > 0
	Remove
endif


# Load all samples are maximally prolonged vowels for the computation of MPT. 

mptLst = Create Strings as file list: "mptList", "'input_directory$'/mpt*.wav"
noMPT = Get number of strings

mpt = 0

for currMPT from 1 to noMPT
	select mptLst
	currMPT$ = Get string: currMPT
	Open long sound file: "'input_directory$'/'currMPT$'"
	curr = Get total duration
	if curr > mpt
		mpt = curr
	endif
endfor

removeObject: mptLst

# Load all samples that should be analysed the softest phonations for determination of IMIN. 

imLst = Create Strings as file list: "imList", "'input_directory$'/im*.wav"
noIm = Get number of strings

imCount = 1

while imCount <= noIm
	select imLst
	currIm$ = Get string: imCount
	if imCount == 1
		outIm = Read from file: "'input_directory$'/'currIm$'"
		Rename: "im"
	else 
		currIm = Read from file: "'input_directory$'/'currIm$'"
		select outIm
		plus currIm
		Concatenate
		Rename: "im"
		removeObject: currIm
	endif
	imCount = imCount + 1
endwhile

removeObject: imLst


# Load all samples that should be analysed the highest pitch for determination of FOHIGH.

fhLst = Create Strings as file list: "fhList", "'input_directory$'/fh*.wav"
noFH = Get number of strings

fhCount = 1

while fhCount <= noFH
	select fhLst
	currFH$ = Get string: fhCount
	if fhCount == 1
		outFH = Read from file: "'input_directory$'/'currFH$'"
		Rename: "fh"
	else 
		currFH = Read from file: "'input_directory$'/'currFH$'"
		select outFH
		plus currFH
		Concatenate
		Rename: "fh"
		removeObject: currFH
	endif
	fhCount = fhCount + 1
endwhile

removeObject: fhLst


# Load all sustained vowels samples to determine JITTER PP05.

ppqLst = Create Strings as file list: "ppqList", "'input_directory$'/ppq*.wav"
noPPQ = Get number of strings

ppqCount = 1

while ppqCount <= noPPQ
	select ppqLst
	currPPQ$ = Get string: ppqCount
	if ppqCount == 1
		outPPQ = Read from file: "'input_directory$'/'currPPQ$'"
		Rename: "ppq"
	else 
		currPPQ = Read from file: "'input_directory$'/'currPPQ$'"
		select outPPQ
		plus currPPQ
		Concatenate
		Rename: "ppq"
		removeObject: currPPQ
	endif
	ppqCount = ppqCount + 1
endwhile

removeObject: ppqLst


# END first addition of Fredrik Karlsson 2021-09-10 from original script


#mpt = maximum_phonation_time. # Now computed from samples
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

if apply_calibration == 0 
	To Intensity... 60 0.0 yes 
else 
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
durationVowel = Get total duration
durationStart = durationVowel - 3 

if durationVowel > 3
	Extract part... durationStart durationVowel rectangular 1 no
	Rename... ppq2
elsif durationVowel <=3
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

Text... 0 left 0.5 half Script: Youri Maryn, PhD
12
Select inner viewport... 0.5 7.5 0 0.5 
Axes... 0 1 0 1 
Text... 0 left 0.5 half ##DYSPHONIA SEVERITY INDEX (DSI) IN PRAAT, v.02.01# 

Font size... 8 

Select inner viewport... 0.5 7.5 0 0.5 
Axes... 0 1 0 3
Text... 1 right 2.3 half %%'Name_of_patient$'%
Text... 1 right 1.5 half %%°'date_of_birth$'%
Text... 1 right 0.7 half %%'assessment_date$'%
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

# BEGIN second addition of Fredrik Karlsson 2021-09-10 from original script

#Now store the results
if generate_PDF_files == 1
	Save as PDF file: "'output_directory$'/'speaker_ID$'_'assessment_date$'.pdf"
endif
outTab = Create Table with column names: "outTab", 1, "ID"
Set string value: 1, "ID", speaker_ID$
Append column: "Maximum phonation time"
Set numeric value: 1, "Maximum phonation time", 'mpt:2'
Append column: "Softest intensity of voiced speech"
Set numeric value: 1, "Softest intensity of voiced speech", 'minimumIntensity:2'
Append column: "Maximum fundamental frequency"
Set numeric value: 1, "Maximum fundamental frequency", 'maximumF0:2'
Append column: "Jitter ppq5"
Set numeric value: 1, "Jitter ppq5", 'jitterPpq5:2'

Save as comma-separated file: output_file$

# END second addition of Fredrik Karlsson 2021-09-10 from original script


select all
minus Sound im
minus Sound fh
minus Sound ppq
#Remove
