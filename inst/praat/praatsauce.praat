###############
# PraatSauce
###############

# Copyright (c) 2021-2023 Fredrik Karlsson
# Based on shellSause.praat in praatsauce 
# which is Copyright (c) 2018-2019 James Kirby

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


# Portions of PraatSauce are based on

# spectralTiltMaster.praat
# version 0.0.5
# copyright 2009-2010 Timothy Mills
# <mills.timothy@gmail.com>

# VoiceSauce 
# version 1.31
# http://www.seas.ucla.edu/spapl/voicesauce/


form File and measures
    sentence inputfile /Users/frkkan96/Desktop/a2.wav
    real beginTime 0
    real endTime 0
    integer channel 1
	 comment 1= n equidistant points
    comment 2=   option every n milliseconds
    natural measure 2
	 natural Points 5
    boolean resample_to_16k 1    
    boolean pitchTracking 1
    boolean formantMeasures 1
    boolean spectralMeasures 1
    positive windowLength 0.025
    positive windowPosition 0.5
    positive maxFormantHz 5000
    positive spectrogramWindow 0.005
    boolean useExistingPitch 0 
    positive f0min 50
    positive f0max 300
    real timeStep 0
    integer maxNumFormants 5
    positive preEmphFrom 50
    boolean formantTracking 1
    positive F1ref 500
    positive F2ref 1500
    positive F3ref 2500
    boolean useExistingFormants 0
    boolean useBandwidthFormula 0
    sentence outputfile /Users/frkkan96/Desktop/spectral_measures.txt
endform

###
### Make sure that spectral measures can be calculated if selected, override user settings
###

#if spectralMeasures
#	pitchTracking = 1
#	formantMeasures = 1
#endif

spectralMeasures = 1
pitchTracking = 1
formantMeasures = 1

###
## Build up the output table based on the user's choices

outTab = Create Table with column names: "outTable", 0, "t"
#Append column: "t_ms"

## Add header columns for selected measures
if pitchTracking 
	select outTab
	Append column: "f0"
endif

if formantMeasures
	select outTab
	Append column: "F1"
	Append column: "F2"
	Append column: "F3"
	Append column: "B1"
	Append column: "B2"
	Append column: "B3"
endif

if spectralMeasures

    specm$[1] = "H1u"
    specm$[2] = "H2u"
    specm$[3] = "H4u"
    specm$[4] = "H2Ku"
    specm$[5] = "H5Ku"
    specm$[6] = "A1u"
    specm$[7] = "A2u"
    specm$[8] = "A3u"
    specm$[9] = "H1H2u"
    specm$[10] = "H2H4u"
    specm$[11] = "H1A1u"
    specm$[12] = "H1A2u"
    specm$[13] = "H1A3u"
    specm$[14] = "H2KH5Ku"
    specm$[15] = "H1c"
    specm$[16] = "H2c"
    specm$[17] = "H4c"
    specm$[18] = "A1c"
    specm$[19] = "A2c"
    specm$[20] = "A3c"
    specm$[21] = "H1H2c"
    specm$[22] = "H2H4c"
    specm$[23] = "H1A1c"
    specm$[24] = "H1A2c"
    specm$[25] = "H1A3c"
    specm$[26] = "CPP"
    specm$[27] = "HNR05"
    specm$[28] = "HNR15"
    specm$[29] = "HNR25"
    specm$[30] = "HNR35"
    for i from 1 to 30
        select outTab
        Append column: specm$[i]
    endfor
endif




## Load Sound
Read from file: inputfile$
## Note that Sound is not resampled b/c later we use To Formant (burg)... 
## which resamples to twice the frequency of maxFormant (so 10k-11k usually)
soundID = selected("Sound")

## If selected, downsample
## VoiceSauce only really does this b/c it is faster for STRAIGHT
## Given how the Burg algorithm works, input will be resampled 
## at the formant estimation stage no matter what.
if resample_to_16k == 1
    Resample... 16000 50
    resampledID = selected("Sound")
    select 'soundID'
    Remove
    select 'resampledID'
    Rename... 'basename$'
    soundID = selected("Sound")
endif  

## If selected, extract the channel of interest
## if e.g. you have audio on channel 1 and EGG on channel 2
if channel
    Extract one channel... channel
    monoID = selected("Sound")
    select 'soundID'
    Remove
    select 'monoID'
    Rename... 'basename$'
    soundID = selected("Sound")
endif

## Code from praat_formant_burg, as an alterantive way of handling analysis of a portion of a file
#selectObject: soundID
#dur = Get total duration

## Check that start or end times should be condidered, and that they are within 
## ok limits.
#if  ( beginTime > 0.0 or endTime > 0.0 ) and (beginTime >= 0.0 and endTime <= dur)
#	
#	selectObject: sound
#	# Preserve times so that start end end record may be computed later
#	soundPart = Extract part: beginTime, endTime, windowShape$, relativeWidth, 1
#	selectObject: sound
#	Remove
#	soundID = soundPart
#endif 

# The praatsauce way
interval_start = max(0, 'beginTime')
    
select soundID
soundEnd = Get end time
soundStart = Get start time

interval_end = endTime
if endTime == 0	
	interval_end = soundEnd
endif

## Determine how many timepoints we're measuring at
if measure = 1
   timepoints = points
elsif measure = 2
   timepoints = round(((interval_end - interval_start)*1000)/points)
endif
        
#Not really needed for an unsupervised script like this 
        manualCheck = 0



############################################
## Load/create Pitch, Formant, etc. objects
############################################

###
## Not sure of the best way to do this. This way seems to be the cleanest,
## because the objects are alway sonly loaded or created once. 
##
## However, it is possible that they are created redundantly, because if
## a given file doesn't have any intervals of interest, nothing will
## be measured. 
###

if pitchTracking
		# if you want to load an existing object from disk...
	if useExistingPitch
    	if fileReadable ("'inputdir$''basename$'.Pitch")
        	Read from file... 'inputdir$''basename$'.Pitch
        else
        	exit Cannot load Pitch object <'basename$'.Pitch>.
        endif
    # else create
    else
    	select 'soundID'
        #To Pitch... 0 'f0min' 'f0max'
        ## TODO April 2019: add this as a user option
        noprogress To Pitch (ac)... 0 'f0min' 15 0 0.03 0.45 0.01 0.35 0.14 'f0max'
		## This will result in two Pitch objects, but that's OK (I hope)...
 		# Interpolate
		## This is maybe nice for some applications but hallucinates f0 in clearly voiceless regions!!
    endif
	## ... since only the second one will get referred to from now on
    pitchID = selected("Pitch")
endif

if formantMeasures
	# if you want to load an existing object from disk...
    if useExistingFormants
    	# Load existing Formant object if available and selected
        if fileReadable ("'inputdir$''basename$'.Formant")
        	Read from file... 'inputdir$''basename$'.Formant
        else
        	exit Cannot load Formant object <'basename$'.Formant>.
        endif
	## else create 
    else
    	select 'soundID'
        noprogress To Formant (burg)... timeStep maxNumFormants maxFormantHz windowLength preEmphFrom
	endif

	formantID = selected("Formant")

	# moved from formantMeasures.praat 2020-08-18
    if formantTracking = 1
        # Tracking cleans up the tracks a little.  The original Formant object is then discarded.
        ## mar 19: should really make the number of tracks a user parameter
        ## also need to tune it possibly for each frame, because if the Formant
        ## object has fewer values than the numTracks parameter, the command
        ## will fail.
        minFormants = Get minimum number of formants
        if 'minFormants' = 2
            noprogress Track... 2 f1ref f2ref f3ref 3850 4950 1 1 1
        else
            noprogress Track... 3 f1ref f2ref f3ref 3850 4950 1 1 1
        endif
        trackedFormantID = selected("Formant")
        select 'formantID'
        Remove
        formantID = trackedFormantID
    endif

endif

if spectralMeasures
    # TODO: add option to read from disk as above
    ### Create Harmonicity objects ###

    ## here we use a 1 period window; the Praat 
    ## default of 4.5 periods per window produces
    ## much less accurate estimates
    select 'soundID'
    Filter (pass Hann band): 0, 500, 100
    Rename... 'basename$'_500
    noprogress To Harmonicity (cc): 0.01, f0min, 0.1, 1.0
    hnr05ID = selected ("Harmonicity")
    select 'soundID'
    noprogress Filter (pass Hann band): 0, 1500, 100
    Rename... 'basename$'_1500
    noprogress To Harmonicity (cc): 0.01, f0min, 0.1, 1.0
    hnr15ID = selected ("Harmonicity")
    select 'soundID'
    noprogress Filter (pass Hann band): 0, 2500, 100
    Rename... 'basename$'_2500
    noprogress To Harmonicity (cc): 0.01, f0min, 0.1, 1.0
    hnr25ID = selected ("Harmonicity")
    select 'soundID'
    noprogress Filter (pass Hann band): 0, 3500, 100
    Rename... 'basename$'_3500
    noprogress To Harmonicity (cc): 0.01, f0min, 0.1, 1.0
    hnr35ID = selected ("Harmonicity")
    ### (end create Harmonicity objects ###
endif

			
################################
### Pitch tracking
################################

if pitchTracking
    select pitchID
    plus soundID

    execute pitchTracking.praat 'interval_start' 'interval_end' 'windowPosition' 'windowLength' 'manualCheck' 1 'measure' 'timepoints' 'points'
   
    ### Save output Matrix
    select Matrix pitchaverages
    pitchResultsID = selected("Matrix")
endif
### (end of pitch tracking)

###################### 
### Formant measures
###################### 

if formantMeasures
    select 'soundID'
    
    plus 'formantID'
    execute formantMeasures.praat 'interval_start' 'interval_end' 'windowPosition' 'windowLength' 'useBandwidthFormula' 'useExistingFormants' 'inputdir$' 'basename$' 'timeStep' 'maxNumFormants' 'maxFormantHz' 'preEmphFrom' 'f1ref' 'f2ref' 'f3ref' 'spectrogramWindow' 'measure' 'timepoints' 'points' 1 20 1 
    formantID = selected("Formant")
    select Matrix FormantAverages
    formantResultsID = selected("Matrix")
endif
### (end of formant measures)

###################################################################################### 
### Spectral corrections (including H1*, H2*, H4, A1*, A2*, A3* from Iseli et al.)
###################################################################################### 

if spectralMeasures
    select 'soundID'
    plus 'formantID'
    plus 'pitchID'
    plus 'hnr05ID'
    plus 'hnr15ID'
    plus 'hnr25ID'
    plus 'hnr35ID'
    execute spectralMeasures.praat 'interval_start' 'interval_end' 'windowPosition' 'windowLength' 'useBandwidthFormula' 'inputdir$' 'maxDisplayHz' 'measure' 'timepoints' 'points' 'f0min' 'f0max'

    ## Assign ID to output matrix
    select Matrix IseliMeasures
    iseliResultsID = selected("Matrix")

endif 
### (end of spectralMagnitude measure)

## Store all results
    select outTab


    

for t from 1 to timepoints

	Append row
   #As we are appending, the current row number is also an index to the last row
   row = Get number of rows
   #Into this row, we start copying data from the output matricies

   if pitchTracking
		select 'pitchResultsID'
     mspoint = Get value in cell... t 1
     currentPitch = Get value in cell... t 2
		select outTab
     Set numeric value... 'row' t 'mspoint:6'
		Set numeric value... 'row' f0 'currentPitch'
	endif
        
   if formantMeasures 
   	select 'formantResultsID'
     	mspoint = Get value in cell... t 1
		select outTab
     Set numeric value... 'row' t 'mspoint:6'
            
		for formant from 2 to 4
			select 'formantResultsID'
			currentFormant = Get value in cell... t 'formant'
			fI = formant - 1
			fn$ = "F'fI'" 
			select outTab
			Set numeric value... 'row' 'fn$' 'currentFormant:3'
		endfor
            
		for bandwidth from 5 to 7
			select 'formantResultsID'
			currentBandwidth = Get value in cell... t 'bandwidth'
			formant = bandwidth - 4
			bN$ = "B'formant'"
			select outTab
			Set numeric value... 'row' 'bN$' 'currentBandwidth:3'
		endfor
	endif

	if spectralMeasures
		select 'iseliResultsID'
     mspoint = Get value in cell... t 1
		select outTab
     Set numeric value... 'row' t 'mspoint:6'
		for measurement from 2 to 31
			select 'iseliResultsID'
			aMeasure = Get value in cell... t 'measurement'
			mindex = measurement -1
			nname$ = specm$['mindex']
			select outTab 
			Set numeric value... 'row' 'nname$' 'aMeasure:3'
		endfor
	endif

endfor

## clean up
select all
minus outTab

nocheck minus 'pitchID'
nocheck minus 'formantID'
nocheck minus 'hnr05ID'
nocheck minus 'hnr15ID'
nocheck minus 'hnr25ID'
nocheck minus 'hnr35ID'
Remove
select outTab
Save as comma-separated file: outputfile$

