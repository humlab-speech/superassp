# formantMeasures.praat
# version 0.2.2
# James Kirby <j.kirby@ed.ac.uk>

# based in part on

# spectralTiltMaster.praat
# version 0.0.5
# copyright 2009-2010 Timothy Mills

# These scripts are released under the GNU General Public License version 3.0 
# The included file "gpl-3.0.txt" or the URL "http://www.gnu.org/licenses/gpl.html" 
# contains the full text of the license.

# It is designed to work as part of the "praatsauce" script suite,
# which can be obtained from:
#
#      https://github.com/kirbyj/praatsauce
#
# This script is released under the GNU General Public License version 3.0 
# The included file "gpl-3.0.txt" or the URL "http://www.gnu.org/licenses/gpl.html" 
# contains the full text of the license.
#
# This script takes a Sound object and an associated TextGrid, and
# determines the first three formant frequencies.  It can work in a 
# fully-automated fashion, but users are advised that formant 
# tracking errors can be common.  Manual checking of the output is
# recommended.

include getbw_HawksMiller.praat

form Parameters for formant measurement
    comment Seclection start and end times
    positive startTime 0
	  positive endTime 0
    comment Window parameters
    positive windowPosition 0.5
    positive windowLength 0.025
    boolean useBandwidthFormula 1
    boolean useExistingFormants 0
    text inputdir /home/username/data/ 
    text basename myFile
    comment Leave timeStep at 0 for auto.
    real timeStep 0
    integer maxNumFormants 5
    positive maxFormantHz 5500
    positive preEmphFrom 50
    positive f1ref 500
    positive f2ref 1500
    positive f3ref 2500
    positive spectrogramWindow 0.005
    positive measure 2
    positive timepoints 3
    positive timestep 1
    boolean formantTracking 1
    positive smoothWindowSize 20
    positive smoother 1
endform

###
### First, check that proper objects are present and selected.
###
numSelectedSound = numberOfSelected("Sound")

numSelectedFormant = numberOfSelected("Formant")
if (numSelectedSound<>1  or numSelectedFormant <>1)
 exit Select only one Sound object, one TextGrid object and one Formant object.
endif
name$ = selected$("Sound")
soundID = selected("Sound")
formantID = selected("Formant")
### (end object check)



###
### Third, decide what times to measure at.
###
d = startTime + (timestep/1000)
## If equidistant points: compute based on number of points
if measure = 1
    diff = (endTime - startTime) / (timepoints+1)
## If absolute: take a measurement every timepoints/1000 points
elsif measure = 2
    diff = timestep / 1000
endif
for point from 1 to timepoints
    mid'point' = d
    d = d + diff
endfor
### (end time point selection)

###
### Fourth, build Matrix object to hold results
### column 1 holds time of measurement
### (relative to distance from startTime)
### columns 2-4 hold Fe, F2, F3
###
Create simple Matrix... FormantAverages timepoints 7 0
matrixID = selected("Matrix")
### (end of build Matrix object)


###
### Store measurements for each timepoint
for i from 1 to timepoints

    select 'formantID'
    f1 = Get value at time... 1 mid'i' Hertz Linear
    f2 = Get value at time... 2 mid'i' Hertz Linear
    f3 = Get value at time... 3 mid'i' Hertz Linear
 
    # bandwidths
    if useBandwidthFormula = 1
        selectObject: "Pitch 'basename$'"
        n_f0md = Get value at time... mid'i' Hertz Linear
        select 'formantID'
        @getbw_HawksMiller(n_f0md, f1) 
        bw1 = getbw_HawksMiller.result
        @getbw_HawksMiller(n_f0md, f2)
        bw2 = getbw_HawksMiller.result
        @getbw_HawksMiller(n_f0md, f3)
        bw3 = getbw_HawksMiller.result
    else
        bw1 = Get bandwidth at time... 1 mid'i' Hertz Linear
        bw2 = Get bandwidth at time... 2 mid'i' Hertz Linear
        bw3 = Get bandwidth at time... 3 mid'i' Hertz Linear
    endif
 
    # can't have undefineds in your Praat Matrices, sorry              
    if f1 = undefined
       f1 = 0
    endif
    if f2 = undefined
       f2 = 0
    endif
    if f3 = undefined
       f3 = 0
    endif
    if bw1 = undefined
       bw1 = 0
    endif
    if bw2 = undefined
       bw2 = 0
    endif
    if bw3 = undefined
       bw3 = 0
    endif
 
    ###
    ### Write to Matrix
    ### 
    select 'matrixID'
    
    # set first col to ms time
    Set value... i 1 mid'i'
 
    # record formants and bandwidths
    Set value... i 2 'f1'
    Set value... i 3 'f2'
    Set value... i 4 'f3'
    Set value... i 5 'bw1'
    Set value... i 6 'bw2'
    Set value... i 7 'bw3'
### End the 'for' loop over timepoints
endfor


select 'formantID'
