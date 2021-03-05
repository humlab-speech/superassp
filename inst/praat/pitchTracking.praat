# pitchTracking.praat
# version 0.0.1


# # Copyright (c) 2021-2023 Fredrik Karlsson
# Based on pitchTracking.praat in praatsauce 
# which is a heavily modified (2011-2017 by James Kirby <j.kirby@ed.ac.uk>)
# version of a script that is copyright 2009-2010 Timothy Mills <mills.timothy@gmail.com>
# 

# This script collects pitch tracks according to the methods use in praatsauce
#
#      https://github.com/kirbyj/praatsauce
#
# but modified in a way that makes it easier to apply it only to a single file and time window.
#
# This script is released under the GNU General Public License version 3.0 
# The included file "gpl-3.0.txt" or the URL "http://www.gnu.org/licenses/gpl.html" 
# contains the full text of the license.

form Parameters for f0 measurement
 comment Start and end times of processing
 real startTime 0
 real endTime 0
 real windowPosition
 positive windowLength
 boolean manualCheck 1
 boolean outputToMatrix 0
 positive measure 2
 positive timepoints 10
 positive timestep 1
endform


###
### Decide what times to measure at.
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
### First, check that proper objects are present and selected.
###
numSelectedSound = numberOfSelected("Sound")
numPitch = numberOfSelected("Pitch")
if (numSelectedSound<>1 or numPitch<>1)
 exit Select only one Sound and one Pitch object.
endif
name$ = selected$("Sound")
soundID = selected("Sound")
pitchID = selected("Pitch")
### (end object check)

### ... not needed code removed from praatsauce

###
### Fourth, build Matrix object to hold results
### each row represents a timepoint
### second column represents absolute timepoint of the measurement (distance from startTime)
### third column represents an f0 measurement
if outputToMatrix
	#writeInfoLine: timepoints
   Create simple Matrix... pitchaverages timepoints 3 0
    matrixID = selected("Matrix")
endif
### (end build Matrix object)

###
### Fifth, store a measurement at each timepoint.
###
for i from 1 to timepoints

	if outputToMatrix
		select 'pitchID'
		f0 = Get value at time... mid'i' Hertz Linear
		timestamp = mid'i'
		if f0 = undefined
			f0 = 0 
		endif
		select 'matrixID'

        # find time of measurement, relative to startTime
        #absPoint = mid'i' - startTime
        #absPoint = round( (mid'i' - startTime)*1000 ) / 1000

        # set first value to ms time
		#Set value... i 1 absPoint
		Set value... i 1 mid'i'
        # set second value to f0
        Set value... i 2 'f0'
	else
        printline "'name$''tab$''f0'"
	pause
	endif
endfor
### (end measurement loop)

select 'soundID'


