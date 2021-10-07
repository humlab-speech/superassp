###############TREMOR#################
# coninout.praat is a Praat[6.1.47] script (http://www.praat.org/) 
# that serves as a procedure within tremor.praat.
###############TREMOR#################
# Author: Markus Brückl (markus.brueckl@tu-berlin.de)
# Copyright 2011-2021 Markus Brückl
# License: GNU GPL v3 (http://www.gnu.org/licenses/gpl.html)
######################################
#
######################################
# In: Sound object (quasi-stationary sustained phonation) in
# Out: tremor measurements as text in Info window
######################################

procedure cinout .startTime, .endTime, .selectionOffset, .selectionLength, .windowType$, .windowWidth



conList = Create Strings as file list: "conList", conPath$+"*.wav"
fileN = Get number of strings
for ifile to fileN
   selectObject: conList
   filename$ = Get string: ifile
   name$ = filename$ - ".wav"
   sndID = Read from file: conPath$ + filename$

	soundEnd = Get end time

	startAt  = .startTime + .selectionOffset

	endAt = .endTime
	selEnd = startAt + .selectionLength
	if selEnd < .endTime 
		endAt = selEnd
	endif



	select sndID
	sound = Extract part: startAt, endAt, .windowType$, .windowWidth, 0
	Rename: "aaa"
	name$ = "aaa"


   slength = Get total duration

   call ftrem
   call atrem

if fileReadable path_and_name_of_result_csv$
   appendFileLine: path_and_name_of_result_csv$, 
..."'name$''tab$'
...'ftrm:3''tab$'
...'ftrc:3''tab$'
...'fmodN''tab$'
...'ftrf:3''tab$'
...'ftri:3''tab$'
...'ftrp:3''tab$'
...'ftrcip:3''tab$'
...'ftrps:3''tab$'
...'ftrHNR:2''tab$'
...'atrm:3''tab$'
...'atrc:3''tab$'
...'amodN''tab$'
...'atrf:3''tab$'
...'atri:3''tab$'
...'atrp:3''tab$'
...'atrcip:3''tab$'
...'atrps:3''tab$'
...'atrHNR:2'"

else
   writeFileLine: path_and_name_of_result_csv$, 
..."soundname'tab$'
...FCoM'tab$'
...FTrC'tab$'
...FMon'tab$'
...FTrF [Hz]'tab$'
...FTrI [%]'tab$'
...FTrP'tab$'
...FTrCIP'tab$'
...FTrPS'tab$'
...FCoHNR[dB]'tab$'
...ACoM'tab$'
...ATrC'tab$'
...AMoN'tab$'
...ATrF [Hz]'tab$'
...ATrI [%]'tab$'
...ATrP'tab$'
...ATrCIP'tab$'
...ATrPS'tab$'
...ACoHNR[dB]"

   appendFileLine: path_and_name_of_result_csv$, 
..."'name$''tab$'
...'ftrm:3''tab$'
...'ftrc:3''tab$'
...'fmodN''tab$'
...'ftrf:3''tab$'
...'ftri:3''tab$'
...'ftrp:3''tab$'
...'ftrcip:3''tab$'
...'ftrps:3''tab$'
...'ftrHNR:2''tab$'
...'atrm:3''tab$'
...'atrc:3''tab$'
...'amodN''tab$'
...'atrf:3''tab$'
...'atri:3''tab$'
...'atrp:3''tab$'
...'atrcip:3''tab$'
...'atrps:3''tab$'
...'atrHNR:2'"
endif

   removeObject: sndID
endfor
removeObject: conList
endproc