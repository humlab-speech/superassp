###############TREMOR#################
# coninout.praat is a Praat[6.1.47] script (http://www.praat.org/) 
# that serves as a procedure within tremor.praat.
###############TREMOR#################
# Author: Markus Brückl (markus.brueckl@tu-berlin.de)
# Copyright 2011-2021 Markus Brückl
# License: GNU GPL v3 (http://www.gnu.org/licenses/gpl.html)
# Modified for single file use by Fredrik Karlsson 2021-10
######################################
#
######################################
# In: Sound object (quasi-stationary sustained phonation) in
# Out: tremor measurements as text in Info window
######################################

procedure cinout .startTime, .endTime, .selectionOffset, .selectionLength, .windowType$, .windowWidth



	outTab = Create Table with column names: "outTab", 1, {"Start Time","End Time", "Selection Start","Selection End","FCoM","FTrC","FMon","FTrF [Hz]","FTrI [%]","FTrP","FTrCIP","FTrPS","FCoHNR[dB]","ACoM","ATrC","AMoN","ATrF [Hz]","ATrI [%]","ATrP","ATrCIP","ATrPS","ACoHNR[dB]" }

	sndID = Read from file: "'path_of_sound_to_be_analyzed$'"

	soundEnd = Get end time

	startAt  = .startTime + .selectionOffset


	
	if .endTime == 0
		.endTime = soundEnd
	endif

 	endAt = .endTime

	selEnd = startAt + .selectionLength

	if .selectionLength == 0
		selEnd = .endTime
	endif 

	if selEnd < .endTime 
		endAt = selEnd
	endif

	select sndID
	sndID = Extract part: startAt, endAt, .windowType$, .windowWidth, 0
	Rename: "aaa"
	name$ = "aaa"

	slength = Get total duration

	call ftrem
	call atrem


	select outTab 
	Set numeric value: 1, "Start Time", '.startTime'
	Set numeric value: 1, "End Time", '.endTime'
	Set numeric value: 1, "Selection Start", 'startAt'
	Set numeric value: 1, "Selection End", 'endAt'

	if ftrm <> undefined
		Set numeric value: 1, "FCoM", 'ftrm:3'
	endif
	if ftrc <> undefined
		Set numeric value: 1, "FTrC", 'ftrc:3'
	endif
	if fmodN <> undefined
		Set numeric value: 1, "FMon", 'fmodN:3'
	endif
	if ftrf <> undefined
		Set numeric value: 1, "FTrF [Hz]", 'ftrf:3'
	endif
	if ftri <> undefined
		Set numeric value: 1, "FTrI [%]", 'ftri:3'
	endif
	if ftrp <> undefined
		Set numeric value: 1, "FTrP", 'ftrp:3'
	endif
	if ftrcip <> undefined
		Set numeric value: 1, "FTrCIP", 'ftrcip:3'
	endif
	if ftrps <> undefined
		Set numeric value: 1, "FTrPS", 'ftrps:3'
	endif
	if ftrHNR <> undefined
		Set numeric value: 1, "FCoHNR[dB]", 'ftrHNR:3'
	endif
	if atrm <> undefined
		Set numeric value: 1, "ACoM", 'atrm:3'
	endif
	if atrc <> undefined
		Set numeric value: 1, "ATrC", 'atrc:3'
	endif
	if amodN <> undefined
		Set numeric value: 1, "AMoN", 'amodN:3'
	endif
	if atrf <> undefined
		Set numeric value: 1, "ATrF [Hz]", 'atrf:3'
	endif
	if atri <> undefined
		Set numeric value: 1, "ATrI [%]", 'atri:3'
	endif
	if atrp <> undefined
		Set numeric value: 1, "ATrP", 'atrp:3'
	endif
	if atrcip <> undefined
		Set numeric value: 1, "ATrCIP", 'atrcip:3'
	endif
	if atrps <> undefined
		Set numeric value: 1, "ATrPS", 'atrps:3'
	endif
	if atrHNR <> undefined
		Set numeric value: 1, "ACoHNR[dB]", 'atrHNR:3'
	endif

	Save as comma-separated file: path_and_name_of_result_csv$
	removeObject: sndID
endproc