###############TREMOR#################
# signleruninout.praat is a Praat[6.0.14] script (http://www.praat.org/) 
# that serves as a procedure within console_tremor.praat.
###############TREMOR#################
# This is a version of the original 'runinout.praat' script 
# which was created by Markus Brückl (markus.brueckl@tu-berlin.de)
# but modified by Fredrik Karlsson (fredrik.k.karlsson@umu.se) to 
# allow professing of a single file for use in the Visible Speech (VISP) platform (visp.humlab.umu.se)
# 
# License: GNU GPL v3 (http://www.gnu.org/licenses/gpl.html)
######################################

######################################
# Sound (.wav) in, results (.csv) out
######################################

procedure singleruninout: .path_of_sound_to_be_analyzed$, .path_and_name_of_result_csv$, .startTime, .endTime, .selectionOffset, .selectionLength, .windowType$, .windowWidth


outTab = Create Table with column names: "outTab", 1, { "FCoM", "FTrC", "FTrF", "FTrI", "FTrP", "FTrCIP", "FTrPS", "ACoM", "ATrC", "ATrF", "ATrI", "ATrP", "ATrCIP", "ATrPS" }

sound = Read from file: "'.path_of_sound_to_be_analyzed$'"

soundEnd = Get end time

startAt  = .startTime + .selectionOffset

endAt = .endTime
selEnd = startAt + .selectionLength
if selEnd < .endTime 
	endAt = selEnd
endif



select sound
sound = Extract part: startAt, endAt, .windowType$, .windowWidth, 0
Rename: "aaa"
name$ = "aaa"

slength = Get total duration

call ftrem
call atrem

select outTab 
if ftrm <> undefined 
	Set numeric value: 1, "FCoM", 'ftrm:3'
endif
if ftrc <> undefined
	Set numeric value: 1, "FTrC", 'ftrc:3'
endif
if ftrf <> undefined
	Set numeric value: 1, "FTrF", 'ftrf:3'
endif
if ftri <> undefined
	Set numeric value: 1, "FTrI", 'ftri:3'
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
if atrm <> undefined
	Set numeric value: 1, "ACoM", 'atrm:3'
endif
if atrc <> undefined
	Set numeric value: 1, "ATrC", 'atrc:3'
endif
if atrf <> undefined
	Set numeric value: 1, "ATrF", 'atrf:3'
endif
if atri <> undefined
	Set numeric value: 1, "ATrI", 'atri:3'
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

Save as comma-separated file: .path_and_name_of_result_csv$

select sound
Remove

endproc