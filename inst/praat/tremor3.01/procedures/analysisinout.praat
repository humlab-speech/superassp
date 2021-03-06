###############TREMOR#################
# analysisinout.praat is a Praat[6.0.14] script (http://www.praat.org/) 
# that serves as a procedure within tremor.praat.
###############TREMOR#################
# Author: Markus Brückl (markus.brueckl@tu-berlin.de)
# Copyright 2011-2017 Markus Brückl
# License: GNU GPL v3 (http://www.gnu.org/licenses/gpl.html)
######################################

######################################
# Sound (.wav) in, results (.txt) out
######################################
procedure anainout
pause Record/open and select the sound to be analyzed (in 'Praat Objects')!

info$ = Info
name$ = extractWord$(info$, "Object name: ")

slength = Get total duration

call ftrem
call atrem

echo 
...Soundname: 'name$''newline$'
...'newline$'
...frequency contour magnitude (FCoM): 'ftrm:3''newline$'
...amplitude contour magnitude (ACoM): 'atrm:3''newline$'
...'newline$'
...frequency tremor cyclicality (FTrC): 'ftrc:3''newline$'
...amplitude tremor cyclicality (ATrC): 'atrc:3''newline$'
...'newline$'
...frequency tremor frequency (FTrF): 'ftrf:3' Hz'newline$'
...amplitude tremor frequency (ATrF): 'atrf:3' Hz'newline$'
...'newline$'
...frequency tremor intensity index (FTrI): 'ftri:3' %'newline$'
...amplitude tremor intensity index (ATrI): 'atri:3' %'newline$'
...'newline$'
...frequency tremor power index (FTrP): 'ftrp:3''newline$'
...amplitude tremor power index (ATrP): 'atrp:3''newline$'
...'newline$'
...frequency tremor cyclicality intensity product (FTrCIP): 'ftrcip:3''newline$'
...amplitude tremor cyclicality intensity product (ATrCIP): 'atrcip:3''newline$'
...'newline$'
...frequency tremor product sum (FTrPS): 'ftrps:3''newline$'
...amplitude tremor product sum (ATrPS): 'atrps:3'

endproc