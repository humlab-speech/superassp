###############TREMOR#################
# analysisinout.praat is a Praat[6.1.47] script (http://www.praat.org/) 
# that serves as a procedure within tremor.praat.
###############TREMOR#################
# Author: Markus Brückl (markus.brueckl@tu-berlin.de)
# Copyright 2011-2021 Markus Brückl
# License: GNU GPL v3 (http://www.gnu.org/licenses/gpl.html)
######################################

######################################
# Sound (.wav) in, results (.txt) out
######################################
procedure anainout
pause Record/open and select the sound to be analyzed (in 'Praat Objects')!

info$ = Info
sndID = extractNumber(info$, "Object id: ")
name$ = extractWord$(info$, "Object name: ")


slength = Get total duration

call ftrem
call atrem

writeInfoLine: "Soundname:", tab$, name$, newline$
appendInfoLine: "FCoM:'tab$''ftrm:3'", tab$, "frequency contour magnitude"
appendInfoLine: "ACoM:'tab$''atrm:3'", tab$, "amplitude contour magnitude", newline$
appendInfoLine: "FTrC:'tab$''ftrc:3'", tab$, "(maximum) frequency tremor cyclicality"
appendInfoLine: "ATrC:'tab$''atrc:3'", tab$, "(maximum) amplitude tremor cyclicality", newline$
appendInfoLine: "FMoN:'tab$''fmodN'", tab$, "number of frequency modulations above thresholds"
appendInfoLine: "AMoN:'tab$''amodN'", tab$, "number of amplitude modulations above thresholds", newline$
appendInfoLine: "FTrF:'tab$''ftrf:3'Hz", tab$, "(strongest) frequency tremor frequency"
appendInfoLine: "ATrF:'tab$''atrf:3'Hz", tab$, "(strongest) amplitude tremor frequency", newline$
appendInfoLine: "FTrI:'tab$''ftri:3'%", tab$, "frequency tremor intensity index at FTrF"
appendInfoLine: "ATrI:'tab$''atri:3'%", tab$, "frequency tremor intensity index at ATrF", newline$
appendInfoLine: "FTrP:'tab$''ftrp:3'", tab$, "frequency tremor power index at FTrF"
appendInfoLine: "ATrP:'tab$''atrp:3'", tab$, "frequency tremor power index at ATrF", newline$
appendInfoLine: "FTrCIP:'tab$''ftrcip:3'%", tab$, "frequency tremor cyclicality intensity product at FTrF"
appendInfoLine: "ATrCIP:'tab$''atrcip:3'%", tab$, "amplitude tremor cyclicality intensity product at ATrF", newline$
appendInfoLine: "FTrPS:'tab$''ftrps:3'", tab$, "frequency tremor product sum"
appendInfoLine: "ATrPS:'tab$''atrps:3'", tab$, "amplitude tremor product sum", newline$
appendInfoLine: "FCoHNR:'tab$''ftrHNR:2'dB", tab$, "frequency contour harmonicity-to-noise ratio"
appendInfoLine: "ACoHNR:'tab$''atrHNR:2'dB", tab$, "amplitude contour harmonicity-to-noise ratio"


endproc