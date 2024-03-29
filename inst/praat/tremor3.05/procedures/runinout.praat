###############TREMOR#################
# runinout.praat is a Praat[6.1.47] script (http://www.praat.org/) 
# that serves as a procedure within tremor.praat.
###############TREMOR#################
# Author: Markus Brückl (markus.brueckl@tu-berlin.de)
# Copyright 2011-2021 Markus Brückl
# License: GNU GPL v3 (http://www.gnu.org/licenses/gpl.html)
######################################

######################################
# Sounds (.wav) in, results (.txt) out
######################################

procedure rinout
beginPause ("Paths")
#   word ("Path of wav-files to be analyzed", "./")
   word ("Path and name of result txt", "./results/restab_tremor.txt")
endPause ("OK", 1)

#sourcedirec$ = path_of_wav-files_to_be_analyzed$
sourcedirec$ = conPath$
resultdirec$ = path_and_name_of_result_txt$

writeFileLine: resultdirec$, 
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

listID = Create Strings as file list: "list", sourcedirec$ + "*.wav"
numberOfFiles = Get number of strings
for ifile from 1 to numberOfFiles
   select Strings list
   fileName$ = Get string: ifile
   name$ = fileName$ - ".wav"

   sndID = Read from file: sourcedirec$ + name$ +".wav"

   slength = Get total duration

   call ftrem
   call atrem

appendFileLine: resultdirec$, 
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

   removeObject: sndID
endfor

removeObject: listID

endproc