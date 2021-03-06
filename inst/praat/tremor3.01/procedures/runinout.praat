###############TREMOR#################
# runinout.praat is a Praat[6.0.14] script (http://www.praat.org/) 
# that serves as a procedure within tremor.praat.
###############TREMOR#################
# Author: Markus Brückl (markus.brueckl@tu-berlin.de)
# Copyright 2011-2017 Markus Brückl
# License: GNU GPL v3 (http://www.gnu.org/licenses/gpl.html)
######################################

######################################
# Sounds (.wav) in, results (.csv) out
######################################

procedure rinout
beginPause ("Paths")
   word ("Path of sounds to be analyzed", "../")
   word ("Path and name of result csv", "../restab_tremor")
endPause ("OK", 1)

sourcedirec$ = path_of_sounds_to_be_analyzed$
resultdirec$ = path_and_name_of_result_csv$

filedelete 'resultdirec$'.csv
fileappend "'resultdirec$'.csv" 
...soundname'tab$'
...FCoM'tab$'
...FTrC'tab$'
...FTrF [Hz]'tab$'
...FTrI [%]'tab$'
...FTrP'tab$'
...FTrCIP'tab$'
...FTrPS'tab$'
...ACoM'tab$'
...ATrC'tab$'
...ATrF [Hz]'tab$'
...ATrI [%]'tab$'
...ATrP'tab$'
...ATrCIP'tab$'
...ATrPS
...'newline$'

Create Strings as file list... list 'sourcedirec$'*.wav
numberOfFiles = Get number of strings
for ifile from 1 to numberOfFiles
   select Strings list
   fileName$ = Get string... ifile
   name$ = fileName$ - ".wav"

   Read from file... 'sourcedirec$''name$'.wav

   slength = Get total duration

   call ftrem
   call atrem

   fileappend "'resultdirec$'.csv" 
...'name$''tab$'
...'ftrm:3''tab$'
...'ftrc:3''tab$'
...'ftrf:3''tab$'
...'ftri:3''tab$'
...'ftrp:3''tab$'
...'ftrcip:3''tab$'
...'ftrps:3''tab$'
...'atrm:3''tab$'
...'atrc:3''tab$'
...'atrf:3''tab$'
...'atri:3''tab$'
...'atrp:3''tab$'
...'atrcip:3''tab$'
...'atrps:3'
...'newline$'

   select Sound 'name$'
   Remove
endfor

select Strings list
Remove

endproc