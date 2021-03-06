###############TREMOR#################
# getCyclicality.praat is a Praat[6.0.14] script (http://www.praat.org/) 
# that serves as a procedure within tremor.praat.
###############TREMOR#################
# Author: Markus Brückl (markus.brueckl@tu-berlin.de)
# Copyright 2011-2017 Markus Brückl
# License: GNU GPL v3 (http://www.gnu.org/licenses/gpl.html)
######################################

######################################
# Read from 1 FRAME (!) Praat Pitch ojects
######################################

procedure cycli

Save as text file: "./temp"
Read Strings from raw text file: "./temp"

stringN = Get number of strings
sN = 0
freq = 100
for istring from 10 to stringN
   select Strings temp
   tEkst$ = Get string: istring
   tEkst$ = replace$(tEkst$, " ", "",100)
#echo 'istring' 'tEkst$'
#pause
   if startsWith(tEkst$, "maxnCandidates")
      cN = extractNumber(tEkst$, "maxnCandidates=")
      Create simple Matrix: "strengths", 'cN', 1, "0"
      Create simple Matrix: "frequencies", 'cN', 1, "0"
   elsif startsWith(tEkst$, "intensity")
      trm = extractNumber(tEkst$, "intensity=")
   elsif startsWith(tEkst$, "frequency")
      freq = extractNumber(tEkst$, "frequency=")
   elsif startsWith(tEkst$, "strength") and (freq <= maxTr)
      sN +=1
      strength = extractNumber(tEkst$, "strength=")
      select Matrix strengths
      Set value: 'sN', 1, 'strength'
      select Matrix frequencies
      Set value: 'sN', 1, 'freq'
#echo 'sN' 'streng_'sN''
#pause
      endif
   endif
endfor
select Matrix strengths
trc = Get maximum
#pause
select Strings temp
#plus Matrix temp
Remove
filedelete ./temp

endproc