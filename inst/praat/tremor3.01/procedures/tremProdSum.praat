###############TREMOR#################
# tremProdSum.praat is a Praat[6.0.14] script (http://www.praat.org/) 
# that serves as a procedure within tremor.praat.
###############TREMOR#################
# Author: Markus Brückl (markus.brueckl@tu-berlin.de)
# Copyright 2011-2017 Markus Brückl
# License: GNU GPL v3 (http://www.gnu.org/licenses/gpl.html)
######################################

######################################
# caculate tremor (cyclicality and intensity) product sums
######################################

procedure tremProdSum
show = 1
tris = 0
rank = 0
for ifreq from 1 to sN
   select Matrix frequencies
   trFreq = Get value in cell: ifreq, 1

   if trFreq > 0
      rank += 1
      select Matrix strengths
      trStren = Get value in cell: ifreq, 1
      Create simple Matrix: "xy", 1, 1, "'trFreq'"
      To Pitch
      Scale times to: 0, slength

      Rename... trem_norm
      call tremIntIndex

      tris += trStren * tri

      select Matrix xy
      Remove
   endif

endfor
endproc