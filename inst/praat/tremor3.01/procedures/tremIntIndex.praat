###############TREMOR#################
# tremIntIndex.praat is a Praat[6.0.14] script (http://www.praat.org/) 
# that calculates Tremor Intensity Indices [%] within tremor.praat.
###############TREMOR#################
# Author: Markus Brückl (markus.brueckl@tu-berlin.de)
# Copyright 2011-2017 Markus Brückl
# License: GNU GPL v3 (http://www.gnu.org/licenses/gpl.html)
######################################

procedure tremIntIndex
   select Sound trem
   plus Pitch trem_norm
   To PointProcess (peaks)... yes no
   Rename... Maxima
   numberofMaxPoints = Get number of points
   tri_max = 0
   noFMax = 0
   for iPoint from 1 to numberofMaxPoints
      select PointProcess Maxima
      ti = Get time from index... iPoint
      select Sound trem
      tri_Point = Get value at time... Average ti Sinc70
      if tri_Point = undefined
         tri_Point = 0
         noFMax += 1
      endif
      tri_max += abs(tri_Point)
   endfor

if mode = 1 and show = 1
   select Sound trem
   plus PointProcess Maxima
   Edit
   beginPause: ""
      comment: "Normalized and de-declined 'contType$' contour and"
      comment: "maxima according to its 'rank'. modulation frequency"
   endPause: "Continue", 1
elsif mode = 1 and show = 0
   select Sound trem
   plus PointProcess Maxima
   Edit
   beginPause: ""
      comment: "Normalized and de-declined 'contType$' contour and"
      comment: "maxima according to its strongest modulation frequency"
   endPause: "Continue", 1
endif
   
# tri_max:= (mean) procentual deviation of contour maxima from mean contour at trf
   numberofMaxima = numberofMaxPoints - noFMax
   tri_max = 100 * tri_max/numberofMaxima

   select Sound trem
   plus Pitch trem_norm
   To PointProcess (peaks)... no yes
   Rename... Minima
   numberofMinPoints = Get number of points
   tri_min = 0
   noFMin = 0
   for iPoint from 1 to numberofMinPoints
      select PointProcess Minima
      ti = Get time from index... iPoint
      select Sound trem
      tri_Point = Get value at time... Average ti Sinc70
      if tri_Point = undefined
         tri_Point = 0
         noFMin += 1
      endif
      tri_min += abs(tri_Point)
   endfor

if mode = 1 and show = 1
   select Sound trem
   plus PointProcess Minima
   Edit
   beginPause: ""
      comment: "Normalized and de-declined 'contType$' contour and"
      comment: "minima according to its 'rank'. modulation frequency"
   endPause: "Continue", 1
elsif mode = 1 and show = 0
   select Sound trem
   plus PointProcess Minima
   Edit
   beginPause: ""
      comment: "Normalized and de-declined 'contType$' contour and"
      comment: "minima according to its strongest modulation frequency"
   endPause: "Continue", 1
endif

# tri_min:= (mean) procentual deviation of contour minima from mean contour at trf
   numberofMinima = numberofMinPoints - noFMin
   tri_min = 100 * tri_min/numberofMinima

   tri = (tri_max + tri_min) / 2

   select Pitch trem_norm
   plus PointProcess Maxima
   plus PointProcess Minima
   Remove
   
endproc