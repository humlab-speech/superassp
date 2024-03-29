###############TREMOR#################
# tremProdSum.praat is a Praat[6.1.47] script (http://www.praat.org/) 
# that serves as a procedure within tremor.praat.
###############TREMOR#################
# Author: Markus Brückl (markus.brueckl@tu-berlin.de)
# Copyright 2011-2021 Markus Brückl
# License: GNU GPL v3 (http://www.gnu.org/licenses/gpl.html)
######################################

######################################
# caculate tremor (cyclicality and intensity) product sums
######################################

procedure tremProdSum
show = 1
tris = 0
rank = 0
tri = 0

for ifreq from 1 to sN
   selectObject: torID
   trFreq = Get value: ifreq, 2
   trStren = Get value: ifreq, 3
   if trStren > tremthresh and trm > tremMagThresh
      rank += 1
      
      nCanID = Create simple Matrix: "nextCand", 1, 1, "'trFreq'"
      
      pitrem_normID = To Pitch
      Scale times to: 0, slength
      Rename: "trem_norm_'rank'"

      call tremIntIndex
      tris += trStren * tri

#pause 'trStren' 'tab$' 'tri'

      removeObject: nCanID

      if ifreq = 1 and contType$ = "frequency"
         ftrf = trFreq
         ftri = tri
         ftrp = ftri * ftrf/(ftrf+1)
         ftrcip = ftri * ftrc
      elsif ifreq = 1 and contType$ = "amplitude"
         atrf = trFreq
         atri = tri
         atrp = atri * atrf/(atrf+1)
         atrcip = atri * atrc
      endif
   else
      if ifreq = 1 and contType$ = "frequency" and nan_out = 2
         ftrf = undefined
         ftri = undefined
         ftrp = undefined
         ftrcip = undefined
      elsif ifreq = 1 and contType$ = "amplitude" and nan_out = 2
         atrf = undefined
         atri = undefined
         atrp = undefined
         atrcip = undefined
      elsif ifreq = 1 and contType$ = "frequency" and nan_out = 1
         ftrf = 0
         ftri = 0
         ftrp = 0
         ftrcip = 0
      elsif ifreq = 1 and contType$ = "amplitude" and nan_out = 1
         atrf = 0
         atri = 0
         atrp = 0
         atrcip = 0
      endif
   endif
endfor
endproc