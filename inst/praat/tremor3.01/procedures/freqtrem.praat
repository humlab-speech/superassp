###############TREMOR#################
# freqtrem.praat is a Praat[6.0.14] script (http://www.praat.org/) 
# that serves as a procedure within tremor.praat.
###############TREMOR#################
# Author: Markus Brückl (markus.brueckl@tu-berlin.de)
# Copyright 2011-2017 Markus Brückl
# License: GNU GPL v3 (http://www.gnu.org/licenses/gpl.html)
######################################

######################################
# Frequency Tremor Analysis
######################################
procedure ftrem
   To Pitch (cc)... ts minPi 15 yes silThresh voiThresh ocCo ocjCo vuvCo maxPi

if mode = 1
Edit
pause Pitch contour
endif

   numberVoice = Count voiced frames
if numberVoice = 0
   ftrc = undefined
   ftrf = undefined
   ftri = undefined
   ftrp = undefined
else

# because PRAAT only runs "Subtract linear fit" if the last frame is "voiceless" (!?):
# numberOfFrames+1 (1)
   numberOfFrames = Get number of frames
   x1 = Get time from frame number... 1
#   am_F0 = Get mean... 0 0 Hertz

   Create Matrix... ftrem_0 0 slength numberOfFrames+1 ts x1 1 1 1 1 1 0
   for i from 1 to numberOfFrames
      select Pitch 'name$'
      f0 = Get value in frame... i Hertz
      select Matrix ftrem_0
# write zeros to matrix where frames are voiceless
      if f0 = undefined
         Set value... 1 i 0
      else
         Set value... 1 i f0
      endif
   endfor

# remove the linear F0 trend (F0 declination)
   To Pitch
   Subtract linear fit... Hertz
   Rename... ftrem_0_lin

# undo (1)
   Create Matrix... trem 0 slength numberOfFrames ts x1 1 1 1 1 1 0
   for i from 1 to numberOfFrames
      select Pitch ftrem_0_lin
      f0 = Get value in frame... i Hertz
      select Matrix trem
# write zeros to matrix where frames are voiceless
      if f0 = undefined
         Set value... 1 i 0
      else
         Set value... 1 i f0
      endif
   endfor

   To Pitch
   am_F0 = Get mean... 0 0 Hertz

# normalize F0-contour by mean F0
   select Matrix trem
   Formula... (self-am_F0)/am_F0

# since zeros in the Matrix (unvoiced frames) become normalized to -1 but 
# unvoiced frames should be zero (if anything)
# write zeros to matrix where frames are voiceless
   for i from 1 to numberOfFrames
      select Pitch trem
      f0 = Get value in frame... i Hertz
      if f0 = undefined
         select Matrix trem
         Set value... 1 i 0
      endif
   endfor

# to calculate autocorrelation (cc-method):
   select Matrix trem
   To Sound (slice)... 1
# calculate Frequency of Frequency Tremor [Hz]
   To Pitch (cc)... slength minTr 15 yes tremMagThresh tremthresh ocFtrem 0.35 0.14 maxTr
   Rename... trem_norm

   ftrf = Get mean... 0 0 Hertz

# calculate frequency contour magnitude and cyclicality
   call cycli
   ftrm = trm
   ftrc = trc

# calculate Magnitude Indices of Frequency Tremor [%]
   contType$ = "frequency"
   show = 0
   call tremIntIndex
   ftri = tri
   ftrp = ftri * ftrf/(ftrf+1)

# calculate the product of Cyclicality and Intensity Indix (at the strongest found frequency tremor frequency)
   ftrcip = ftri * ftrc

# calculate (by cyclicality) weighted Sum of Intensity Indices at all found frequency tremor frequencies
# equals the sum of FTrCIP-values for all found frequency tremor frequencies
   call tremProdSum
   ftrps = tris

# clean up the Object Window
   select Pitch trem
# uncomment if only frequency tremor is to be analyzed:
#   plus Pitch 'name$'
   plus Matrix ftrem_0
   plus Pitch ftrem_0
   plus Pitch ftrem_0_lin
   plus Matrix trem
   plus Sound trem
   plus Matrix strengths
   plus Matrix frequencies
#   plus Pitch trem_norm
#   plus PointProcess Maxima
#   plus PointProcess Minima
   Remove

endif
endproc