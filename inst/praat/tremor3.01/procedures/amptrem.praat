###############TREMOR#################
# amptrem.praat is a Praat[6.0.14] script (http://www.praat.org/) 
# that serves as a procedure within tremor.praat.
###############TREMOR#################
# Author: Markus Brückl (markus.brueckl@tu-berlin.de)
# Copyright 2011-2017 Markus Brückl
# License: GNU GPL v3 (http://www.gnu.org/licenses/gpl.html)
######################################

######################################
# Amplitude Tremor Analysis
######################################
procedure atrem

select Sound 'name$'
plus Pitch 'name$'
To PointProcess (cc)
numbOfGlotPoints = Get number of points
if numbOfGlotPoints < 3
   atrc = undefined
   atrf = undefined
   atri = undefined
   atrp = undefined
else
   if amplitude_extraction_method = 2
      select Sound 'name$'
      plus PointProcess 'name$'_'name$'
# amplitudes are RMS per period -- not intensity maxima ?? -- no! unclear, Praat help missing
      To AmplitudeTier (period)... 0 0 0.0001 0.02 1.7
      numbOfAmpPoints = Get number of points
   elsif amplitude_extraction_method = 1
# NEW 2014-12-13: corrected misinterpretation of Praat-function "To AmplitudeTier (period)"
      numbOfAmpPoints = numbOfGlotPoints - 1
      Create AmplitudeTier... 'name$'_'name$'_'name$' 0 'slength'
      for iAmpPoint from 1 to numbOfAmpPoints
         select PointProcess 'name$'_'name$'
         perStart = Get time from index... iAmpPoint
         perEnd = Get time from index... iAmpPoint+1
         select Sound 'name$'
         rms = Get root-mean-square... perStart perEnd
# very seldomly (with bad pitch settings) it occurs that perStart and perEnd are nearer 
# than sampling period -> rms would be undefined
         if rms = undefined
            samplPer = Get sampling period
            rms = Get root-mean-square... perStart-samplPer perEnd+samplPer
         endif
         select AmplitudeTier 'name$'_'name$'_'name$'
         Add point... ('perStart'+'perEnd')/2 rms
      endfor
   endif
######################################

if mode = 1
Edit
pause Amplitude contour
endif

# since bad pitch extraction may result in not even one amplitude point
ampPointN = Get number of points
if ampPointN = 0
   atrc = undefined
   atrf = undefined
   atri = undefined
   atrp = undefined
else

# from here on out: prepare to autocorrelate AmplitudeTier-data
# sample AmplitudeTier at (constant) rate ts
# to be able to -- automatically -- read Amp. values...
   Down to TableOfReal

# to enable autocorrelation of the Amp.-contour: ->Matrix->Sound

   Create Matrix... atrem_nlc 0 slength numberOfFrames+1 ts x1 1 1 1 1 1 0
# from here on out: get the mean of (the curve of) the amplitude contour in each frame
   for iframe from 1 to numberOfFrames
      select Pitch 'name$'
      f0 = Get value in frame... iframe Hertz
# determine (the time of) fixed interval borders for the resampled amplitude contour
         t = (iframe-1) * ts + x1
         tl = t - ts/2
         tu = t + ts/2
# get the indices of the amplitude points surrounding around these borders
         select AmplitudeTier 'name$'_'name$'_'name$'
         loil = Get low index from time... tl
         hiil = Get high index from time... tl
         loiu = Get low index from time... tu
         hiiu = Get high index from time... tu
# if the sound is unvoiced the amplitude is not extracted
      if f0 = undefined
         select Matrix atrem_nlc
         Set value... 1 iframe 0
# if the amplitude contour has not begun yet...
      elsif loil = 0
         select Matrix atrem_nlc
         Set value... 1 iframe 0
# ...or is already finished the amplitude is not extracted
      elsif hiiu = numbOfAmpPoints + 1; 
         select Matrix atrem_nlc
         Set value... 1 iframe 0
      else
         select TableOfReal 'name$'_'name$'_'name$'
         lotl = Get value... loil 1; time value of Amp.Point before tl in the PointProcess [s]
         druck_lol = Get value... loil 2; amplitude value before tl in the PointProcess [Pa, ranged from 0 to 1]
         hitl = Get value... hiil 1
         druck_hil = Get value... hiil 2; amplitude value after tl in the PointProcess
         lotu = Get value... loiu 1
         druck_lou = Get value... loiu 2; amplitude value before tu in the PointProcess
         hitu = Get value... hiiu 1; time value after tu in the PointProcess
         druck_hiu = Get value... hiiu 2; amplitude value after tu in the PointProcess
# caculate (linearly interpolated) pressure/amplitude at the borders
         druck_tl = ((hitl-tl)*druck_lol + (tl-lotl)*druck_hil) / (hitl-lotl)
         druck_tu = ((hitu-tu)*druck_lou + (tu-lotu)*druck_hiu) / (hitu-lotu)

         nPinter = hiiu - 1 - loil; = loiu - loil; = hiiu - hiil; number of amp.-points between tl and tu
         if nPinter = 0; loil = loiu; hiil = hiiu
            druck_mean = (druck_tl + druck_tu) / 2
         else
            tlinter = tl
            plinter = druck_tl
            sumtdruck = 0
            for iinter from 1 to nPinter
               tuinter = Get value... loil+iinter 1
               puinter = Get value... loil+iinter 2
               deltat = tuinter - tlinter
               tdruck_iinter = deltat*(plinter+puinter)/2
               sumtdruck += tdruck_iinter
               tlinter = tuinter
               plinter = puinter
            endfor
            deltat = tu - tlinter
            tdruck_iinter = deltat*(plinter+druck_tu)/2
            sumtdruck += tdruck_iinter
            druck_mean = sumtdruck / ts
         endif

         select Matrix atrem_nlc
         Set value... 1 iframe druck_mean
      endif
   endfor

# because PRAAT classifies frequencies in Pitch objects <=0 as "voiceless" and 
# therefore parts with extreme INTENSITIES would be considered as "voiceless"
# (irrelevant) after "Subtract linear fit" (1)
# "1" is added to the original Pa-values (ranged from 0 to 1) -- not to the voiceless parts
   select Matrix atrem_nlc
   for i from 1 to numberOfFrames+1
      grms =  Get value in cell... 1 i
      if grms > 0
         Set value... 1 i grms+1
      endif
   endfor

# remove the linear amp.-trend (amplitude declination)
   To Pitch
   Rename... hilf_lincorr

   Subtract linear fit... Hertz
   Rename... atrem
   am_Int = Get mean... 0 0 Hertz
   am_Int = am_Int - 1

# undo (1)... and normalize Amp. contour by mean Amp.
   To Matrix
   for i from 1 to numberOfFrames+1
      grms =  Get value in cell... 1 i
      if grms > 0
         Set value... 1 i (grms-1-am_Int)/am_Int
      endif
   endfor

# remove last frame, undo (2)
   Create Matrix... trem 0 slength numberOfFrames ts x1 1 1 1 1 1 0
   for iframe from 1 to numberOfFrames
      select Matrix atrem
      spring = Get value in cell... 1 iframe
      select Matrix trem
      Set value... 1 iframe spring
   endfor

# to calculate autocorrelation (cc-method)
   To Sound (slice)... 1
   To Pitch (cc)... slength minTr 15 yes tremMagThresh tremthresh ocAtrem 0.35 0.14 maxTr
   Rename... trem_norm

   atrf = Get mean... 0 0 Hertz

# calculate amplitude contour magnitude and cyclicality
   call cycli
   atrm = trm
   atrc = trc

# calculate Magnitude Indices of Amplitude Tremor [%]
   contType$ = "amplitude"
   show = 0
   call tremIntIndex
   atri = tri
   atrp = atri * atrf/(atrf+1)

# calculate the product of Cyclicality and Intensity Indix (at the strongest found amplitude tremor frequency)
   atrcip = atri * atrc

# calculate (by cyclicality) weighted Sum of Intensity Indices at all found amplitude tremor frequencies
# equals the sum of ATrCIP-values for all found amplitude tremor frequencies
   call tremProdSum
   atrps = tris

endif
endif

# clean up the Object Window
   select Pitch 'name$'
   plus PointProcess 'name$'_'name$'
if numbOfGlotPoints >= 3
   plus AmplitudeTier 'name$'_'name$'_'name$'
if ampPointN > 0
   plus TableOfReal 'name$'_'name$'_'name$'
   plus Matrix atrem_nlc
   plus Pitch hilf_lincorr
   plus Pitch atrem
   plus Matrix atrem
   plus Matrix trem
   plus Sound trem
   plus Matrix strengths
   plus Matrix frequencies
#   plus Pitch trem_norm
#   plus PointProcess Maxima
#   plus PointProcess Minima
endif
endif
   Remove

endproc