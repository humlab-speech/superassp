# synthTrem4096.praat is a Praat [6.0.29] script that generates 
# 4096 (may take some time!) wav-files (4096*282 KB!)
# into the directory, where this script resides.
# These files contain fully synthetical, tremulous /a/ vowels.
# These tremulous sounds served as test sounds in 
# Brückl, Ibragimova & Bögelein:
# Acoustic Tremor Measurement: Comparing Two Systems.
# 10th MAVEBA, Florence, 2017.
# Author: Markus Brückl (markus.brueckl@tu-berlin.de)
# License: GNU GPL > v3 (http://www.gnu.org/licenses/gpl.html)

slength = 3; sound duration [s]
amf = 200; mean frequency [Hz]
ama = 0.5; mean amplitude [Pa]
sts = 0.005; sythesis_time_step [s]

for ftstep from 0 to 3
   ftrf = 2.5+ftstep*3.5; frequency tremor frequency [Hz]
   for atstep from 0 to 3
      atrf = 3+atstep*4; amplitude tremor frequency [Hz]
      for ftistep from 0 to 3
         ftri = 0.01+ftistep*0.035; frequency tremor intensity index
         for atistep from 0 to 3
            atri = 0.02+atistep*0.045; amplitude tremor intensity index
            for fdstep from 0 to 3
               decf = fdstep*5; (linear) frequency decline [Hz/s]
               for adstep from 0 to 3
                  deca = adstep*0.05; (linear) amplitude decline [Pa/s]


nPoints = slength / sts; number of synthesis target points
bf = amf + (slength*decf)/2; frequency at beginning [Hz]
ba = ama + (slength*deca)/2; amplitude at beginning [Pa]


Create PitchTier... "tone" 0 slength
Add point... slength/2 amf

# source synthesis, Model: Rosenberg (1971); Klatt & Klatt (1990)
To Sound (phonation)... 48000 0.5 0.05 0.7 0.03 3 4 "no"

# filter synthesis
Create FormantGrid... "a" 0 slength 10 600 1100 60 50
Add formant point... 1 1.6 820
Remove formant points between... 1 1.49 1.51
Add formant point... 2 1.6 1350
Remove formant points between... 2 1.49 1.51
Add formant point... 3 1.6 2900
Remove formant points between... 3 1.49 1.51
Add formant point... 4 1.6 4000
Remove formant points between... 4 1.49 1.51
Add formant point... 5 1.6 4300
Remove formant points between... 5 1.49 1.51

select Sound tone
plus FormantGrid a
Filter

# tremor generation via overlap-and-add resynthesis
# frequency tremor generation:
To Manipulation... sts 75 600

Create PitchTier... "tremor" 0 slength
for iPoint from 0 to nPoints
   tPoint = iPoint * sts
   pitch = bf + ftri * amf * sin(ftrf*2*pi*tPoint) - decf * tPoint
   Add point... tPoint pitch
endfor

select Manipulation tone_filt
plus PitchTier tremor
Replace pitch tier

select Manipulation tone_filt
Get resynthesis (overlap-add)
Rename... tremor_filt


# amplitude tremor generation:
Create AmplitudeTier... "tremor" 0 slength
for iPoint from 0 to nPoints
   tPoint = iPoint * sts
   press = ba + atri * ama * sin(atrf*2*pi*tPoint) - deca * tPoint
   Add point... tPoint press
endfor

select Sound tremor_filt
plus AmplitudeTier tremor

# since "Multiply" not only multiplies, but scales (the maximum to 0.9),
# the definition of a mean amplitude is only meaningful in the context of 
# an amplitude decline -> set realtive in form
Multiply

# output tremor intensities in %
ftri = ftri*100
atri = atri*100

soundname$ = "'ftrf:1'_'atrf:0'_'ftri:0'_'atri:0'_'decf:0'_'deca:2'"
soundname$ = replace$(soundname$,".","-",4)

Save as WAV file... ./'soundname$'.wav

# undo tremor intensities in %
ftri = ftri/100
atri = atri/100

#Play

select PitchTier tone
plus Sound tone
plus FormantGrid a
plus Sound tone_filt
plus Manipulation tone_filt
plus PitchTier tremor
plus Sound tremor_filt
plus AmplitudeTier tremor
plus Sound tremor_filt_amp
endif

Remove


               endfor
            endfor
         endfor
      endfor
   endfor
endfor