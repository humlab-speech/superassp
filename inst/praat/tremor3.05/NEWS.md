# new in version...

## 3.05
- corrected the run and console output of ACoHNR (which was wrongly the same as for FCoHNR)
- written a short parameter explanation in README.md
- improved the console mode / added an argument for sound-file input

## 3.04
- re-written all procedures in the new syntax, now using object IDs for object selection and removal

## 3.03
- corrected the newly implemented contour harmonicity measures
- inserted a missing "selectObject", because without runtime errors occur sometimes

## 3.02
- output of contour harmonicities FCoHNR and ACoHNR
- implementation of a choice to output indeterminate values
- correction of tremor measures that were erroneous/inconsistent due to unexplainable (to me) and yet unexplained (see Praat-Users-List, Msg #2401, #5727, and #8679) automatic sortation / ranking order of pitch candidates within a (single frame) Pitch object
- implementation of the "console" mode to run tremor.praat

## 3.01
- invention of new measures that combine cyclicality and intensity

## 2.07
- output of measures for tremor cyclicality (measures of tremor periodicity, technically speaking: the contours' auto-correlation coefficients, or as Praat names it in Pitch objects: the frequency's "strength")
- output of measures of the contours' (mean?) magnitudes (technically speaking: the contours' "intensity" as given by the Pitch object that results from analyzing the contours like a (real) sound with respect to their cyclicality.

## 2.06
- improved the resampling at a constant rate (time step) of the amplitude contour (per period)

## 2.05
- removed a bug causing the script to stop if the sound is considered voiceless at the beginnig or at the end -- e.g. because of wrong pitch range

## 2.04
- removed a bug in amptrem.praat causing a (the more) raised "zero"-level in  amplitude contour (the lower the mean sound intensity was) and therefore raised ATrI and ATrP values

## 2.03
- second way to extract an AmplitudeTier from Sound & PointProcess: RMS per pitch period

## 2.02
- modularization of the script into separate files

## 2.01
- the version presented in 2012 at Interspeech

## 1.00
- the version used in Brückl (2011)