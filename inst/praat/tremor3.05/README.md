# tremor.praat, version 3.05
tremor.praat, version 3.05 was built under Praat version 6.1.47.
tremor.praat is licensed under GNU GPL version 3 (see gnu_gpl_v3.txt).


## To run tremor.praat
Praat is required. Praat can be downloaded from http://www.praat.org.
Open tremor.praat within Praat and select *Run -> Run* from the *Script* window menu.


## Input
tremor.praat is thought for analyzing sustained vowels (optimally with an /a/-quality) that were produced as quasi-stationary (the speaker should intend to hold pitch and loudness constant for several seconds) as possible and that are available as Praat *Sound* objects or as .wav-files.


## Output
tremor.praat extracts the following measures...
...that capture frequency modulations:
1. frequency contour magnitude (FCoM),
2. (maximum) frequency tremor cyclicality (FTrC),
3. number of frequency modulations above thresholds (FMoN),
4. (strongest) frequency tremor frequency (FTrF),
5. frequency tremor intensity index (FTrI) at FTrF,
6. frequency tremor power index (FTrP) at FTrF,
7. frequency tremor cyclicality intensity product (FTrCIP) at FTrF,
8. frequency tremor product sum (FTrPS),
9. frequency contour harmonicity-to-noise ratio (FCoHNR),

...as well as amplitude modulations:
10. amplitude contour magnitude (ACoM),
11. (maximum) amplitude tremor cyclicality (ATrC),
12. number of amplitude modulations above thresholds (AMoN),
13. (strongest) amplitude tremor frequency (ATrF),
14. amplitude tremor intensity index (ATrI),
15. amplitude tremor power index (ATrP),
16. amplitude tremor cyclicality intensity product (ATrCIP),
17. amplitude tremor product sum (ATrPS), and
18. amplitude contour harmonicity-to-noise ratio (ACoHNR).


## Details
- CoM: ranging from 0 to 1, with the higher values indicating more reliable extraction of the following measures, technically the same as Praat's *intensity* within a Pitch object (created on a contour)
- TrC: ranging from 0 to 1, with the higher values indicating more cyclic/periodic modulations, technically the same as Praat's *strength* within a Pitch object (created on a contour), i.e. the auto-correlation coefficient at a certain frequency
- MoN: the number of modulations that show CoM-values greater than the contour magnitude threshold and TrC-values greater than the tremor cyclicality threshold
- TrF [Hz]: the frequency of the strongest (with highest TrC) modulation
- TrI [%, relative to the mean]: the magnitude/intensity of the strongest (with highest TrC) modulation, see references for more detail
- TrP: the "power" (in trems of perception, not in terms of physics) of the strongest (with highest TrC) modulation, technically TrP is TrI weighted by TrF, see references for more detail
- TrCIP [%]: an index weighting the TrI by the TrC (by simple multiplication)
- TrPS: the sum of TrCIP ofall MoN modulations
- CoHNR [dB]: overall harmonicity (to noise) of the contour

## Program modes
tremor.praat can be run in different modes:
1. Analysis mode: This mode is recommended as tremor praat is still experimental, missing a thorough pre-processing of the sounds. In this mode you can look into the processing steps and change the (initial pitch extraction) arguments if you think that something went wrong. A sound object is needed in the *Praat Objects* window. Results are given into the *Praat Info* window.
2. Simple mode: The default mode: A sound object is needed in the *Praat Objects* window. Results are given into the *Praat Info* window.
3. Run mode: This mode is to process more than one sound (file) with same arguments for all sounds using the Praat GUI. Results are output to a tab-separated text file (that is found by default in the tremor3.05 sub-directory *results*).
4. Console mode: This mode is to process more than one sound (file) with same arguments for all sounds from a console/shell/other programs. For calling tremor.praat for example from a Linux console you may type:
`/usr/bin/praat --run tremor3.05/tremor.praat 4 "./sounds/" 0.015 60 350 0.03 0.3 0.01 0.35 0.14 2 1.5 15 0.01 0.15 0.01 0.01 2`
This command would process all .wav-files that are found in the tremor3.05 sub-directory *sounds* using the default/standard arguments, just like set by default in the tremor3.05 form.
The result file (tab-separated text) can be found in the tremor3.05 sub-directory *results*. You can change this by changing the path or the name within coninout.praat (found in the sub-directory *procedures*).
Further details for calling Praat scripts can be found here: https://www.fon.hum.uva.nl/praat/manual/Scripting_6_9__Calling_from_the_command_line.html


## References
- Markus Brückl (2012): [Vocal Tremor Measurement Based on Autocorrelation of Contours.](http://www.isca-speech.org/archive/interspeech_2012/i12_0715.html) In: Proceedings of the [ISCA](http://www.isca-speech.org/) Conference [Interspeech '12](http://www.isca-speech.org/archive/interspeech_2012/), Portland (OR), 715–718.
- Markus Brückl, Alain Ghio, François Viallet (2015): [Measurement of Tremor in the Voices of Speakers with Parkinson's Disease.](http://icnlsp.org/IMG/pdf/-10.pdf) In: Proceedings of the [1st International Conference on Natural Language and Speech Processing (ICNLSP 2015)](http://icnlsp.org/spip.php?rubrique20), Algiers, 44–48.
- Markus Brückl, Elvira Ibragimova, Silke Bögelein (2017): [Acoustic Tremor Measurement: Comparing Two Systems.](http://www.brykl.de/tremMeasSys.pdf) In: [Proceedings](http://digital.casalini.it/9788864536071) of the [10th International Workshop on Models and Analysis of Vocal Emissions for Biomedical Applications (MAVEBA 2017)](http://maveba.dinfo.unifi.it/), Florence, 19–22.


## Please note 
that tremor.praat is still experimental: its reliability of extraction of tremor measures is strongly dependant on an accurate pitch extraction (first step, see *Analysis mode*).


## Please note further 
that there was at least one Praat version, namely 6.1.13, that yields wrong values for the tremor cyclicality (FTrC and AtrC) and derived measures (FTrCIP and ATrCIP, FTrPS and ATrPS) due to a not correct assignment of the maximum value of Matrix objects.
Version 6.1.05 is the last proofed version before with which tremor.praat is working properly and version 6.1.47 is the first proofed version after this "matrix bug".