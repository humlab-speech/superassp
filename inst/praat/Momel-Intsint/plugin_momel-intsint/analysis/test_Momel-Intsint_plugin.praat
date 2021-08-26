#praat script
script_name$ = "test_Momel-Intsint_plugin.praat"
version$ = "2011-02-12"
date$ = date$()

#author Daniel Hirst
#email daniel.hirst@lpl-aix.fr

#purpose Test the different scripts available in the Momel-Intsint plugin
#	The script will successively
#	- Create recording directories for each .wav file in the parent directory "test1"
#	- Detect f0 for each file with optimised max and min f0
#	- Calculate Momel targets for each file
#	- Calculate Intsint labels for each file

pause Create recording directories for each .wav file in directory "test"
#	- Create recording directories for each .wav file in directory "test"
execute create_recording_directories.praat test .wav

pause Detect f0 for each sound in directory "test"
#	- Detect f0 for each file with optimised max and min f0
execute detect_f0.praat test .wav .Pitch yes 60 750 0.01 no yes

pause Calculate Momel targets for each file in directory "test"
#	- Calculate Momel targets for each file
execute calculate_momel_targets_extracts.praat  test .wav .Pitch .hz .momel .momel_auto .PitchTier
... no yes 60 750 0.01 "30 60 750 1.04 20 5 0.05" "-25 0.25 0.05" "#" "sound"

pause Press continue to calculate Intsint labels
#	- Calculate Intsint labels for each file
execute calculate_intsint_labels.praat test Momel Intsint IntsintMomel .momel .intsint .wav .TextGrid .PitchTier
... no yes 5
