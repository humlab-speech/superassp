#Praat script: Create microprosodic profile


form Create microprosodic profile
	sentence file_path /Users/dan/Documents/Corpus/English/Eurom1-EN/example
	word folder Eurom1-example
	word f0_extension .hz
	word momel_extension .PitchTier
	word microprosodic_profile_extension .mpp
	positive time_step 0.01
endform

minimum_f0 = 60
maximum_f0 = 700
file_name$ = file_path$ + "/" + folder$ + "/" + folder$
min_f0_file$ = file_name$+".min_f0"
max_f0_file$ = file_name$+".max_f0"

if fileReadable(min_f0_file$)
	min_f0$ < 'min_f0_file$'
	min_f0 = 'min_f0$'
else
	min_f0 = minimum_f0
endif
if fileReadable(max_f0_file$)
	max_f0$ < 'max_f0_file$'
	max_f0 = 'max_f0$'
else
	max_f0 = maximum_f0
endif

f0_file$ = file_name$+f0_extension$
myMatrix = Read Matrix from raw text file... 'f0_file$'
nRows = Get number of rows
duration = (nRows + 4)*time_step
myPitchTier = Create PitchTier... 'folder$' 0 duration
for iPoint from 1 to nRows
	select myMatrix
	time = (iPoint+4)*time_step
	pitch = Get value in cell... iPoint 1
	select myPitchTier
	Add point... time pitch
endfor

momel_file$ = file_name$ + momel_extension$
myMomelPitchTier = Read from file... 'momel_file$'
Interpolate quadratically... 4 Semitones
myMomelPitch = To Pitch... time_step min_f0 max_f0
myProfile = Create PitchTier... 'folder$' 0 duration
for iPoint from 1 to nRows
	time = (iPoint+4)*time_step
	select myMatrix
	pitch = Get value in cell... iPoint 1
	select myMomelPitch
	momelPitch = Get value in frame... iPoint Hertz
	microprosody = pitch/momelPitch
	select myProfile
	Add point... time microprosody
endfor

select myProfile
Write to text file... 'file_name$''microprosodic_profile_extension$'
plus myMatrix
plus myPitchTier
plus myMomelPitchTier
plus myMomelPitch
Remove



