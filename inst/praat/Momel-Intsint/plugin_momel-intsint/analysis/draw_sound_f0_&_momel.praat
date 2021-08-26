#Praat script draw_sound_f0_&_momel.praat

#author:  Daniel Hirst <daniel.hirst@lpl-aix.fr>
#version:	2013-04-15

form Draw sound, f0 and momel
	comment Picture size (in cms)
	positive width 25
	positive height 15
	comment Extract
	real start 0
	real end 0 (= all)
	positive pitch_floor 60
	real pitch_ceiling 0 (= use automatic)
	word Pitch_colour Green
	word Spline_colour Red
	word Target_colour Black
	word Sound_colour Black
	positive Target_size 3 (=mm)
	boolean Draw_momel yes
endform

debug=0

cms = 1/2.54
nSounds = numberOfSelected("Sound")
nPitchTiers = numberOfSelected("PitchTier")

if draw_momel
	if nSounds = 0 or nPitchTiers = 0
		pause Select one Sound and one PitchTier
	endif
else
	if nSounds = 0
		pause Select one Sound
	endif
endif

mySound = selected("Sound")
name$ = selected$("Sound")
if draw_momel
	myPitchTier = selected("PitchTier")
endif

select mySound
duration = Get total duration

if end=0
	end = duration
endif

if pitch_ceiling = 0
	To Pitch... 0.01 60 700
	q25 = Get quantile... 0 0 0.25 Hertz
	q75 = Get quantile... 0 0 0.75 Hertz
	pitch_floor = q25 * 0.75
	pitch_ceiling = q75 *2
	Remove
endif
select mySound
myPitch = To Pitch... 0.01 pitch_floor pitch_ceiling
min_f0 = Get minimum... 0 0 Hertz None
max_f0 = Get maximum... 0 0 Hertz None

if draw_momel
	select myPitchTier
	nTargets = Get number of points
	for iTarget to nTargets
		hz = Get value at index... iTarget
		if hz < min_f0
			min_f0 = hz
		elsif hz > max_f0
			max_f0 = hz
		endif
	endfor
endif
			### Draw figure ###

Erase all
Select outer viewport... 0 width*cms 0 height*cms*0.8

Font size... 18
Solid line

	#Draw Pitch#

'pitch_colour$'
Line width... 2
min_f0 = floor(min_f0/25)*25
max_f0 = ceiling(max_f0/25)*25
select myPitch
Draw... start end min_f0 max_f0 no

if draw_momel
	#Draw Quadratic Spline#

	select myPitchTier
	myInterpolated = Copy... 'name$'
	Interpolate quadratically... 10 Hz
	myQSpline = To Pitch... 0.01 min_f0 max_f0

	'spline_colour$'
	Draw... start end min_f0 max_f0 no

	#Draw Targets#

	'target_colour$'
	Line width... 1
	select myPitchTier

	for iTarget to nTargets
		time = Get time from index... iTarget
		hz = Get value at index... iTarget

		if debug
			pause 'time:3' 'hz:0'
		endif

		if time > start and time < end
			Draw circle (mm)... time hz target_size
		endif
		endfor
endif

Marks left every... 1 50 yes yes yes
Text left... yes f_0 (Hz)

	#Draw Sound#

'sound_colour$'
select mySound
Select outer viewport... 0 width*cms height*cms*0.6 height*cms
Line width... 1
Draw... start end 0 0 no curve

#Add title etc.
Select outer viewport... 0 width*cms 0 height*cms
Text top... yes 'name$' 
Marks bottom every... 1 0.5 yes yes no
Text bottom... yes time (secs.) ['start:3'..'end:3']

select myPitch
if draw_momel
	plus myInterpolated
	plus myQSpline
endif
Remove
	
select mySound
if draw_momel
	plus myPitchTier
endif

#versionhistory

#2013-04-15		Made draw momel optional
#2011-11-02		First version


