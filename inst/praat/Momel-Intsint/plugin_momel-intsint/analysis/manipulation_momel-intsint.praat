#praat script
 
script_name$ = "manipulate_momel_intsint.praat"
date$ = date$()
version$ = "[2015-03-19]"
#author: Daniel Hirst
#email: daniel.hirst@lpl.univ-aix.fr

#purpose: Create modification file with PitchTier derived from Insint coding
#requires: One Sound file and one TextGrid with same name (given in Form)
#		   TextGrid should have Momel, Intsint and IntsintMomel tiers


form Define parameters
	comment Leave empty to select with the browser
	sentence Sound_directory 
	word Sound_extension .wav
	word TextGrid_extension .TextGrid
	comment Tier names:
	word Momel_tier Momel
	word Intsint_tier Intsint
	word IntsintMomel_tier IntsintMomel
	positive time_step 0.01
	natural minimum_f0 60
	natural maximum_f0 750
endform

#2011-02-19	Added check for version of Praat
minimum_version = 5204
minimum_version$ = "5.2.04"
if 'praatVersion' < minimum_version
	exit You are using version 'praatVersion$' of Praat. 
	... 'newline$'This plugin requires at least version 'minimum_version$'
	... 'newline$'Please update your version of Praat from http://www.praat.org
endif


if windows
	path_separator$ = "\"
else
	path_separator$ = "/"
endif


if sound_directory$ = ""
	sound_directory$ = chooseDirectory$("Select your parent directory")
endif
index = rindex(sound_directory$, path_separator$)
sound_name$ = right$(sound_directory$, length(sound_directory$)-index)

clearinfo


sound_file$ = sound_directory$+path_separator$+sound_name$+sound_extension$
textGrid_file$ = sound_directory$+path_separator$+sound_name$+textGrid_extension$
min_f0_file$ = sound_directory$+path_separator$+sound_name$+".min_f0"
max_f0_file$ = sound_directory$+path_separator$+sound_name$+".max_f0"

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

mySound = Read from file... 'sound_file$'
myTextGrid = Read from file... 'textGrid_file$'


select myTextGrid
duration = Get total duration
nTiers = Get number of tiers

for iTier from 1 to nTiers
	tier$ = Get tier name... iTier
	if tier$ = momel_tier$
		momel_tier = iTier
	elsif tier$ = intsint_tier$
		intsint_tier = iTier
	elsif tier$ = intsintMomel_tier$
		intsintMomel_tier = iTier
	endif
endfor

printline Momel = tier 'momel_tier'; Intsint = tier 'intsint_tier'; IntsintMomel = tier 'intsintMomel_tier'
nTargets = Get number of points... momel_tier
myMomel = Create PitchTier... "Momel" 0 duration
myIntsintMomel = Create PitchTier... "IntsintMomel" 0 duration

for iTarget from 1 to nTargets
	select myTextGrid
	time = Get time of point... momel_tier iTarget
	momel$ = Get label of point... momel_tier iTarget
	momel = 'momel$'
	intsintMomel$ = Get label of point... intsintMomel_tier iTarget
	intsintMomel = 'intsintMomel$'
	select myMomel
	Add point... time momel
	select myIntsintMomel
	Add point... time intsintMomel
endfor

select mySound
To Manipulation... 'time_step' 'min_f0' 'max_f0'
momel_manipulation = selected("Manipulation")

Rename... Momel
Copy... IntsintMomel

intsint_manipulation = selected("Manipulation")
Edit
plus myIntsintMomel
Replace pitch tier
select momel_manipulation
Edit
plus myMomel
Replace pitch tier

exit

#2015-03-19	Added check for system to specify path_separator$ as "\" for Win-dos.
#2011-02-19	Added check for version of Praat
#2011-02-14	No longer necessary to specify file system or Working Directory
#			Recording directory can be selected via the browser (or specified as parameter)
#			Removed 'self modifying' feature which is no longer necessary.
# 2007:03:31  saving arguments as default now optional
# 2007:02:24  script now self-modifying - new values are written as defaults
# 2006:12:12 working directory and system parameters defined by script set_working_directory and read from files in plugin
# 2006:11:01 First version
