#praat script correct_momel_targets

script_name$ = "correct_momel_targets.praat"
version$ = "[2011-02-12]"
date$ = date$()

#author: Daniel Hirst
#email: daniel.hirst@lpl-aix.fr

#purpose: treatment of recording directory 'name$' containing
# Sound file ('name$'.wav) and momel file ('name$'.momel)
# outputs corrected momel targets to 'name$.PitchTier' 
# and to 'name$'.momel overwriting earlier file

form Correct Momel targets
	comment To analyse several files consecutively
	comment  - click "Apply" for each one except the last  
	comment - then click "OK" for the last one 
	comment Leave empty to select with the browser:
	sentence sound_directory 
	word Sound_extension .wav
	word Momel_extension .momel
	word PitchTier_extension .PitchTier
	comment Default values
	natural minimum_f0 65
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

clearinfo
path_separator$ = "/"

if sound_directory$ = ""
	sound_directory$ = chooseDirectory$("Select your sub directory")
endif
index = rindex(sound_directory$, path_separator$)
name$ = right$(sound_directory$, length(sound_directory$)-index)


file_name$ = sound_directory$+path_separator$+name$
sound_file$ = file_name$+sound_extension$
momel_file$ = file_name$+momel_extension$
pitchTier_file$ = file_name$+pitchTier_extension$
log_file$ = file_name$+ ".log"

if fileReadable(sound_file$)
	printline treating file 'name$'
	mySound = Read from file... 'sound_file$'
	duration = Get total duration
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

	if  fileReadable(pitchTier_file$)
		myPitchTier = Read from file... 'pitchTier_file$'
		call treatment
		select mySound
		plus myPitchTier
		Remove
	elsif fileReadable(momel_file$)
		call create_pitchTier
		call treatment
		select mySound
		plus myPitchTier
		Remove
	else
		printline No file 'name$' with ['momel_extension$'] extension
	endif
else
	printline No file 'name$' with ['sound_extension$'] extension
endif

	
procedure create_pitchTier
	Create PitchTier... 'name$' 0 duration
	myPitchTier = selected("PitchTier")
	Read Strings from raw text file... 'momel_file$'
	myStrings = selected ("Strings")
	nStrings = Get number of strings
	for iString from 1 to nStrings
		select myStrings
		string$ = Get string... iString
		ms = extractNumber(string$,"")
		secs = ms/1000
		f0 = extractNumber(string$, " ")
		select myPitchTier
		Add point... secs f0
	endfor
endproc

procedure treatment
	select mySound
	myManipulation = To Manipulation... 0.01 min_f0 max_f0
	plus myPitchTier
	Replace pitch tier
	select myManipulation
	Edit
	pause Correct Momel targets then switch back to linear interpolation
	select myManipulation
	Extract pitch tier
	Write to text file... 'pitchTier_file$'
	nTargets = Get number of points
	filedelete 'momel_file$'
	for iTarget from 1 to nTargets
		time = Get time from index... iTarget
		time_ms = time*1000
		target = Get value at index... iTarget
		line$ = "'time_ms:0''tab$''target:0''newline$'"
		line$ >> 'momel_file$'
	endfor
	plus myManipulation
	Remove

	fileappend 'log_file$' 'name$''momel_extension$' and 'name$''pitchTier_extension$' corrected by 'script_name$' version 'version$' on 'date$' by 'script_name$' version 'version$''newline$'
endproc


#2011-02-19	Added check for version of Praat
#2011-02-12	No longer necessary to specify file system or Working Directory
#			File system is known from predefined variable and Parent directory can be 
#			selected via the browser (or specified as parameter)
#			Pitch max raised to 2.5*q75 to allow for expressive speech
#			Removed 'self modifying' feature which is no longer necessary.
# 2010:03:06 corrected bug in pause
# 2007:06:17  saving new values is optional
# 2007:02:24  script is self-modifying - new values are written as defaults
# 2007:01:31 corrected errors in log file
# 2006:12:20 treat a single file rather than batch
# 2006:12:12 working directory and system parameters defined by script set_working_directory and read from files in plugin
# 2006:11:20 skips "." and ".."
# 2006:10:29 first version
 
