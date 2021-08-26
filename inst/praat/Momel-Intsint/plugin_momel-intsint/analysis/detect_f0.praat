#Praat script
script_name$ = "detect_f0.praat"
version$ = "[2011-02-05]"
date$ = date$()

#author: 	Daniel Hirst <daniel.hirst@lpl-aix.fr>


# Purpose: batch detection of f0 of sound files in parent directory
#	each sound and analysis file in subdirectory of working directory is stored in folder 'name$'
#	default extensions:
#	- sound file : [.wav]
#	- f0 file: [.hz] (one value Hz per line, 0 for unvoiced)
#	date, script and version stored in 'name$'.log
#   parameters saved as 'name$'.pitchStep, ~.min_f0, ~.max_f0

form Select parent directory for treatment
	comment To analyse several directories  
	comment  - click "Apply" for each one (except the last, for which click "OK")
	comment 
	comment Leave empty to select with the browser
	sentence Parent_directory 
	word Sound_extension .wav
	word Pitch_extension .Pitch
	boolean automatic_max_and_min 1
	natural Minimum_f0  60
	natural Maximum_f0 750
	real pitch_step 0.01
	comment Overwrite existing pitch files ?
	boolean overwrite no
	comment Print info on analysis
	boolean verbose yes
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

if parent_directory$ = ""
	parent_directory$ = chooseDirectory$("Select your parent directory")
endif

myFolders = Create Strings as directory list... myList 'parent_directory$'

nFolders = Get number of strings

if verbose
	clearinfo
	printline 'parent_directory$' contains 'nFolders' folders
endif

for iFolder to nFolders
	select myFolders
	name$ = Get string... iFolder
	if name$ != "." and name$ != ".."
		if verbose
			printline ['iFolder'] 'name$'
		endif
		file_path$ = parent_directory$+path_separator$+name$+path_separator$+name$
		sound_file$ = file_path$+sound_extension$
		pitch_file$ = file_path$+pitch_extension$
		if fileReadable(sound_file$)
			if not(fileReadable(pitch_file$)) or overwrite
				call treatment
			elsif verbose
				printline file ['pitch_extension$'] already exists for 'name$'
			endif
		else
			printline cannot read 'name$''sound_extension$'
		endif
	endif
endfor

select myFolders
Remove

procedure treatment
	mySound = Read from file... 'sound_file$'

	if automatic_max_and_min
		call calculate_min_max_f0
	else
		min_f0 = minimum_f0
		max_f0 = maximum_f0
	endif

	if verbose
		printline treating file 'name$' with min = 'min_f0', max = 'max_f0'
	endif

	select mySound
	myPitch = To Pitch... 'pitch_step' 'min_f0' 'max_f0'
	Write to text file... 'pitch_file$'
	plus mySound
	Remove

# save parameters and log
	pitch_step$ = "'pitch_step'"
	min_f0$ = "'min_f0'"
	max_f0$ = "'max_f0'"
	pitchStepFile$ = file_path$+".pitch_step"
	minF0File$ = file_path$+".min_f0"
	maxF0File$ = file_path$+".max_f0"
	logfile$ = file_path$+".log"
	pitch_step$ > 'pitchStepFile$'
	min_f0$ > 'minF0File$'
	max_f0$ > 'maxF0File$'
	fileappend "'logfile$'" 'name$''pitch_extension$' and 'f0_extension$' created by 'script_name$' version 'version$' on 'date$''newline$'
endproc

procedure calculate_min_max_f0
#  estimate of newMinF0 as 0.75 * quantile 25 (first quartile)
# and newMaxF0 as 4 * newMinF0

#  rounded to higher (resp. lower) 10
	To Pitch... 'pitch_step' 'minimum_f0' 'maximum_f0'
	.q25 = Get quantile... 0.0 0.0 0.25 Hertz
	min_f0 = 10*floor((0.75*.q25)/10)
	max_f0 = min_f0*4
	Remove
endproc

#2015-01-23	Set max_f0 to 4*min_f0
#2011-02-19	Added check for version of Praat
#2011-02-12	No longer necessary to specify file system or Working Directory
#			File system is known from predefined variable and Parent directory can be 
#			selected via the browser (or specified as parameter)
#			Pitch max raised to 2.5*q75 to allow for expressive speech
#			Removed 'self modifying' feature which is no longer necessary.
#2011-02-05 	save Pitch file (.Pitch) and F0 file (.Hz)
#2008-07-06	corrected default parameters
#2007-10-04 	logfile$ quoted to allow spaces in name
#2007-06-17 	subdirectory declared as sentence to allow spaces in name
#			saving values of arguments as new standards is optional
#		   	only removes objects created by the script
#2007-02-24 	script is self-modifying - new values are written as defaults
#2007-01-07 	tidied up
#2007-01-04 	default pitch extension changed (again!) to .hz to avoid confusion with Praat .Pitch files
#2006-12-12 	working directory, system, path separator 
#			defined by script set_working_directory.praat 
#			values read from parameter files stored in plugin
#2006-10-31 	sound and pitch extension put as parameters in form
#2006-06-02 	changed extension of f0 file to .pitch to avoid confusion with
#           		.f0 files in Eurom1 which contain time and f0 couples
#2006-05-27 	corrected format - no time values in .f0 file
#2006-04-19 	first version - adapted from batch.praat version 2006-04-18
