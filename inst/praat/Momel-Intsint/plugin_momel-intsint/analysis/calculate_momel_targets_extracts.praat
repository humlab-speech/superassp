#praat script
script_name$ = "calculate_momel_targets_extracts.praat"
version$ = "[2013-03-21]"
date$ = date$()

#author: Daniel Hirst
#email: daniel.hirst@lpl.univ-aix.fr

#purpose: batch analysis of pitch files
#	The parent directory contains at least one recording directory (eg text01, text02…)
#	Each recording directory contains at least one sound file (eg text01.wav) and one pitch file (eg text01.hz)

#this script calculates momel targets <ms, Hz> (calling an external C program)
#default extensions: f0 : ".hz"; momel targets: ".momel"
#before calling momel, the script splits the sound file into pause separated segments
#this improves the target detection before and after pauses (cf Hirst et al. 2007 - Interspeech)

#before running this script, the operating system and the working directory 
#containing the subdirectory should be selected using the script
#		set_working_directory.praat

form Select Parent directory for treatment
	comment Leave empty to select with the browser
	sentence Parent_directory 
	word Sound_extension .wav
	word Pitch_extension .Pitch
	word F0_extension .hz
	word Momel_extension .momel
	word Momel_auto_extension .momel_auto
	word PitchTier_extension .PitchTier
	boolean Overwrite_existing_target_points 0
	boolean Display_info_on_treatment 1
	comment default values: Only change these if you know what you're doing!
	natural left_Pitch_range 60
	natural right_Pitch_range 750
	positive pitch_step 0.01
	sentence momel_parameters 30 60 750 1.04 20 5 0.05
	sentence silence_parameters -25 0.25 0.05 
	word silence_label #
	word sound_label sound
endform

debug=0

minimum_version = 5204
minimum_version$ = "5.2.04"
if 'praatVersion' < minimum_version
	exit You are using version 'praatVersion$' of Praat. 
	... 'newline$'This plugin requires at least version 'minimum_version$'
	... 'newline$'Please update your version of Praat from http://www.praat.org
endif

if parent_directory$ = ""
	parent_directory$ = chooseDirectory$("Select your sub directory")
endif

overwrite = overwrite_existing_target_points
verbose = display_info_on_treatment
clearinfo

if macintosh
	momel$ = "./momel_osx_intel"
#	momel$ = "./momel_osx_ppc"		
elsif unix
	momel$ = "momel_linux"
elsif windows
	momel$ = "momel_win"
else
	exit Sorry your system does not seem to be supported. Contact daniel.hirst@lpl-aix.fr
endif

silence_tier = 1
path_separator$ = "/"

myFolders = Create Strings as directory list... myFolders 'parent_directory$'
nFolders = Get number of strings

if verbose
	printline 'parent_directory$' contains 'nFolders' folders
endif

for iFolder to nFolders
	select myFolders
	name$ = Get string... iFolder
	if verbose
		printline ['iFolder'] 'name$'
	endif
	file_name$ = parent_directory$+path_separator$+name$+path_separator$+name$
	sound_file$ = file_name$+sound_extension$
	pitch_file$ = file_name$ + pitch_extension$
	momel_file$ = file_name$+momel_extension$
	momel_auto_file$ = file_name$+momel_auto_extension$
	pitchTier_file$ = file_name$+pitchTier_extension$
	if fileReadable(momel_file$) and not overwrite
		printline 'name$''momel_extension$' already exists
	else
		if fileReadable(sound_file$)
			call treat_sound
		else
			printline 'name$''sound_extension$' is not readable
		endif
	endif
endfor

select myFolders
Remove
exit


procedure treat_sound
	if verbose
		printline treating file 'sound_file$'
	endif
	min_f0_file$ = file_name$+".min_f0"
	max_f0_file$ = file_name$+".max_f0"
	log_file$ = file_name$+".log"

	if fileReadable(min_f0_file$)
		min_f0$ < 'min_f0_file$'
		min_f0 = 'min_f0$'
	else
		min_f0 = left_Pitch_range
	endif
	if fileReadable(max_f0_file$)
		max_f0$ < 'max_f0_file$'
		max_f0 = 'max_f0$'
	else
		max_f0 = right_Pitch_range
	endif

#create temporary folder in analysis folder of plugin if it doesn't already exist
	system_nocheck mkdir temp
#remove any existing files in this folder
	myOld_files = Create Strings as file list... Old_files temp/*
	nOld_files = Get number of strings

	for iOld_file to nOld_files
		old_file$ = Get string... iOld_file
		filedelete temp/'old_file$'
	endfor

#Read Sound and create Pitch, Matrix, TextGrid and PitchTier
	mySound = Read from file... 'sound_file$'
	duration = Get total duration

	if  fileReadable(pitch_file$)
	myPitch = Read from file... 'pitch_file$'
	else
		myPitch = To Pitch... pitch_step min_f0 max_f0
	endif
		
	myMatrix = To Matrix
	select mySound
	myTextGrid = To TextGrid (silences)... 'min_f0' 'pitch_step' 'silence_parameters$' 'silence_label$' 'sound_label$'
	nIntervals = Get number of intervals... silence_tier
if debug
	pause : 'name$' 'duration'
endif
	myPitchTier = Create PitchTier... "'name$'" 0 duration

#treat non-silent extracts of sound
	nExtracts = 0

	for iInterval from 1 to nIntervals
		select myTextGrid
		label$ = Get label of interval... silence_tier iInterval
		if label$ = sound_label$
			nExtracts = nExtracts+1
			call treat_extract
		endif
	endfor

	select myPitchTier
	call write_to_momel
	fileappend "'log_file$'" 'name$''momel_extension$' and 'name$''momel_auto_extension$' 
	...created by 'script_name$' version 'version$' on 'date$''newline$'
	plus myOld_files
	plus myPitch
	plus mySound
	plus myMatrix
	plus myTextGrid
	Remove
endproc

procedure treat_extract
	extract$ = "extract_"+"'nExtracts'"
	f0_extract$ = "temp/"+extract$+f0_extension$
	momel_extract$ = "temp/"+extract$+momel_extension$
	if verbose
		printline Treating  'extract$'
	endif
	start = Get start point...  silence_tier iInterval
	iStart = ceiling(start / pitch_step)
	if iStart <1
		iStart = 1
	endif
	end = Get end point... silence_tier iInterval
	iEnd = floor(end / pitch_step)
	if verbose
		printline Extract 'nExtracts' (interval 'iInterval'): ['iStart':'iEnd']
	endif

#write values of pitch to temporary f0 file
	select myMatrix
	nColumns = Get number of columns
	if iEnd > nColumns
		iEnd = nColumns
	endif	
	for iPitch from iStart to iEnd
		pitch = Get value in cell... 1 iPitch
		fileappend "'f0_extract$'" 'pitch''newline$'
	endfor

#call momel
	system 'momel$' >"'momel_extract$'" 'momel_parameters$' <"'f0_extract$'"
#Add targets from extract to PitchTier
	myStrings = Read Strings from raw text file... 'momel_extract$'
	nStrings = Get number of strings

	for iString from 1 to nStrings
		select myStrings
		string$ = Get string... iString
		ms = extractNumber(string$,"")
		if ms = undefined
			printline String ['iString'] ('string$') doesn't contain a number
		else
			secs = ms/1000
			f0 = extractNumber(string$, " ")
			if f0 > max_f0
				f0 = max_f0
			elsif f0 < min_f0
				f0 = min_f0
			endif
			time = secs+start
			if time <0
				time = 0
			elsif time > duration
				time = duration
			endif	
			select myPitchTier
			Add point... time f0
		endif
	endfor ; iString

	select myStrings
	Remove
endproc

procedure write_to_momel
	filedelete "'momel_file$'"
	filedelete "'momel_auto_file$'"
	filedelete "'pitchTier_file$'"
	select myPitchTier
	Write to text file... 'pitchTier_file$'
	nTargets = Get number of points

	for iTarget from 1 to nTargets
		time = Get time from index...  iTarget
		time_ms = time*1000
		target = Get value at index... iTarget
		fileappend "'momel_file$'" 'time_ms:' 'target:''newline$'
		fileappend "'momel_auto_file$'" 'time_ms:' 'target:''newline$'
	endfor ; iTarget
endproc

# [2013-03-21]	Corrected bug that crashed when file name contains a space
# [2011-09-23]	Replaced mac_osx version of momel by intel version "momel_osx_intel"
# [2011-02-19]	Added check for version of Praat
# [2011-02-12]	No longer necessary to specify file system or Working Directory
#			File system is known from predefined variable and Parent directory can be 
#			selected via the browser (or specified as parameter)
#			Removed 'self modifying' feature which is no longer necessary.
#			No longer need to keep separate f0 and Pitch files 
#			f0 only used for f0_extract as input to Momel
# [2008-07-13]	corrected a 'bug' which assumed Windows understood the system command 'rm -r'
#			replaced this with Praat function filedelete in loop
#			modified numbering of extracts
#			now only removes objects created by script
# [2007-09-29] 	tidied up form to make it smaller
# [2007-03-07]	saving arguments as defaults is now optional
# [2007-02-24] 	script now self-modifying - new arguments are written as defaults
# [2007-02-05] 	added list subdirectories to form
# [2007-01-22] 	Renumbered versions of script
#			Remove intermediate files during treatment to avoid Praat error from too many objects 
#			Added constraint on targets: cannot be >max_f0 or <min_f0
#			Corrected indexing error when first interval contains sound
# [2007-01-07]		Momel targets calculated on pause separated portions of signal
# [2007-01-04] 	default pitch extension changed to .hz to avoid confusion with Praat .Pitch files
# [2006-12-12] 	working directory and system parameters defined by script set_working_directory and read from files in plugin
# [2006-12-06] 	path_separator and momel version selection from choice of system in form
# [2006-06-02] 	plugin version
# [2004-11-17] 	Cyril Auran's Praat implementation
