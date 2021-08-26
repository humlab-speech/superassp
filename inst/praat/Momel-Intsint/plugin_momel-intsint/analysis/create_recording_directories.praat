#praat script 

script_name$ = "create_recording_directories.praat"
version$ = "[2011-02-19]"
date$ = date$()

#author: Daniel Hirst
#email: daniel.hirst@lpl-aix.fr

#purpose: Convert directory of sound files into subdirectories of recording folders
# 		one for each sound with same name

#before running this script, the operating system and the working directory 
#	should be selected using the script:
#		set_working_directory.praat

form Select directory for treatment
	comment Leave empty to select with the browser
	sentence Parent_directory 
	word Sound_extension .wav
endform

#2011-02-19	Added check for version of Praat
minimum_version = 5204
minimum_version$ = "5.2.04"
if 'praatVersion' < minimum_version
	exit You are using version 'praatVersion$' of Praat. 
	... 'newline$'This plugin requires at least version 'minimum_version$'
	... 'newline$'Please update your version of Praat from http://www.praat.org
endif

path_separator$ = "/"
clearinfo

if parent_directory$ = ""
	parent_directory$ = chooseDirectory$("Select your directory of sound files")
endif

Create Strings as file list... sounds 'parent_directory$'/*'sound_extension$'
mySounds = selected("Strings")
nSounds = Get number of strings

for iSound to nSounds
	select mySounds
	sound$ = Get string... iSound
	if not(startsWith(sound$, "."))
		sound_file$ = parent_directory$+path_separator$+sound$
		recording$ = sound_file$-sound_extension$
		printline ['iSound'] 'recording$'
		Read from file... 'sound_file$'
		filedelete 'sound_file$'
		createDirectory (recording$)
		Write to WAV file... 'recording$''path_separator$''sound$'
		Remove
	endif
endfor
	
select mySounds
Remove
exit

#version history
#[2011-02-19]		Added check for version of Praat
#				Replaced system call mkdir by Praat function createDirectory
#[2011-02-12]		No longer necessary to specify file system or Working Directory
#				File system is known from predefined variable and Parent directory can be 
#				selected via the browser (or specified as parameter)
#[2007-11-06]		Added check to avoid treating backup files beginning with "."
#[2007-10-29]
#[2007-09-28]
#[2007-05-27]
#[2007-05-09]