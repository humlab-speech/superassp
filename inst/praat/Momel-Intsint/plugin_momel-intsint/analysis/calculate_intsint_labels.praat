#praat script
script_name$ = "calculate_intsint_labels.praat"
version$ = "[2011-02-19]"

#author Daniel Hirst
#email daniel.hirst@lpl-aix.fr

#purpose batch treatment of folders name$ each containing
#name$.momel and name$.TextGrid or name$.wav or name$.PitchGrid
#Calculate intsint code from Momel targets 
#

form Select Parent directory for treatment
	comment Leave empty to select with the browser
	sentence Parent_directory 
	word Momel_tier Momel
	word Intsint_tier Intsint
	word IntsintMomel_tier IntsintMomel
	word Momel_extension .momel
	word Intsint_extension .intsint
	word Sound_extension .wav
	word TextGrid_extension .TextGrid
	word PitchTier_extension .PitchTier
	comment Overwrite existing Intsint file and TextGrid tiers?
	boolean overwrite 0
	comment Print info during treatment?
	boolean verbose 1
	comment default values (only change if you know why!)
	natural Intsint_header 5
endform

minimum_version = 5204
minimum_version$ = "5.2.04"
if 'praatVersion' < minimum_version
	exit You are using version 'praatVersion$' of Praat. 
	... 'newline$'This plugin requires at least version 'minimum_version$'
	... 'newline$'Please update your version of Praat from http://www.praat.org
endif

epsilon_time = 0.0001
clearinfo
path_separator$ = "/"

if parent_directory$ = ""
	parent_directory$ = chooseDirectory$("Select your sub directory")
endif

myFolders = Create Strings as directory list... myList 'parent_directory$'
nFolders = Get number of strings
if verbose
	printline 'parent_directory$' contains 'nFolders' folders.
endif

for iFolder from 1 to nFolders
	select myFolders
	name$ = Get string... iFolder
	if name$ != "." and name$ != ".."
		if verbose
			printline ['iFolder'] 'name$'
		endif
		file_name$ = parent_directory$+path_separator$+name$+path_separator$+name$
		momel_file$ = file_name$+momel_extension$
		intsint_file$ = file_name$+intsint_extension$
		textGrid_file$ = file_name$+textGrid_extension$
		sound_file$ = file_name$+sound_extension$
		pitchTier_file$ = file_name$+pitchTier_extension$
		if fileReadable(momel_file$)
			if not(fileReadable(intsint_file$)) or overwrite
				call treatment
			else
				printline 'name$''intsint_extension$' already exists
			endif
		else
			printline 'name$''momel_extension$' is not readable
		endif
	endif
endfor

select myFolders
Remove
exit

procedure treatment
	system perl intsint.pl "'momel_file$'"
	call treat_TextGrid
	if verbose
		printline momel intsint intsintMomel tiers: 'momel_tier' 'intsint_tier' 
		...'intsintMomel_tier'
	endif
	if fileReadable (intsint_file$)
		myIntsint = Read Strings from raw text file... 'intsint_file$'
		nStrings = Get number of strings
		nTargets = nStrings-intsint_header
		last_time = 0
		for iTarget from 1 to nTargets
			iString = iTarget+intsint_header
			select myIntsint
			string$ = Get string... iString
			time = extractNumber(string$, "")
			if time = last_time
				time = time + epsilon_time
			endif
			last_time = time
			intsint$ = extractWord$(string$," ")
			momel = extractNumber(string$,intsint$)
			intsintMomel = extractNumber(string$," 'momel' ")
			select myTextGrid
			Insert point... momel_tier time 'momel'
			Insert point... intsint_tier time 'intsint$'
			Insert point... intsintMomel_tier time 'intsintMomel'	
		endfor 
		Write to text file... 'textGrid_file$'
		plus myIntsint
		Remove
	else
		exit There is no 'intsint_file$'
		... 'newline$'
		... 'newline$' Perhaps Perl is not installed on your computer.
		... 'newline$' You can download Perl from: 
		... 'newline$'
		... 'newline$''tab$'http://www.activestate.com/downloads
	endif
endproc

procedure treat_TextGrid
	momel_tier = 0
	intsint_tier = 0
	intsintMomel_tier = 0
	if fileReadable(textGrid_file$) and not overwrite
#read existing TextGrid
		Read from file... 'textGrid_file$'
		myTextGrid = selected("TextGrid")
		nTiers = Get number of tiers

		for iTier from 1 to nTiers
		if verbose
			print tier 'iTier' 
		endif
		tier_name$ = Get tier name... iTier
		if verbose
			print ['tier_name$']
		endif
			if tier_name$ = momel_tier$
				if verbose
					print  ....changing
				endif
				momel_tier = iTier
				Remove tier... iTier
				Insert point tier... iTier 'momel_tier$'
			elsif tier_name$ = intsint_tier$
				if verbose
					print  ....changing
				endif
				intsint_tier = iTier
				Remove tier... iTier
				Insert point tier... iTier 'intsint_tier$'		
			elsif tier_name$ = intsintMomel_tier$
				if verbose
					print  ....changing
				endif
				intsintMomel_tier = iTier
				Remove tier... iTier
				Insert point tier... iTier 'intsintMomel_tier$'		
			endif
				if verbose
					printline
				endif
		endfor 
		call add_tiers
	elsif fileReadable(sound_file$)
#create new TextGrid
		Read from file... 'sound_file$'
		mySound = selected("Sound")
		To TextGrid... "'momel_tier$' 'intsint_tier$' 'intsintMomel_tier$'" 'momel_tier$' 
		...'intsint_tier$' 'intsintMomel_tier$'
		momel_tier = 1
		intsint_tier = 2
		intsintMomel_tier = 3
		nTiers = 3
		myTextGrid = selected ("TextGrid")
		select mySound
		Remove
		select myTextGrid
	elsif fileReadable(pitchTier_file$)
		Read from file... 'pitchTier_file$'
		duration = Get total duration
		Remove
		Create TextGrid... 0 duration "'momel_tier$' 'intsint_tier$' 'intsintMomel_tier$'"  'momel_tier$' 
		...'intsint_tier$' 'intsintMomel_tier$'
		momel_tier = 1
		intsint_tier = 2
		intsintMomel_tier = 3
		nTiers = 3
		myTextGrid = selected ("TextGrid")
		select myTextGrid
	else 
		printline Cannot create TextGrid without 'name$''sound_extension$' or 'name$''pitchTier_extension$'
	endif
endproc

procedure add_tiers
	if momel_tier = 0
		nTiers = nTiers+1
		momel_tier = nTiers
		Insert point tier... momel_tier 'momel_tier$'
		if verbose
			printline inserting tier 'nTiers'
		endif
	endif
	if intsint_tier = 0
		nTiers = nTiers+1
		intsint_tier = nTiers
		Insert point tier... intsint_tier 'intsint_tier$'
		if verbose
			printline inserting tier 'nTiers'
		endif
	endif
	if intsintMomel_tier = 0
		nTiers = nTiers+1
		intsintMomel_tier = nTiers
		Insert point tier... intsintMomel_tier 'intsintMomel_tier$'
		if verbose
			printline inserting tier 'nTiers'
		endif
	endif
endproc


#2013-03-21  	Added possiblity to create TextGrid from duration of PitchTier
#2011-02-19	Added check for version of Praat
#2011-02-12	No longer necessary to specify file system or Working Directory
#			File system is known from predefined variable and Parent directory can be 
#			selected via the browser (or specified as parameter)
#			Removed 'self modifying' feature which is no longer necessary.
# 2007:10:02   check directories to exclude "." and ".." (on Windows)
#			ensure time of target is different to that of preceding target
# 2007:06:17   subdirectory declared as sentence and filenames quoted to allow spaces in names
# 2007:03:31  	saving arguments as default is optional
# 2007:02:24  	script self-modifying - new arguments are written as default
# 2006:01:10		corrected errors in procedure create TextGrid
# 2006:01:09		corrected errors in procedure add_tiers
# 2006:12:12 	working directory and system parameters defined by script set_working_directory 
#				and read from files in plugin
# 2006:11:01 	first version
