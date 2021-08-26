#praat script
script_name$ = "reset_test_data.praat"
version$ = "2011-02-19"
date$ = date$()

#author Daniel Hirst
#email daniel.hirst@lpl-aix.fr

#purpose Reset the test data in the Momel-Intsint plugin

temporary_folder$ = "temp"
test_folder$ = "test"
nSounds=4
sound1$ = "English"
sound2$ = "French"
sound3$ = "Korean"
sound4$ = "Chinese"

list = Create Strings as file list... list 'temporary_folder$'/*
nFiles = Get number of strings
for iFile to nFiles
	file$ = Get string... iFile
	filedelete 'temporary_folder$'/'file$'
endfor
Remove


for iSound to nSounds
	sound$ = sound'iSound'$
	if fileReadable(test_folder$+"/"+sound$+"/"+sound$+".wav")
		Read from file... 'test_folder$'/'sound$'/'sound$'.wav
		Save as WAV file... 'test_folder$'/'sound$'.wav
		Remove
		list = Create Strings as file list... list 'test_folder$'/'sound$'/*
			nFiles = Get number of strings
			for iFile to nFiles
				file$ = Get string... iFile
				filedelete 'test_folder$'/'sound$'/'file$'
			endfor
			filedelete 'test_folder$'/'sound$'
		Remove
	endif
endfor

#Version history
#2011-02-19	First version
