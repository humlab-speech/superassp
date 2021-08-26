#praat script
script_name$ = "momel_single_file.praat"
version$ = "2011-09-23"
date$ = date$()

#author: Daniel Hirst
#email: daniel.hirst@lpl.univ-aix.fr

#purpose: single file call to momel

form Momel single file
	comment Do not change the following unless know what you are doing!
	sentence momel_parameters 30 60 750 1.04 20 5 0.05
endform

path_separator$ = "/"
if macintosh
	momel$ = "./momel_osx_intel"
#	momel$ = "./momel_osx_ppc"	### uncomment this and comment preceding for PPC version
elsif unix
	momel$ = "momel_linux"
elsif windows
	momel$ = "momel_win"
endif

pitch_extension$ = ".hz"
momel_extension$ =  ".momel"
output_folder$ = "temp"
createDirectory(output_folder$)

nPitches = numberOfSelected("Pitch")
if nPitches != 1
	pause Please select one pitch object
endif

myPitch = selected("Pitch")
name$ = selected$("Pitch")
duration = Get total duration

clearinfo
file_path$ = output_folder$+path_separator$+name$
pitch_file$ = file_path$+pitch_extension$
momel_file$ = file_path$+momel_extension$
myMatrix = To Matrix
Transpose
Save as headerless spreadsheet file... 'pitch_file$'
plus myMatrix
Remove
	select myPitch
endif

system 'momel$' > "'momel_file$'" 'momel_parameters$' < "'pitch_file$'"
Read Matrix from raw text file... 'momel_file$'
myMatrix = selected()
nRows = Get number of rows
Create PitchTier... 'name$' 0 'duration'
myPitchTier = selected()

for iRow from 1 to nRows
	select myMatrix
	time_ms = Get value in cell... iRow 1
	hz = Get value in cell... iRow 2
	select myPitchTier
	time_s = time_ms/1000
	Add point... time_s hz
endfor

select myMatrix
Remove
filedelete 'momel_file$'
filedelete 'pitch_file$'
select myPitchTier

exit

#2011-09-23  Replaced momel_osx  for Mac by momel_oxs_intel
#2011-02-21  Removed two lines left over from earlier version
#2011-02-19	No longer necessary to specify system
#			system is known from predefined variable
#			Removed 'self modifying' feature which is no longer necessary.
#2008:07:13	applies to selected Pitch object only now, uses temporary files in "temp"
#			added selection of System and autosave of default parameter
#2007:06:16 	adapted to apply to selected file
#2007:01:04 	working!
