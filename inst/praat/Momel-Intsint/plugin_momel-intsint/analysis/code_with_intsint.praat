# Praat script: code_with_intsint.praat
# Author: Daniel Hirst <daniel.hirst@lpl-aix.fr>
# Version: [2013-03-21]

# Purpose: code Momel pitch targets with INTSINT labels
# Requires: one PitchTier (selected) containing Momel target points

momel_extension$ = ".momel"
intsint_extension$ = ".intsint"
pitchTier_extension$ = ".PitchTier"
textGrid_extension$ = ".TextGrid"

myPitchTier = selected("PitchTier")
nTargets = Get number of points

name$ = selected$("PitchTier")
folder$ = "temp/"+name$

createDirectory (folder$)
intsint_file$ = folder$+"/"+name$+intsint_extension$
momel_file$ = folder$+"/"+name$+momel_extension$
pitchTier_file$ = folder$+"/"+name$+pitchTier_extension$
textGrid_file$ = folder$+"/"+name$+textGrid_extension$

Save as text file... 'pitchTier_file$'

for iTarget to nTargets
	time = Get time from index... iTarget
	time_ms = time * 1000
	hz = Get value at index... iTarget
	output_line$ = "'time_ms:3'"+tab$+"'hz'"+newline$
	output_line$ >> 'momel_file$'
endfor

execute calculate_intsint_labels.praat temp Momel Intsint IntsintMomel .momel .intsint .wav .TextGrid .PitchTier
... no yes 5

Read from file... 'textGrid_file$'
Read Strings from raw text file... 'intsint_file$'
deleteFile(intsint_file$)
deleteFile(momel_file$)
deleteFile(pitchTier_file$)
deleteFile(textGrid_file$)
deleteFile(folder$)

# Version history:

# [2013-03-21]	First version
