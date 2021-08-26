#praat script
scriptName$ = "make_phoneme_table.praat"
version$ = "2007:07:14"

#author Daniel Hirst
#email daniel.hirst@lpl-aix.fr

#purpose make table of phonemes and mean durations

form MakePhonemeTable
	sentence Directory /Users/daniel/Documents/Corpus/English/Eurom1-EN/files
	sentence Subdirectory fa
	word Means_table_extension .means
	word TextGrid_extension .TextGrid
	word Phoneme_tier phoneme
endform
clearinfo

if subdirectory$ != ""
	path$ = directory$+"/"+subdirectory$
else
	path$ = directory$
endif

myTable = Create Table with column names... 'subdirectory$' 1 phoneme duration
iTablePhoneme = 1

Create Strings as directory list... myList 'path$'
myStrings = selected("Strings")

nFolders = Get number of strings
printline folder 'path$' contains 'nFolders' recording folders

for iFolder to nFolders
	select myStrings
	name$ = Get string... iFolder
	textGrid_file$ = path$+"/"+name$+"/"+name$+textGrid_extension$
	if fileReadable(textGrid_file$)
		call treatment
	else
		printline sorry, can't read TextGrid file for 'name$'
	endif
endfor

select myStrings
Remove
select myTable
#first row is empty - can't create table with 0 rows
if iTablePhoneme > 1
	Remove row... 1
endif

Down to TableOfReal... phoneme
myTableOfReal = selected("TableOfReal")
To TableOfReal (means by row labels)... no
table_file$ = path$ + "/" + subdirectory$ + "-phonemes.means"
filedelete 'table_file$'
Write to headerless spreadsheet file... 'table_file$'
plus myTable
plus myTableOfReal
Remove

procedure treatment
	phoneme_tier = 0
	Read from file... 'textGrid_file$'
	myTextGrid = selected("TextGrid")
	nTiers = Get number of tiers
	printline treating 'phoneme_tier$' tier of file 'name$'

	for iTier to nTiers
		tier$ = Get tier name... iTier
		if tier$ = phoneme_tier$
			phoneme_tier = iTier
		endif
	endfor

	if phoneme_tier
		nPhonemes = Get number of intervals... phoneme_tier
		printline file 'name$' contains 'nPhonemes' phonemes
		for iPhoneme from 1 to nPhonemes
			iTablePhoneme = iTablePhoneme+1
			select myTextGrid
			phoneme$ = Get label of interval... phoneme_tier iPhoneme
			start_phoneme = Get starting point... phoneme_tier iPhoneme
			end_phoneme = Get end point... phoneme_tier iPhoneme
			duration = end_phoneme - start_phoneme
			select myTable
			Append row
			Set string value... iTablePhoneme phoneme 'phoneme$'
			Set numeric value... iTablePhoneme duration 'duration'
		endfor
		select myTextGrid
		Remove
	else
		printline No 'phoneme_tier$' tier for 'name$'
	endif
endproc

#2007:07:14	adapted to structure with working directory, subdirectory and recording directories
# 2007:02:21	modified
# 2005:05:12 	first version
