#Praat script automatic_min_max_f0.praat
#
#version:	2016-11-21
#author 	Daniel Hirst
#email:	daniel.hirst@lpl-aix.fr

#purpose: calculate F0 using two passes with automatic estimation of optimal max and min
#		min f0 is 0.75 * 1st quartile of f0 distribution
#		max f0 is min f0 * max pitch span
#		values rounded down/up to nearest 10.

#requires: one Sound file (selected)

form Automatic min/max f0
	positive Maximum_pitch_span 1.5 (= octaves)
endform

clearinfo
nSounds = numberOfSelected("Sound")

if nSounds != 1
	pause Please select one Sound then press 'Continue'
endif

mySound = selected()
To Pitch... 0.01 60 750
q25 = Get quantile... 0 0 0.25 Hertz
Remove
select mySound
min_f0 = floor(q25 * 0.75/10)*10
max_f0 = ceiling((min_f0 * 2^maximum_pitch_span)/10)*10
printline detecting f0 with minimum pitch 'min_f0:' maximum pitch 'max_f0:'
To Pitch... 0.01 min_f0 max_f0

# Version History
# [2016-11-21]		added  variable pitch-span to control range
# [2011-02-16]		max_pitch raised to 2.5*q75 to allow for expressive speech (was 1.5*q75)
# [2008-07-11]		