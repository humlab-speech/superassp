minimum_version = 5204
minimum_version$ = "5.2.04"
if 'praatVersion' < minimum_version
	exit You are using version 'praatVersion$' of Praat. 
	... 'newline$'This plugin requires at least version 'minimum_version$'
	... 'newline$'Please update your version of Praat from http://www.praat.org
endif