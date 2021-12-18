## writetable.praat: write output file and PointProcess object

## James Kirby <j.kirby@ed.ac.uk>
## 23 February 2017


procedure writetable
   

 
    ## Note that because Praat Matrix objects can't contain characters,
    ## cols 2-5 of rows that have been 'discarded' are set to 0.
    ## This could potentially be handled more gracefully here, or can
    ## be dealt with at the analysis stage (e.g., by replacing 0s with 
    ## NAs in your statistical analysis software).
    

	## if you have at least one period:
	if nb_periods > 0
		exitScript: "Found some periods"
    	for i from 1 to nb_periods-1
			myFileLine$ = ""

    	    ## DEGG
    	    select Matrix 'name$'_degg
    	    currPeriod = Get value in cell... i 1
    	    pstart = Get value in cell... i 2
    	    pend = Get value in cell... i 3
    	    f0 = Get value in cell... i 4
    	    degg_oq = Get value in cell... i 5

			
			#myFileLine$ = myFileLine$ + name$ + "," + lingVars$ + printIntervalLabel$ + ",'currPeriod','pstart','pend','f0','degg_oq'"

    	    ## Howard
			select Matrix 'name$'_howard
			howard_oq = Get value in cell... i 5
			#myFileLine$ = myFileLine$ + ",'howard_oq'"

			## If you want to add other measures, just continue to append them here

			## Now, write to file
			#appendFileLine: "'directory$''outfile$'", "'myFileLine$'"

			select outTab
			Append row
			row = Get number of rows

			if currPeriod <> undefined
				Set numeric value: row, "period", currPeriod
			endif


		endfor

	else
		exitScript: "No periods were found. Please check your settings."
	endif

  	## Save PointProcess object
   	select PointProcess 'name$'_degg_both
    Save as text file... 'directory$''name$'_degg_both.PointProcess
	select outTab
	Save as semicolon-separated file: "'outfile'"
endproc
