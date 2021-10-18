## getoq.praat

## James Kirby <j.kirby@ed.ac.uk>
## 23 February 2017

## Given an EGG signal (maybe smoothed, maybe not)... compute Oq two ways

#########################
## NOTES AND ISSUES
#########################

## 1. Both the dEGG and Howard methods involve the determination of *periods*, 
## measured from (dEGG) closing peak to closing peak. So the first point in the 
## PointProcess needs to be a closing peak, as does the last. These are computed
## from the dEGG signal, so the only difference between the methods is in the 
## Oq (not f0) computation;, in particular, with respect to how the start of the
## opening phase is determined.  It is very much up to the USER to insure that 
## the periods have been correctly selected. Both the FIRST and LAST points in
## the PointProcess should therefore be CLOSING (i.e. positive) peaks. 

## Note that it is not enough to simply remove the point of uncertain dEGG 
## opening peaks, because the algorithm expects closing-opening peak pairs. 
## Therefore, there is a separate step where whole periods are removed.

## 2. At present, nothing intelligent is done about multiple peaks: they
## are not detected, nor is there an option given to do anything about them.
## Maxima and minima are currently determined using the built-in Praat function
## To PointProcess (periodic, peaks)... 

## 3. It is possible to discard particular periods manually, a la peakdet. 
## This entails keeping track of the periods when adding/deleting points, and 
## then entering these in the next step. If periods are discarded, e.g. due to 
## the presence of multiple opening peaks, they are removed from BOTH methods' 
## matrices (or more precisely, set to 0).

## This is likely to be problematic only if the Oq of particular, individual
## periods is ultimately of interest. In most cases this probably won't 
## matter, because you will interpolate between the missing values at the analysis 
## stage. However, this ultimately depends on your use case.

## One advantage of this method - i.e., ensuring that there are points where 
## (you think) there should be pulses, then explicitly removing whole periods
## at a later stage - is that this preserves information by indicating in the 
## output file  that there 'was' a pulse there, but that it was manually removed 
## (i.e., set to 0). 

## Dependencies
include splitstring.praat
include smooth.praat
include peakdet.praat
include degg.praat
include howard.praat
include writelns.praat
include plotoq.praat
include exclude.praat

procedure getoq: .manualCheck

    ## Clear Picture window
    Erase all

    ## Extract Lx signal
    name$ = selected$ ("Sound", 1)
    Extract one channel... eggChan
  
    ## Rename
    ch$ = selected$ ("Sound", 1)
    select Sound 'name$'
    Remove
    select Sound 'ch$'
    Rename... 'name$'
    name$ = selected$ ("Sound", 1)

    ## Create standarized version of file name for plotting purposes
    @splitstring: name$, separator$
    plotName$ = ""
    for i from 1 to (splitstring.strLen - 1)
    	plotName$ = plotName$ + splitstring.array$[i] + "-"
    endfor
    plotName$ = plotName$ + splitstring.array$[splitstring.strLen]
   
    ## Filter
    Copy: "'name$'_filtered"
    Filter (pass Hann band)... passFrequency 0 smoothHz

    ## Smooth the filtered signal
    Copy: "'name$'_fsmooth"
    @smooth: k
    Formula... 'smooth.formula$'

    ## Get opening and closing peaks based on dEGG signal
    @peakdet

    #############
    ## Main loop
    #############
    
    ## Set flag 
    .findOQ = 1

    while .findOQ <> 0

        ###############################
        ## Find peaks
        ###############################
       
        ## find total number of opening and closing peaks
        select PointProcess 'name$'_degg_both

        ## if we could assume the entire file was relevant than we could just say
        # nb_peaks = Get number of points
        ## but given we are potentially interested only in a sub-region...
        
        ## get the first point *following the onset* of the region
        first_point  = Get high index... start_time
        ## get the last point *preceding the offset* of the region
        last_point = Get low index... end_time
        nb_peaks = (last_point - first_point) + 1

        ## get number of *periods* (close->close)
        nb_periods = (nb_peaks / 2) - 1

        ####################################
		## only continue if there is a 
		## non-zero number of periods!!
        ####################################

		if nb_periods > 0

			####################################
			## Get Oq using dEGG-only method 
			####################################

			@degg

			####################################
			## Get Oq using Howard's method 
			####################################
			
			@howard

			#####################
			## Other procedures
			#####################
			
			#@skewness

			#####################
			## Plot and check
			#####################
		   
			if .manualCheck <> 0
				@plotoq: plotName$

				beginPause: "Manual check options"
					comment: "Do you want to manually add/delete any points? (1=yes)"
					integer: ".manualCheck", .manualCheck
				endPause: "Continue", 1
		
				################################
				## Add/remove points if desired 
				################################
				
				if .manualCheck == 1
					## Remove existing Matrix objects since new ones will be created
					select Matrix 'name$'_degg
					plus Matrix 'name$'_howard
					Remove

					## Now add/remove points
					select PointProcess 'name$'_degg_both
					plus Sound 'name$'_degg
					View & Edit
					editor: "PointProcess 'name$'_degg_both"
						Zoom: start_time, end_time
					endeditor
					pause Add missing peaks/delete spurious points
					editor: "PointProcess 'name$'_degg_both"
						Close
					endeditor
				endif
			   
				################################
				## Exclude periods, if desired
				################################
				
				if .manualCheck == 0
					## Call procedure to allow user to exclude periods
					@exclude 
					## We're done with this file. Set flag and write results to disk
					.findOQ = 0
				endif

			else
				## Not doing manual check; get me out of this loop!
				.findOQ = 0
			endif 

		endwhile
  
		##################################
		## Write to file and save objects
		##################################
		@writelns

		###################
		## Clean up 
		###################
		select Matrix 'name$'_degg
		plus Matrix 'name$'_howard
		Remove

	## otherwise just write the single line
	else
		@writelns
	endif

endproc 
