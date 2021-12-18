## shelldet.praat: wrapper script to get Oq values for multiple files in a single directory from command line
#
## James Kirby <j.kirby@ed.ac.uk>
## last update: 17 August 2018

## If something goes wrong, you can stop the script and pick up where you 
## left off, by noting the last file in the Strings list that was correctly 
## processed. However be sure to rename your output file, or rename it 
## something new when you restart, or you will overwrite your previous data.

#########################
## USER-DEFINED VARIABLES
#########################

## - directory: path to EGG files
## - textgrid: path to TextGrids, if applicable (default: same as directory$)
## - outfile: name of output file (saved in directory$)
## - extension: file extension for EGG files (.wav, .egg...)
## - eggChan: channel number of EGG signal
## - intervalTier, intervalLabel, intervalNum: used to specify a specific
##   portion of the file to edit/extract from. if intervalNum <> 0, this 
##   will take precedence over intervalLabel. if intervalLabel == "" and
##   intervalNum == 0, entire file will be processed.

## - minF0, maxF0: minimum and maximum pitch values
## - k: used to calculate window size for smoothing.
##   k = 0 is same as no smoothing.
##   k = 2 > 5-point window; k = 3 > 7-point window; etc.
## - threshold: Howard's method threshold (default: 3/7)
form File info
#   comment Full path to EGG files
   sentence Input_File /Users/frkkan96/Desktop/egg/Lx.wav
	real BeginTime 0.0
	real EndTime 0.0	
    integer eggChan 2
#    comment Minimum and maximum f0 thresholds
    integer minF0 75
    integer maxF0 600
#    comment k: Smoothing window size parameter (points on each side)
    integer k 10
#    comment Threshold for Howard's method
    real    threshold 3/7
#    comment Filter frequency cutoff
    integer passFrequency 40
#    comment Filter cutoff smoothing
    integer smoothHz 20
#    comment Manually edit points and periods?
    boolean manualCheck 0
#    comment Use existing PointProcess files, if available?
    boolean useExistingPP 0
#    comment Invert signal (if your EGG has closed=down for some reason)
    boolean invertSignal 0
	 word WindowShape rectangular
	 real RelativeWidth 1.0
#    comment Name of output file (written to same path as EGG files)
    word outfile /Users/frkkan96/Desktop/Lx.csv
endform

## including getoq.praat includes everything else
include getoq.praat


## Create output table 

outTab = Create Table with column names: "outTab", 0, "label period start end egg_f0 degg_oq howard_oq"

if fileReadable(input_File$)
	sound = Read from file: input_File$
	soundFile$ = selected$("Sound")
	plotName$ = soundFile$
else
	exitScript("Could not read file 'input_File$'")
endif

selectObject: sound
dur = Get total duration

# Check that start or end times should be condidered, and that they are within 
# ok limits.
if  ( beginTime > 0.0 or endTime > 0.0 ) and (beginTime >= 0.0 and endTime <= dur)
	
	selectObject: sound
	# Preserve times so that start end end record may be computed later
	soundPart = Extract part: beginTime, endTime, windowShape$, relativeWidth, 1
	selectObject: sound
	Remove
	sound = soundPart
endif 

#Prepare for running the rest of the script
start_time = beginTime
end_time = endTime
intervalNum = 0
intervalLabel$ = ""
directory$ = ""

## invert signal if necessary
if invertSignal
  Formula... -self
endif

@getoq: manualCheck
  

#select all
#minus Strings list

#Remove

#clearinfo

#printline All done.
