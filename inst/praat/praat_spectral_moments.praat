form Compute the spectral moments 1-4 
sentence SoundFile /Users/frkkan96/Desktop/a1.wav
real BeginTime 0.0
real EndTime 0.0
real WindowLength 0.005
real Maximum_frequency_(Hz) 0.0
real Time_step 0.005
real Frequency_step 20.0
real Power 2
word WindowShape Gaussian1
real RelativeWidth 1.0
sentence TrackOut /Users/frkkan96/Desktop/a1.mom
endform


if fileReadable(soundFile$)
	sound = Read from file... 'soundFile$'
else
	exitScript("Could not read file 'soundFile$'")
endif

selectObject: sound
dur = Get total duration
sr = Get sampling frequency

#Pre-compute some things based on the sound, if needed
if maximum_frequency == 0.0
	maximum_frequency = sr / 2
endif

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

outTab = Create Table with column names: "outTable", 0, "Time CenterOfGravity SD Skewness Kurtosis"

selectObject: sound
noprogress To Spectrogram: windowLength, maximum_frequency, time_step, frequency_step, "Gaussian"
spect = selected ("Spectrogram")
noFrames = Get number of frames

for frame from 1 to noFrames

	select spect

	currTime = Get time from frame number: frame
	spec = To Spectrum (slice): currTime
	

	cog = Get centre of gravity: power
	sd = Get standard deviation: power
	skew = Get skewness: power
	kurt = Get kurtosis: power
	selectObject: spec
	Remove
	selectObject: outTab
	Append row
	row = Get number of rows
	Set numeric value: row, "Time", currTime
	Set numeric value: row, "CenterOfGravity", cog
	Set numeric value: row, "SD", sd
	Set numeric value: row, "Skewness", skew
	Set numeric value: row, "Kurtosis", kurt

endfor

select outTab
Save as comma-separated file: trackOut$
Remove
