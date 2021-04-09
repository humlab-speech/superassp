form Compute the intensity of a signal
sentence SoundFile /Users/frkkan96/Desktop/kaa_yw_pb.wav
real BeginTime 0.0
real EndTime 0.0
real Time_step 0.0
real Minimal_f0_Frequency 50.0
boolean Subtract_mean 1
word WindowShape Gaussian1
real RelativeWidth 1.0
sentence TrackOut /Users/frkkan96/Desktop/kaa_yw_pb.IntensityTab
endform

if fileReadable(soundFile$)
	sound = Read from file... 'soundFile$'
else
	exitScript("Could not read file 'soundFile$'")
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

selectObject: sound
noprogress To Intensity: minimal_f0_Frequency, time_step, subtract_mean
Down to IntensityTier
Down to TableOfReal
To Table: "dummy"
Remove column: "dummy"
Save as comma-separated file: trackOut$
Remove
