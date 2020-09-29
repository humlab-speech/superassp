form Compute a formant track
sentence SoundFile /Users/frkkan96/Desktop/kaa_yw_pb.wav
real BeginTime 0.0
real EndTime 0.0
real Time_step 0.0
real NO_formants 5.0
real MaxHzFormant 5500.0
real WindowLength 0.025
real PreEmph 50.0
sentence TrackOut /Users/frkkan96/Desktop/kaa_yw_pb.FormantTab
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
	soundPart = Extract part: beginTime, endTime, "Hanning", 1.0, 1
	selectObject: sound
	Remove
	sound = soundPart
endif 

selectObject: sound
noprogress To Formant (burg): time_step, nO_formants, maxHzFormant, windowLength, preEmph
Down to Table: 1, 1, 10, 1, 3, 1, 3, 1
Save as comma-separated file: trackOut$
Remove
