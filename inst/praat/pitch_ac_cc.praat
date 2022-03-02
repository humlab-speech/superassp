form Compute f0 using the autocorrelation and cross-correlation  
	sentence SoundFile /Users/frkkan96/Desktop/kaa_yw_pb.wav
	real BeginTime 0.0
	real EndTime 0.0
	real Time_step 0.005
	real Minimum_f0 75.0
	real Maximum_f0 600
	boolean Very_accurate 1
#	real Maximum_period_factor 1.3
#   real Maximum_amplitude_factor 1.6
	real Silence_threshold 0.03
	real Voicing_threshold 0.45
	real Octave_cost 0.01
	real Octave_jump_cost 0.35
	real Voiced/unvoiced_cost 0.14
	word WindowType Gaussian1
	real RelativeWidth 1.0
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
	soundPart = Extract part: beginTime, endTime, windowType$, relativeWidth, 1
	selectObject: sound
	Remove
	sound = soundPart
endif 

selectObject: sound
noprogress To Pitch (cc)... 0.0 'minimum_f0' 15 'very_accurate' 'silence_threshold' 'voicing_threshold' 'octave_cost' 'octave_jump_cost' 'voiced/unvoiced_cost' 'maximum_f0'
ccMatrix = To Matrix

selectObject: sound
noprogress To Pitch (ac)... 0.0 'minimum_f0' 15 'very_accurate' 'silence_threshold' 'voicing_threshold' 'octave_cost' 'octave_jump_cost' 'voiced/unvoiced_cost' 'maximum_f0'
acMatrix = To Matrix

select ccMatrix
plus acMatrix
Merge (append rows)
Transpose
To TableOfReal
To Table: "time"
Set column label (index): 2, "cc"
Set column label (index): 3, "ac"


Save as comma-separated file: trackOut$
#Remove

