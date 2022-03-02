form Compute f0 using the autocorrelation and cross-correlation  
	sentence SoundFile /Users/frkkan96/Desktop/kaa_yw_pb.wav
	real BeginTime 0.0
	real EndTime 0.0
	real Time_step 0.005
	real Window_length_(s) 0.040
	real Minimum_f0 75.0
	real Maximum_f0 600
	boolean Very_accurate 1
 	natural Number_of_cancidates 15
	real Silence_threshold 0.03
	real Voicing_threshold 0.45
	real Octave_cost 0.01
	real Octave_jump_cost 0.35
	real Voiced/unvoiced_cost 0.14
	boolean Only_correlation_methods 1
	real Minimum_filter_frequency_(Hz) 70.0
	real Maximum_filter_frequency_(Hz) 5000.0
	natural Number_of_filters 250
	real Maximum_frequency_components 1250.0
	natural Maximum_number_of_subharmonics 15
	real Compression_factor 0.84
	natural Number_of_points_per_octave 48
	word WindowType Gaussian1
	real RelativeWidth 1.0
	sentence TrackOut /Users/frkkan96/Desktop/kaa_yw_pb.f0Tab
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
noprogress To Pitch (cc)... 'time_step' 'minimum_f0' 'number_of_cancidates' 'very_accurate' 'silence_threshold' 'voicing_threshold' 'octave_cost' 'octave_jump_cost' 'voiced/unvoiced_cost' 'maximum_f0'
ccMatrix = To Matrix

selectObject: sound
noprogress To Pitch (ac)... 'time_step' 'minimum_f0' 'number_of_cancidates' 'very_accurate' 'silence_threshold' 'voicing_threshold' 'octave_cost' 'octave_jump_cost' 'voiced/unvoiced_cost' 'maximum_f0'
acMatrix = To Matrix

if only_correlation_methods == 0

	selectObject: sound
	noprogress To Pitch (SPINET)... 'time_step' 'window_length' 'minimum_filter_frequency' 'maximum_filter_frequency' 'number_of_filters' 'maximum_f0' 'number_of_cancidates'
	spinetMatrix = To Matrix

	selectObject: sound
	noprogress To Pitch (shs)... 'time_step' 'minimum_f0' 'number_of_cancidates' 'maximum_frequency_components' 'maximum_number_of_subharmonics' 'compression_factor' 'maximum_f0' 'number_of_points_per_octave'
	shsMatrix = To Matrix
endif 

select ccMatrix
plus acMatrix
Merge (append rows)

if only_correlation_methods == 0
	plus spinetMatrix
	Merge (append rows)
	plus shsMatrix
	Merge (append rows)
endif 

Transpose
To TableOfReal
To Table: "time"
Set column label (index): 2, "cc"
Set column label (index): 3, "ac"

if only_correlation_methods == 0

	Set column label (index): 4, "spinet"
	Set column label (index): 5, "shs"

endif

Save as comma-separated file: trackOut$
#Remove

