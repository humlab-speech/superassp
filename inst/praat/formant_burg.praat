form Compute a formant track
sentence SoundFile /Users/frkkan96/Desktop/kaa_yw_pb.wav
real BeginTime 0.0
real EndTime 0.0
real Time_step 0.005
real Number_of_formants 5.0
real MaxHzFormant 5500.0
real WindowLength 0.025
real Pre_emphasis 50.0
boolean Track_formants 1
natural Number_of_tracks 3
natural Reference_F1 550
natural Reference_F2 1650
natural Reference_F3 2750
natural Reference_F4 3850
natural Reference_F5 4950
real Frequency_cost 1.0
real Bandwidth_cost 1.0
real Transition_cost 1.0 
word WindowShape Gaussian1
real RelativeWidth 1.0
word Spectrogram_window_shape Gaussian
real Spectrogram_resolution 40.0
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
	soundPart = Extract part: beginTime, endTime, windowShape$, relativeWidth, 1
	selectObject: sound
	Remove
	sound = soundPart
endif 

selectObject: sound
noprogress To Formant (burg): time_step, number_of_formants, maxHzFormant, windowLength, pre_emphasis
if track_formants == 1
	Track: number_of_tracks, reference_F1, reference_F2, reference_F3, reference_F4, reference_F5, frequency_cost, bandwidth_cost, transition_cost
	number_of_formants = min (number_of_formants, number_of_tracks)
endif

formantTab = Down to Table: 1, 1, 10, 1, 3, 1, 3, 1
noRows = Get number of rows
for f from 1 to number_of_formants
	lab$ = "L'f'(dB)"
	Append column: lab$
endfor

# Just make sure that we dont get edge effects
maxHz = maxHzFormant + 2000.0

selectObject: sound
noprogress To Spectrogram: windowLength, maxHz, time_step, spectrogram_resolution,  spectrogram_window_shape$
spectrogram = selected ("Spectrogram")

for r from 1 to noRows
	for f from 1 to number_of_formants
		outLab$ = "L'f'(dB)"
		inLab$ = "F'f'(Hz)"
		selectObject: formantTab
		currFreq = Get value: r, inLab$
		if currFreq != undefined
			currTime = Get value: r, "time(s)"
			selectObject:  spectrogram
			currIntensity = Get power at: currTime, currFreq
		else
			currIntensity = undefined
		endif
		selectObject: formantTab
		Set numeric value: r, outLab$, currIntensity
	endfor
endfor


Save as comma-separated file: trackOut$
Remove

