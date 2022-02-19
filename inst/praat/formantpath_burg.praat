form Compute a formant track using the iterative FormantPath functionality of Praat
sentence SoundFile /Users/frkkan96/Desktop/kaa_yw_pb.wav
real BeginTime 0.0
real EndTime 0.0
real Time_step 0.005
real Number_of_formants 5.0
real MaxHzFormant 5500.0
real WindowLength 0.025
real Pre_emphasis 50.0
real Ceiling_step_size 0.05
natural Number_of_steps_each_direction 4
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
noprogress To FormantPath (burg): time_step, number_of_formants, maxHzFormant, windowLength, pre_emphasis, ceiling_step_size, number_of_steps_each_direction
Extract Formant
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

