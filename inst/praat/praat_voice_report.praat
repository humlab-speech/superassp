

form Compute Praat voice report parameters
	sentence SoundFile /Users/frkkan96/Desktop/aaa.wav
	real StartTime 0.0
	real EndTime 0.0
	real SelectionOffset 0.0
	real SelectionLength 2.0
	real Minimum_f0 75.0
	real Maximum_f0 600
	real Maximum_period_factor 1.3
   real Maximum_amplitude_factor 1.6
	real Silence_threshold 0.03
	real Voicing_threshold 0.45
	real Octave_cost 0.01
	real Octave_jump_cost 0.35
	real Voiced/unvoiced_cost 0.14
	word WindowType Gaussian1
	real WindowWidth 1.0
	sentence OutFile /Users/frkkan96/Desktop/aaa.csv
endform

stopwatch

currSound = Read from file: soundFile$
soundEnd = Get end time

startAt  = startTime + selectionOffset

endAt = endTime
selEnd = startAt + selectionLength
if selEnd < endTime 
	endAt = selEnd
endif



select currSound
currSound = Extract part: startAt, endAt, windowType$, windowWidth, 0
noprogress To PointProcess (periodic, cc)... 'minimum_f0' 'maximum_f0'
currPP = selected("PointProcess")
select currSound
noprogress To Pitch (cc)... 0.0 'minimum_f0' 15 1 'silence_threshold' 'voicing_threshold' 'octave_cost' 'octave_jump_cost' 'voiced/unvoiced_cost' 'maximum_f0'
currPitch = selected ("Pitch")
plus currPP
plus currSound

voiceReport$ = Voice report... 0.0 0.0 'minimum_f0' 'maximum_f0' 'maximum_period_factor' 'maximum_amplitude_factor' 'silence_threshold' 'voicing_threshold'
#Pitch
medianPitch = extractNumber (voiceReport$, "Median pitch: ")
meanPitch = extractNumber (voiceReport$, "Mean pitch: ")
sdPitch = extractNumber (voiceReport$, "Standard deviation: ")
minPitch = extractNumber (voiceReport$, "Minimum pitch: ")
maxPitch = extractNumber (voiceReport$, "Maximum pitch: ")

#Pulses
numPulse = extractNumber (voiceReport$, "Number of pulses: ")
numPer = extractNumber (voiceReport$, "Number of periods: ")
meanPer = extractNumber (voiceReport$, "Mean period: ")
sdPer = extractNumber (voiceReport$, "Standard deviation of period: ")

#Voicing
fracVoice = extractNumber (voiceReport$, "Fraction of locally unvoiced frames: ")
breaksVoice = extractNumber (voiceReport$, "Number of voice breaks: ")
breakeratioVoice = extractNumber (voiceReport$, "Degree of voice breaks: ")

#Jitter
jitterLocal = extractNumber (voiceReport$, "Jitter (local): ")
jitterLocalAbs = extractNumber (voiceReport$, "Jitter (local, absolute): ")
jitterRap = extractNumber (voiceReport$, "Jitter (rap): ")
jitterPpq5 = extractNumber (voiceReport$, "Jitter (ppq5): ")
jitterDdp = extractNumber (voiceReport$, "Jitter (ddp): ")

#Shimmer
shimmerLocal = extractNumber (voiceReport$, "Shimmer (local): ")
shimmerLocalAbs = extractNumber (voiceReport$, "Shimmer (local, dB): ")
shimmerApq3 = extractNumber (voiceReport$, "Shimmer (apq3): ")
shimmerApq5 = extractNumber (voiceReport$, "Shimmer (apq5): ")
shimmerApq11 = extractNumber (voiceReport$, "Shimmer (apq11): ")
shimmerDda = extractNumber (voiceReport$, "Shimmer (dda): ")

#Harmonicity
meanAutocor = extractNumber (voiceReport$, "Mean autocorrelation: ")
meanNHR = extractNumber (voiceReport$, "Mean noise-to-harmonics ratio: ")
meanHNR = extractNumber (voiceReport$, "Mean harmonics-to-noise ratio: ")

#Intensity things

minPitch = meanPitch
if meanPitch = undefined
	minPitch = 20
endif
select currSound
noprogress To Intensity... 'minPitch' 0 0
intMean = Get mean... 0 0 energy
intMedian = Get quantile... 0 0 0.50
intSD = Get standard deviation... 0 0


outTab = Create Table with column names: "outTable", 1, { "Start Time","End Time","Selection start","Selection end","Median pitch", "Mean pitch", "Standard deviation", "Minimum pitch", "Maximum pitch", "Number of pulses", "Number of periods", "Mean period", "Standard deviation of period", "Fraction of locally unvoiced frames", "Number of voice breaks", "Degree of voice breaks", "Jitter (local)", "Jitter (local, absolute)", "Jitter (rap)", "Jitter (ppq5)", "Jitter (ddp)", "Shimmer (local)", "Shimmer (local, dB)", "Shimmer (apq3)", "Shimmer (apq5)", "Shimmer (apq11)", "Shimmer (dda)", "Mean autocorrelation", "Mean noise-to-harmonics ratio", "Mean harmonics-to-noise ratio" }

Set numeric value: 1, "Start Time" , 'startTime'
Set numeric value: 1, "Selection start", 'startAt'

if selectionLength == 0.0
	Set numeric value: 1, "End Time", 'soundEnd'
	Set numeric value: 1, "Selection end", 'soundEnd'
else 
	Set numeric value: 1, "End Time", 'endTime'
	Set numeric value: 1, "Selection end", 'endAt'
endif


if medianPitch <> undefined 
	Set numeric value: 1, "Median pitch", 'medianPitch'
endif
if meanPitch <> undefined 
	Set numeric value: 1, "Mean pitch", 'meanPitch'
endif
if sdPitch <> undefined 
	Set numeric value: 1, "Standard deviation", 'sdPitch'
endif
if minPitch <> undefined 
	Set numeric value: 1, "Minimum pitch", 'minPitch'
endif
if maxPitch <> undefined 
	Set numeric value: 1, "Maximum pitch", 'maxPitch'
endif
if numPulse <> undefined 
	Set numeric value: 1, "Number of pulses", 'numPulse'
endif
if numPer <> undefined 
	Set numeric value: 1, "Number of periods", 'numPer'
endif
if meanPer <> undefined 
	Set numeric value: 1, "Mean period", 'meanPer'
endif
if sdPer <> undefined 
	Set numeric value: 1, "Standard deviation of period", 'sdPer'
endif
if fracVoice <> undefined 
	Set numeric value: 1, "Fraction of locally unvoiced frames", 'fracVoice'
endif
if breaksVoice <> undefined 
	Set numeric value: 1, "Number of voice breaks", 'breaksVoice'
endif
if breakeratioVoice <> undefined 
	Set numeric value: 1, "Degree of voice breaks", 'breakeratioVoice'
endif
if jitterLocal <> undefined 
	Set numeric value: 1, "Jitter (local)", 'jitterLocal'
endif
if jitterLocalAbs <> undefined 
	Set numeric value: 1, "Jitter (local, absolute)", 'jitterLocalAbs'
endif
if jitterRap <> undefined 
	Set numeric value: 1, "Jitter (rap)", 'jitterRap'
endif
if jitterPpq5 <> undefined 
	Set numeric value: 1, "Jitter (ppq5)", 'jitterPpq5'
endif
if jitterDdp <> undefined 
	Set numeric value: 1, "Jitter (ddp)", 'jitterDdp'
endif
if shimmerLocal <> undefined 
	Set numeric value: 1, "Shimmer (local)", 'shimmerLocal'
endif
if shimmerLocalAbs <> undefined 
	Set numeric value: 1, "Shimmer (local, dB)", 'shimmerLocalAbs'
endif
if shimmerApq3 <> undefined 
	Set numeric value: 1, "Shimmer (apq3)", 'shimmerApq3'
endif
if shimmerApq5 <> undefined 
	Set numeric value: 1, "Shimmer (apq5)", 'shimmerApq5'
endif
if shimmerApq11 <> undefined 
	Set numeric value: 1, "Shimmer (apq11)", 'shimmerApq11'
endif
if shimmerDda <> undefined 
	Set numeric value: 1, "Shimmer (dda)", 'shimmerDda'
endif
if meanAutocor <> undefined 
	Set numeric value: 1, "Mean autocorrelation", 'meanAutocor'
endif
if meanNHR <> undefined 
	Set numeric value: 1, "Mean noise-to-harmonics ratio", 'meanNHR'
endif
if meanHNR <> undefined 
	Set numeric value: 1, "Mean harmonics-to-noise ratio", 'meanHNR'
endif

Save as semicolon-separated file: "'outFile$'"
