

form Ange attribut
	sentence SoundFile /Users/frkkan96/Desktop/aaa.wav
	real StartTime 0.0
	real EndTime 0.0
	real SelectFrom 0.0
	real SelectionLength 2.0
	word WindowType Hanning
	real WindowWidth 1.0
endform


currSound = Read from file: soundFile$
soundEnd = Get end time

endAt = endTime
selEnd = selectFrom + selectionLength
if selEnd < endTime 
	endAt = selEnd
endif

startAt  = max(startTime, selectFrom)

select currSound
currSound = Extract part: startAt, endAt, windowType$, windowWidth, 0
noprogress To PointProcess (periodic, cc)... 30 600
currPP = selected("PointProcess")
select currSound
noprogress To Pitch (cc)... 0.0 30 15 1 0.03 0.45 0.01 0.35 0.14 600
currPitch = selected ("Pitch")
plus currPP
plus currSound

voiceReport$ = Voice report... 0.0 0.0 30 600 1.3 1.6 0.03 0.45
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
shimmerLocalAbs = extractNumber (voiceReport$, "Shimmer (local, absolute): ")
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



#Pitch
out$ = "'startAt';'endAt';'startTime';'endTime';'medianPitch';'meanPitch';'sdPitch';'minPitch';'maxPitch'"
header$ = "Selection start;Selection end;Vowel start;Vowel end;Median Pitch; Mean Pitch; Pitch SD;Min Pitch;Max Pitch"
#Pulse
out$ = out$ + ";'numPulse';'numPer';'meanPer';'sdPer'"
header$ = header$ + ";No Pulses;No Periods;Mean period;Period SD"
#Voice 
out$ = out$ + ";'fracVoice';'breaksVoice';'breakeratioVoice'"
header$ = header$ + ";Frac local unvoiced frames;Voice breaks;Degree voice breaks"
#Jitter
out$ = out$ + ";'jitterLocal';'jitterLocalAbs';'jitterRap';'jitterPpq5';'jitterDdp'"
header$ = header$ + ";Jitter (local);Jitter (local, absolute);Jitter (rap);Jitter (ppq5);Jitter (ddp)"
#Shimmer
out$ = out$ + ";'shimmerLocal';'shimmerLocalAbs';'shimmerApq3';'shimmerApq5';'shimmerApq11';'shimmerDda'"
header$ = header$ + ";Shimmer (local);S‘himmer (local, absolute);Shimmer (apq3);Shimmer (apq5);Shimmer (apq11);Shimmer (dda)"
#Harmonicity
out$ = out$ + ";'meanAutocor';'meanNHR';'meanHNR'"
header$ = header$ + ";Mean Autocorrelation;Mean noise-to-harmonics ratio;Mean harmonics-to-noise ratio"
#Intensity
out$ = out$ + ";'intMean';'intMedian';'intSD'"
header$ =header$ + ";Mean intensity;Median intensity;Intensity standard deviation"

writeInfoLine: "'header$'\n"
appendInfoLine: "'out$'", "\n"
