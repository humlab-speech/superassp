soundFile = "/Users/frkkan96/Desktop/a1.wav"
minF = 70.0
maxF =400.0
voiced_voiceless_threshold = 0.1
windowShift = 10.0
beginTime = 1.2
endTime = 1.6
dimensions = 1

duration = endTime - beginTime

if duration < (windowShift / 1000) :
	duration = None


import pyworld as pw
import librosa as lr
import numpy as np



x, fs = lr.load(soundFile,
	dtype=np.float64,
	offset= beginTime,
	duration= duration
	)

f0, t = pw.harvest(x,
	fs,
	f0_floor=minF,
	f0_ceil=maxF,
	frame_period=windowShift )


sp = pw.cheaptrick(x,f0, t,fs, f0_floor=minF)

sl = pw.code_spectral_envelope(sp,fs,dimensions)