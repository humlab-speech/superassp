import numpy as np
import pysptk as sp
#from scipy.io import wavfile
import librosa as lr
x, fs = lr.load("/Users/frkkan96/Desktop/a1.wav",dtype=np.float64, offset= 0.5, duration = 15.0)

#fs, x = wavfile.read(pysptk.util.example_audio_file())
#x, f = lr.lood(pysptk.util.example_audio_file())
f0_swipe = sp.swipe(x, fs=fs, hopsize=80, min=60, max=200, otype="f0")

