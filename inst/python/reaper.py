import numpy as np
#from scipy.io import wavfile
import librosa as lr
import pyreaper
x, fs = lr.load("/Users/frkkan96/Desktop/a1.wav",dtype=np.float64, offset= 0.5, duration = 15.0)

max_16bit = 2**15
raw_x = x * max_16bit

# now change the data type
raw_x = raw_x.astype(np.int16)

#fs, x = wavfile.read(pysptk.util.example_audio_file())
#x, f = lr.lood(pysptk.util.example_audio_file())
#f0_swipe = sp.swipe(x, fs=fs, hopsize=80, min=60, max=200, otype="f0")
pm_times, pm, f0_times, f0, corr = pyreaper.reaper(x=raw_x, fs=fs, minf0 = 80, maxf0 = 200, do_high_pass=True, do_hilbert_transform= False, inter_pulse =0.01, frame_period=0.005, unvoiced_cost =0.9)

print(f0_times)

