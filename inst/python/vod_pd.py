#!/Users/frkkan96/opt/anaconda3/envs/pyannote/bin/python

import sys
import pandas as pd

wavfile = sys.argv[1]

# instantiate pretrained speaker diarization pipeline
from pyannote.audio import Pipeline
pipeline = Pipeline.from_pretrained("pyannote/voice-activity-detection")

# apply pretrained pipeline
diarization = pipeline(wavfile)


# print the result
for turn, _, speaker in diarization.itertracks(yield_label=True):
    print(f"{turn.start:.4f};{turn.end:.4f};{wavfile}")
