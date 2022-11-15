#!/Users/frkkan96/opt/anaconda3/envs/pyannote/bin/python

import sys
import pandas as pd
import os

wavfile = sys.argv[1]

# instantiate pretrained speaker diarization pipeline
from pyannote.audio import Pipeline
pipeline = Pipeline.from_pretrained("pyannote/speaker-diarization")

arr = []

# apply pretrained pipeline
diarization = pipeline(wavfile)

# print the result
for turn, _, speaker in diarization.itertracks(yield_label=True):
    #print(f"start={turn.start:.1f}s stop={turn.end:.1f}s speaker_{speaker}")
    arr.append( {"start":turn.start, "end": turn.end, "speaker":speaker, "file":wavfile} )

outpd = pd.DataFrame(arr)


csvFile = os.path.splitext(wavfile)[0]+'_dia.csv'


outpd.to_csv(csvFile,sep=",", columns=['start','end','speaker','file'],encoding='utf-8',index=False)