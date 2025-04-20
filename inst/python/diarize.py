#!/Users/frkkan96/opt/anaconda3/envs/pyannote/bin/python

import sys
import pandas as pd
import os

wavfile = sys.argv[1]

# instantiate pretrained speaker diarization pipeline
from pyannote.audio import Pipeline


def diarize(wavfile, session, level, attribute, YOUR_AUTH_TOKEN):

    pipeline = Pipeline.from_pretrained("pyannote/speaker-diarization-3.1", use_auth_token=YOUR_AUTH_TOKEN)

    # apply pretrained pipeline
    diarization = pipeline(wavfile)


    diadf = pd.DataFrame(
        data={
            "start":[turn.start for turn,_,speaker in diarization.itertracks(yield_label=True)],
            "end":[turn.end for turn,_,speaker in diarization.itertracks(yield_label=True)],
            "labels":[speaker for turn,_,speaker in diarization.itertracks(yield_label=True)]
            }
        )
    diadf["bundle"] = wavfile.removesuffix(".wav")
    diadf["session"] = session
    diadf["level"] = level
    diadf["attribute"] = attribute
    diadf['start_item_seq_idx'] = diadf.index

    return diadf


