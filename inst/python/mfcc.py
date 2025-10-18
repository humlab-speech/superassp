#!/usr/bin/env python3
"""
MFCC feature extraction using torchaudio.

This script extracts Mel-Frequency Cepstral Coefficients (MFCCs) from audio
using torchaudio. MFCCs are widely used in speech recognition, speaker 
identification, and audio analysis.
"""

import torch
import torchaudio
import torchaudio.transforms as T
import numpy as np
import math
import sys
import json


def mfcc(
    soundFile,
    windowShift=10.0,
    windowSize=25.0,
    n_mfcc=13,
    n_mels=40,
    fmin=0.0,
    fmax=None,
    beginTime=0.0,
    endTime=0.0
):
    """
    Extract MFCC features from audio.
    
    Parameters
    ----------
    soundFile : str
        Path to audio file
    windowShift : float
        Frame shift in milliseconds (default: 10.0)
    windowSize : float
        FFT window size in milliseconds (default: 25.0)
    n_mfcc : int
        Number of MFCC coefficients (default: 13)
    n_mels : int
        Number of mel filterbanks (default: 40)
    fmin : float
        Minimum frequency in Hz (default: 0.0)
    fmax : float
        Maximum frequency in Hz (default: None = sample_rate/2)
    beginTime : float
        Start time in seconds (default: 0.0)
    endTime : float
        End time in seconds (default: 0.0 = end of file)
    
    Returns
    -------
    dict
        Dictionary with keys:
        - mfcc: numpy array of shape [n_frames, n_mfcc]
        - sample_rate: original sample rate
        - n_frames: number of frames
        - n_mfcc: number of coefficients
    """
    # Get audio metadata
    metadata = torchaudio.info(soundFile)
    sample_rate = metadata.sample_rate
    
    # Calculate frame offset and duration
    if beginTime > 0 and endTime > beginTime:
        start_sample = math.floor(beginTime * sample_rate)
    else:
        start_sample = 0
    
    if endTime > 0 and endTime > beginTime:
        end_sample = math.ceil(endTime * sample_rate)
        n_samples = end_sample - start_sample
    else:
        n_samples = -1
    
    # Load audio
    waveform, _ = torchaudio.load(soundFile, frame_offset=start_sample, num_frames=n_samples)
    
    # Calculate FFT parameters
    hop_length = int(sample_rate * windowShift / 1000)
    n_fft = int(sample_rate * windowSize / 1000)
    # Round to next power of 2 for efficiency
    n_fft = 2 ** math.ceil(math.log2(n_fft))
    
    # Set fmax if not specified
    if fmax is None or fmax <= 0:
        fmax = sample_rate / 2.0
    
    # Create MFCC transform
    mfcc_transform = T.MFCC(
        sample_rate=sample_rate,
        n_mfcc=n_mfcc,
        melkwargs={
            'n_fft': n_fft,
            'hop_length': hop_length,
            'n_mels': n_mels,
            'f_min': fmin,
            'f_max': fmax
        }
    )
    
    # Compute MFCCs
    mfcc_features = mfcc_transform(waveform)
    
    # Convert to numpy (shape: [channel, n_mfcc, time])
    mfcc_np = mfcc_features.numpy()
    
    # Average over channels if stereo
    if mfcc_np.shape[0] > 1:
        mfcc_np = mfcc_np.mean(axis=0, keepdims=True)
    
    # Transpose to [time, n_mfcc] for easier handling in R
    mfcc_np = mfcc_np[0].T
    
    return {
        'mfcc': mfcc_np.tolist(),
        'sample_rate': float(sample_rate),
        'n_frames': mfcc_np.shape[0],
        'n_mfcc': mfcc_np.shape[1]
    }


if __name__ == '__main__':
    # When called from R, parameters are passed via command line arguments
    # Format: python mfcc.py <json_params>
    if len(sys.argv) > 1:
        params = json.loads(sys.argv[1])
        result = mfcc(**params)
        print(json.dumps(result))
    else:
        print("Error: No parameters provided", file=sys.stderr)
        sys.exit(1)
