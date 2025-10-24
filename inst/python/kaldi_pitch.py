#!/usr/bin/env python3
"""
Kaldi-style pitch tracking using torchaudio.

This script provides pitch tracking using normalized cross-correlation function (NCCF) 
and median smoothing, similar to the Kaldi ASR toolkit's pitch tracker.

Uses torchaudio.functional.detect_pitch_frequency which replaced the deprecated
compute_kaldi_pitch (removed in torchaudio 2.9+).
"""

import torch
import torchaudio
import torchaudio.functional as F
import numpy as np
import math
import sys
import json


def kaldi_pitch(
    waveform=None,
    sample_rate=None,
    soundFile=None,
    windowShift=10.0,
    windowSize=30,
    minF=85.0,
    maxF=400.0,
    beginTime=0.0,
    endTime=0.0
):
    """
    Extract F0 using Kaldi-style pitch tracker.

    Parameters
    ----------
    waveform : torch.Tensor, optional
        Pre-loaded audio waveform (preferred method)
    sample_rate : int, optional
        Sample rate of the waveform (required if waveform is provided)
    soundFile : str, optional
        Path to audio file (legacy method, for backward compatibility)
    windowShift : float
        Frame shift in milliseconds (default: 10.0)
    windowSize : int
        Window length for median smoothing in number of frames (default: 30)
    minF : float
        Minimum F0 in Hz (default: 85.0)
    maxF : float
        Maximum F0 in Hz (default: 400.0)
    beginTime : float
        Start time in seconds (ignored if waveform is provided, default: 0.0)
    endTime : float
        End time in seconds (ignored if waveform is provided, default: 0.0)

    Returns
    -------
    dict
        Dictionary with keys:
        - f0: numpy array of F0 values
        - sample_rate: original sample rate
        - n_frames: number of frames
    """
    # Check if detect_pitch_frequency is available
    if not hasattr(F, 'detect_pitch_frequency'):
        raise AttributeError(
            'torchaudio.functional.detect_pitch_frequency not found. '
            'Please install torchaudio >= 0.13.0'
        )

    # Load audio if not provided
    if waveform is None:
        if soundFile is None:
            raise ValueError("Either waveform or soundFile must be provided")

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
    else:
        # Waveform provided - ensure it's a tensor
        if not isinstance(waveform, torch.Tensor):
            waveform = torch.from_numpy(waveform)

        # Ensure correct shape (add channel dimension if needed)
        if waveform.ndim == 1:
            waveform = waveform.unsqueeze(0)

        if sample_rate is None:
            raise ValueError("sample_rate must be provided when waveform is given")
    
    # Detect pitch
    frame_time = windowShift / 1000.0  # Convert ms to seconds
    pitch = F.detect_pitch_frequency(
        waveform=waveform,
        sample_rate=sample_rate,
        frame_time=frame_time,
        win_length=int(windowSize),
        freq_low=int(minF),
        freq_high=int(maxF)
    )
    
    # Convert to numpy
    f0_values = pitch.numpy()
    
    # Average over channels if stereo
    if f0_values.ndim > 1 and f0_values.shape[0] > 1:
        f0_values = f0_values.mean(axis=0)
    elif f0_values.ndim > 1:
        f0_values = f0_values[0]
    
    return {
        'f0': f0_values.tolist(),
        'sample_rate': float(sample_rate),
        'n_frames': len(f0_values)
    }


if __name__ == '__main__':
    # When called from R, parameters are passed via command line arguments
    # Format: python kaldi_pitch.py <json_params>
    if len(sys.argv) > 1:
        params = json.loads(sys.argv[1])
        result = kaldi_pitch(**params)
        print(json.dumps(result))
    else:
        print("Error: No parameters provided", file=sys.stderr)
        sys.exit(1)
