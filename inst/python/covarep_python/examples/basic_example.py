"""
Example: Basic F0 Tracking and Glottal Analysis

This example demonstrates:
1. F0 estimation using SRH method
2. Glottal inverse filtering using IAIF
3. Basic voice quality parameter extraction
"""

import numpy as np
import matplotlib.pyplot as plt
import soundfile as sf
from covarep.f0 import F0Tracker
from covarep.glottal import iaif, get_vq_params


def example_f0_tracking():
    """
    Example of F0 tracking using SRH method
    """
    # Load test audio (you'll need to provide a wav file)
    audio_file = "test_speech.wav"
    
    try:
        audio, fs = sf.read(audio_file)
    except FileNotFoundError:
        print(f"Please provide a test audio file: {audio_file}")
        # Create synthetic signal for demonstration
        t = np.linspace(0, 1, 16000)
        f0 = 150  # Hz
        audio = np.sin(2 * np.pi * f0 * t)
        fs = 16000
        print("Using synthetic signal for demonstration")
    
    # Create F0 tracker
    tracker = F0Tracker(method='srh', f0_min=50, f0_max=400)
    
    # Estimate F0
    f0, vuv, srh_values, times = tracker.estimate(audio, fs)
    
    # Plot results
    fig, axes = plt.subplots(3, 1, figsize=(12, 8))
    
    # Plot waveform
    t_audio = np.arange(len(audio)) / fs
    axes[0].plot(t_audio, audio, 'k', linewidth=0.5)
    axes[0].set_ylabel('Amplitude')
    axes[0].set_title('Speech Signal')
    axes[0].grid(True, alpha=0.3)
    
    # Plot F0
    axes[1].plot(times, f0, 'b.-', linewidth=1.5, markersize=3)
    axes[1].fill_between(times, 0, 500, where=~vuv, alpha=0.2, color='red', label='Unvoiced')
    axes[1].set_ylabel('F0 (Hz)')
    axes[1].set_title('F0 Contour')
    axes[1].set_ylim([0, 500])
    axes[1].grid(True, alpha=0.3)
    axes[1].legend()
    
    # Plot SRH values
    axes[2].plot(times, srh_values, 'g.-', linewidth=1.5, markersize=3)
    axes[2].set_xlabel('Time (s)')
    axes[2].set_ylabel('SRH Value')
    axes[2].set_title('SRH Criterion')
    axes[2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('f0_tracking_example.png', dpi=150)
    print("F0 tracking example saved to f0_tracking_example.png")
    
    return f0, vuv, times


def example_glottal_analysis():
    """
    Example of glottal inverse filtering using IAIF
    """
    # Load or create test signal
    try:
        audio, fs = sf.read("test_speech.wav")
    except FileNotFoundError:
        # Create synthetic signal
        t = np.linspace(0, 1, 16000)
        f0 = 150
        audio = np.sin(2 * np.pi * f0 * t)
        fs = 16000
        print("Using synthetic signal for demonstration")
    
    # Take a voiced segment (first 30 ms)
    frame_len = int(0.03 * fs)
    if len(audio) > frame_len:
        frame = audio[:frame_len]
    else:
        frame = audio
    
    # Apply IAIF
    g, dg, a, ag = iaif(frame, fs)
    
    # Get voice quality parameters
    vq_params = get_vq_params(frame, fs)
    
    # Plot results
    fig, axes = plt.subplots(3, 1, figsize=(12, 8))
    
    t = np.arange(len(frame)) / fs * 1000  # in ms
    
    # Original signal
    axes[0].plot(t, frame, 'k', linewidth=1)
    axes[0].set_ylabel('Amplitude')
    axes[0].set_title('Original Speech Signal')
    axes[0].grid(True, alpha=0.3)
    
    # Glottal flow
    t_g = np.arange(len(g)) / fs * 1000
    axes[1].plot(t_g, g, 'b', linewidth=1)
    axes[1].set_ylabel('Amplitude')
    axes[1].set_title('Estimated Glottal Flow')
    axes[1].grid(True, alpha=0.3)
    
    # Glottal flow derivative
    axes[2].plot(t_g, dg, 'r', linewidth=1)
    axes[2].set_xlabel('Time (ms)')
    axes[2].set_ylabel('Amplitude')
    axes[2].set_title('Glottal Flow Derivative')
    axes[2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('iaif_example.png', dpi=150)
    print("IAIF example saved to iaif_example.png")
    
    # Print voice quality parameters
    print("\nVoice Quality Parameters:")
    for key, value in vq_params.items():
        if not np.isnan(value):
            print(f"  {key}: {value:.4f}")
        else:
            print(f"  {key}: Not computed (requires additional processing)")
    
    return g, dg, vq_params


if __name__ == "__main__":
    print("=" * 60)
    print("COVAREP Python - Basic Examples")
    print("=" * 60)
    
    print("\n1. F0 Tracking Example")
    print("-" * 60)
    f0, vuv, times = example_f0_tracking()
    print(f"   F0 range: {np.min(f0[vuv]):.1f} - {np.max(f0[vuv]):.1f} Hz")
    print(f"   Voicing: {100*np.sum(vuv)/len(vuv):.1f}%")
    
    print("\n2. Glottal Analysis Example")
    print("-" * 60)
    g, dg, vq_params = example_glottal_analysis()
    
    print("\n" + "=" * 60)
    print("Examples completed successfully!")
    print("=" * 60)
