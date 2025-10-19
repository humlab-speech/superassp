#!/usr/bin/env python3
"""
Example: Basic Voice Analysis

Demonstrates basic usage of the Voice Analysis Toolbox
"""

import sys
sys.path.insert(0, '..')

from voice_analysis import analyze_voice_file
import numpy as np


def main():
    # Analyze a voice recording
    print("Voice Analysis Toolbox - Example")
    print("=" * 60)
    
    # Path to audio file (adjust as needed)
    audio_file = "../../a1.wav"
    
    print(f"\nAnalyzing: {audio_file}")
    print("-" * 60)
    
    # Run analysis with SWIPE F0 algorithm
    measures, F0 = analyze_voice_file(
        audio_file,
        f0_min=50,
        f0_max=500,
        f0_algorithm='SWIPE'
    )
    
    print("\n" + "=" * 60)
    print(f"Results: {len(measures)} measures computed")
    print("=" * 60)
    
    # Display F0 statistics
    F0_valid = F0[F0 > 0]
    print(f"\nF0 Statistics:")
    print(f"  Mean: {np.mean(F0_valid):.2f} Hz")
    print(f"  Std:  {np.std(F0_valid):.2f} Hz")
    print(f"  Range: {np.min(F0_valid):.2f} - {np.max(F0_valid):.2f} Hz")
    
    # Display some key measures
    print(f"\nKey Measures:")
    print(f"  Jitter RAP:     {measures.get('jitter_RAP', np.nan):.6f}")
    print(f"  Shimmer RAP:    {measures.get('shimmer_RAP', np.nan):.6f}")
    print(f"  HNR mean:       {measures.get('HNR_mean', np.nan):.2f} dB")
    print(f"  NHR mean:       {measures.get('NHR_mean', np.nan):.6f}")
    print(f"  PPE:            {measures.get('PPE', np.nan):.6f}")
    print(f"  DFA:            {measures.get('DFA', np.nan):.6f}")
    print(f"  RPDE:           {measures.get('RPDE', np.nan):.6f}")
    
    # Count measures by category
    categories = {
        'Jitter': 0,
        'Shimmer': 0,
        'HNR/NHR': 0,
        'MFCC': 0,
        'Wavelet': 0,
        'GNE': 0,
        'Other': 0
    }
    
    for key in measures.keys():
        if key.startswith('jitter_'):
            categories['Jitter'] += 1
        elif key.startswith('shimmer_'):
            categories['Shimmer'] += 1
        elif key.startswith(('HNR', 'NHR')):
            categories['HNR/NHR'] += 1
        elif key.startswith('MFCC'):
            categories['MFCC'] += 1
        elif key.startswith('wavelet'):
            categories['Wavelet'] += 1
        elif key.startswith('GNE'):
            categories['GNE'] += 1
        else:
            categories['Other'] += 1
    
    print(f"\nMeasures by Category:")
    for cat, count in categories.items():
        print(f"  {cat:12s}: {count:3d}")
    
    print("\n" + "=" * 60)
    print("Analysis complete!")
    print("=" * 60)


if __name__ == '__main__':
    main()
