#!/usr/bin/env python3
"""
Demonstration script for pure Python prosody_measures implementation

This script shows how to use the optimized dysprosody_pure.prosody_measures
function to analyze audio files without requiring compiled C binaries or Perl.

Usage:
    python demo_pure.py <audio_file.wav>
    python demo_pure.py  # analyzes all .wav files in current directory
"""

import sys
import glob
import os
import pandas as pd
from pathlib import Path
from dysprosody_pure import prosody_measures


def analyze_single_file(audio_file):
    """Analyze a single audio file and display results"""
    print(f"\n{'='*70}")
    print(f"Analyzing: {audio_file}")
    print('='*70)

    try:
        result = prosody_measures(audio_file)

        if result is None:
            print("⚠️  File skipped (duration < 1 second)")
            return None

        print(f"\n✓ Successfully extracted {len(result)} features")

        # Display key metadata
        print("\n📊 Key Metadata:")
        metadata_keys = [
            'Duration', 'PitchKey', 'PitchRange', 'PitchMean',
            'IntsIntLabels', 'UniqueIntsInt', 'IntsIntConcentration'
        ]
        for key in metadata_keys:
            if key in result.index:
                print(f"   {key:25s} = {result[key]:.2f}")

        # Display sample prosodic features
        print("\n📈 Sample Prosodic Features:")
        feature_samples = [k for k in result.index if any(x in k for x in ['L2L1', 'SLF', 'C1', 'Spectral'])][:5]
        for key in feature_samples:
            print(f"   {key:25s} = {result[key]:.4f}")

        return result

    except Exception as e:
        print(f"❌ Error analyzing {audio_file}: {e}")
        import traceback
        traceback.print_exc()
        return None


def analyze_directory(directory=".", pattern="*.wav"):
    """Analyze all WAV files in a directory and return DataFrame"""
    wav_files = glob.glob(os.path.join(directory, pattern))

    if not wav_files:
        print(f"No .wav files found in {directory}")
        return None

    print(f"\n🔍 Found {len(wav_files)} WAV file(s)")

    results = {}
    for wav_file in wav_files:
        result = analyze_single_file(wav_file)
        if result is not None:
            basename = os.path.basename(wav_file).removesuffix('.wav')
            results[basename] = result

    if results:
        df = pd.DataFrame(results).T
        print(f"\n✓ Created DataFrame with {len(df)} files × {len(df.columns)} features")
        return df
    else:
        print("\n⚠️  No files were successfully analyzed")
        return None


def main():
    """Main entry point"""
    print("=" * 70)
    print(" Pure Python Prosody Analysis (dysprosody_pure)")
    print("=" * 70)
    print("\nNo compiled binaries or Perl required! 🎉\n")

    if len(sys.argv) > 1:
        # Analyze specific file(s)
        for audio_file in sys.argv[1:]:
            if not os.path.exists(audio_file):
                print(f"❌ File not found: {audio_file}")
                continue

            result = analyze_single_file(audio_file)

            if result is not None:
                # Optionally save to CSV
                output_file = audio_file.replace('.wav', '_prosody.csv')
                result.to_csv(output_file)
                print(f"\n💾 Saved results to: {output_file}")

    else:
        # Batch analysis of all WAV files in current directory
        df = analyze_directory()

        if df is not None:
            # Save batch results
            output_file = "prosody_batch_results.csv"
            df.to_csv(output_file)
            print(f"\n💾 Saved batch results to: {output_file}")

            # Display summary statistics
            print("\n📊 Summary Statistics:")
            summary_cols = ['Duration', 'PitchMean', 'IntsIntLabels', 'UniqueIntsInt']
            available_cols = [c for c in summary_cols if c in df.columns]
            if available_cols:
                print(df[available_cols].describe())


if __name__ == "__main__":
    main()
