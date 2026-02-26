"""
Test Optimizations with Multiple Audio Files

Tests the optimized implementations across diverse audio samples:
- Different speakers
- Different recording conditions
- Different durations
- Different sample rates

Generates comprehensive validation report.
"""

import sys
import os
import time
import numpy as np
import soundfile as sf
from pathlib import Path
import glob

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from covarep.f0 import pitch_srh
from covarep.f0.f0_optimized import pitch_srh_vectorized
from covarep.glottal import iaif
from covarep.glottal.iaif_optimized import iaif_optimized


def test_audio_file(audio_file, test_f0=True, test_iaif=True):
    """
    Test both F0 and IAIF on a single audio file

    Returns
    -------
    results : dict
        Test results
    """
    results = {
        'file': os.path.basename(audio_file),
        'success': True,
        'error': None
    }

    try:
        # Load audio
        audio, fs = sf.read(audio_file)
        audio = audio[:int(5 * fs)]  # Limit to 5 seconds for speed

        results['duration'] = len(audio) / fs
        results['fs'] = fs
        results['samples'] = len(audio)

        # Test F0 tracking
        if test_f0:
            # Original
            start = time.perf_counter()
            f0_orig, vuv_orig, srh_orig, times_orig = pitch_srh(audio, fs)
            time_orig = time.perf_counter() - start

            # Optimized
            start = time.perf_counter()
            f0_opt, vuv_opt, srh_opt, times_opt = pitch_srh_vectorized(audio, fs)
            time_opt = time.perf_counter() - start

            # Validate
            both_voiced = vuv_orig & vuv_opt
            if np.sum(both_voiced) > 0:
                f0_diff = np.abs(f0_orig[both_voiced] - f0_opt[both_voiced])
                results['f0_median_error'] = np.median(f0_diff)
                results['f0_mean_error'] = np.mean(f0_diff)
                results['f0_max_error'] = np.max(f0_diff)
            else:
                results['f0_median_error'] = np.nan
                results['f0_mean_error'] = np.nan
                results['f0_max_error'] = np.nan

            results['f0_time_orig'] = time_orig * 1000  # ms
            results['f0_time_opt'] = time_opt * 1000
            results['f0_speedup'] = time_orig / time_opt if time_opt > 0 else np.inf
            results['f0_vuv_agreement'] = 100 * np.sum(vuv_orig == vuv_opt) / len(vuv_orig)

        # Test IAIF (on first 30ms)
        if test_iaif and len(audio) >= int(0.03 * fs):
            frame = audio[:int(0.03 * fs)]

            # Original
            start = time.perf_counter()
            g_orig, dg_orig, a_orig, ag_orig = iaif(frame, fs)
            time_orig = time.perf_counter() - start

            # Optimized
            start = time.perf_counter()
            g_opt, dg_opt, a_opt, ag_opt = iaif_optimized(frame, fs)
            time_opt = time.perf_counter() - start

            # Validate
            if len(g_orig) == len(g_opt):
                g_diff = np.max(np.abs(g_orig - g_opt))
                dg_diff = np.max(np.abs(dg_orig - dg_opt))
                results['iaif_g_max_diff'] = g_diff
                results['iaif_dg_max_diff'] = dg_diff
            else:
                results['iaif_g_max_diff'] = np.nan
                results['iaif_dg_max_diff'] = np.nan

            results['iaif_time_orig'] = time_orig * 1000  # ms
            results['iaif_time_opt'] = time_opt * 1000
            results['iaif_speedup'] = time_orig / time_opt if time_opt > 0 else np.inf

    except Exception as e:
        results['success'] = False
        results['error'] = str(e)

    return results


def generate_report(all_results, output_file='multi_file_test_report.txt'):
    """
    Generate comprehensive test report
    """
    with open(output_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("COVAREP PYTHON - MULTI-FILE VALIDATION REPORT\n")
        f.write("=" * 80 + "\n\n")

        # Summary
        successful = sum(1 for r in all_results if r['success'])
        failed = len(all_results) - successful

        f.write(f"Total Files Tested: {len(all_results)}\n")
        f.write(f"Successful: {successful}\n")
        f.write(f"Failed: {failed}\n\n")

        if failed > 0:
            f.write("Failed Files:\n")
            for r in all_results:
                if not r['success']:
                    f.write(f"  - {r['file']}: {r['error']}\n")
            f.write("\n")

        # F0 Statistics
        f.write("-" * 80 + "\n")
        f.write("F0 TRACKING PERFORMANCE\n")
        f.write("-" * 80 + "\n\n")

        f0_speedups = [r['f0_speedup'] for r in all_results if r['success'] and 'f0_speedup' in r]
        f0_errors_median = [r['f0_median_error'] for r in all_results if r['success'] and not np.isnan(r.get('f0_median_error', np.nan))]
        f0_errors_mean = [r['f0_mean_error'] for r in all_results if r['success'] and not np.isnan(r.get('f0_mean_error', np.nan))]

        if f0_speedups:
            f.write(f"Speedup Statistics:\n")
            f.write(f"  Mean:   {np.mean(f0_speedups):.2f}x\n")
            f.write(f"  Median: {np.median(f0_speedups):.2f}x\n")
            f.write(f"  Min:    {np.min(f0_speedups):.2f}x\n")
            f.write(f"  Max:    {np.max(f0_speedups):.2f}x\n\n")

        if f0_errors_median:
            f.write(f"Accuracy (Median F0 Error):\n")
            f.write(f"  Mean:   {np.mean(f0_errors_median):.3f} Hz\n")
            f.write(f"  Median: {np.median(f0_errors_median):.3f} Hz\n")
            f.write(f"  Max:    {np.max(f0_errors_median):.3f} Hz\n\n")

        # IAIF Statistics
        f.write("-" * 80 + "\n")
        f.write("IAIF PERFORMANCE\n")
        f.write("-" * 80 + "\n\n")

        iaif_speedups = [r['iaif_speedup'] for r in all_results if r['success'] and 'iaif_speedup' in r]
        iaif_errors_g = [r['iaif_g_max_diff'] for r in all_results if r['success'] and not np.isnan(r.get('iaif_g_max_diff', np.nan))]

        if iaif_speedups:
            f.write(f"Speedup Statistics:\n")
            f.write(f"  Mean:   {np.mean(iaif_speedups):.2f}x\n")
            f.write(f"  Median: {np.median(iaif_speedups):.2f}x\n")
            f.write(f"  Min:    {np.min(iaif_speedups):.2f}x\n")
            f.write(f"  Max:    {np.max(iaif_speedups):.2f}x\n\n")

        if iaif_errors_g:
            f.write(f"Accuracy (Max Glottal Flow Diff):\n")
            f.write(f"  Mean:   {np.mean(iaif_errors_g):.6f}\n")
            f.write(f"  Median: {np.median(iaif_errors_g):.6f}\n")
            f.write(f"  Max:    {np.max(iaif_errors_g):.6f}\n\n")

        # Detailed Results
        f.write("-" * 80 + "\n")
        f.write("DETAILED RESULTS PER FILE\n")
        f.write("-" * 80 + "\n\n")

        for r in all_results:
            if not r['success']:
                continue

            f.write(f"File: {r['file']}\n")
            f.write(f"  Duration: {r['duration']:.2f}s @ {r['fs']} Hz\n")

            if 'f0_speedup' in r:
                f.write(f"  F0:\n")
                f.write(f"    Speedup: {r['f0_speedup']:.2f}x\n")
                f.write(f"    Time (orig): {r['f0_time_orig']:.1f} ms\n")
                f.write(f"    Time (opt):  {r['f0_time_opt']:.1f} ms\n")
                if not np.isnan(r.get('f0_median_error', np.nan)):
                    f.write(f"    Median error: {r['f0_median_error']:.3f} Hz\n")
                f.write(f"    VUV agreement: {r['f0_vuv_agreement']:.1f}%\n")

            if 'iaif_speedup' in r:
                f.write(f"  IAIF:\n")
                f.write(f"    Speedup: {r['iaif_speedup']:.2f}x\n")
                f.write(f"    Time (orig): {r['iaif_time_orig']:.3f} ms\n")
                f.write(f"    Time (opt):  {r['iaif_time_opt']:.3f} ms\n")
                if not np.isnan(r.get('iaif_g_max_diff', np.nan)):
                    f.write(f"    Max diff (g): {r['iaif_g_max_diff']:.6f}\n")

            f.write("\n")

        f.write("=" * 80 + "\n")

    print(f"✓ Saved report: {output_file}")


def main():
    """
    Main function to test multiple files
    """
    print("=" * 70)
    print("MULTI-FILE VALIDATION TEST")
    print("=" * 70)

    # Find audio files in howtos directory
    audio_dir = Path(__file__).parent.parent.parent / 'howtos'
    audio_files = list(audio_dir.glob('*.wav'))

    if not audio_files:
        print("\n⚠ No audio files found in howtos directory")
        print("  Looking for: ../howtos/*.wav")
        return

    print(f"\nFound {len(audio_files)} audio files:")
    for f in audio_files:
        print(f"  - {f.name}")

    # Test each file
    print(f"\n{'─' * 70}")
    print("RUNNING TESTS")
    print(f"{'─' * 70}\n")

    all_results = []

    for i, audio_file in enumerate(audio_files, 1):
        print(f"[{i}/{len(audio_files)}] Testing: {audio_file.name}...")

        results = test_audio_file(str(audio_file))
        all_results.append(results)

        if results['success']:
            if 'f0_speedup' in results:
                print(f"  ✓ F0: {results['f0_speedup']:.2f}x speedup, "
                      f"{results.get('f0_median_error', np.nan):.2f} Hz median error")
            if 'iaif_speedup' in results:
                print(f"  ✓ IAIF: {results['iaif_speedup']:.2f}x speedup")
        else:
            print(f"  ✗ Error: {results['error']}")

    # Generate report
    print(f"\n{'─' * 70}")
    print("GENERATING REPORT")
    print(f"{'─' * 70}\n")

    generate_report(all_results)

    # Summary
    print(f"\n{'=' * 70}")
    print("SUMMARY")
    print(f"{'=' * 70}\n")

    successful = sum(1 for r in all_results if r['success'])
    f0_speedups = [r['f0_speedup'] for r in all_results if r['success'] and 'f0_speedup' in r]
    iaif_speedups = [r['iaif_speedup'] for r in all_results if r['success'] and 'iaif_speedup' in r]

    print(f"Files tested: {len(all_results)}")
    print(f"Successful: {successful}\n")

    if f0_speedups:
        print(f"F0 Tracking:")
        print(f"  Average speedup: {np.mean(f0_speedups):.2f}x")
        print(f"  Range: {np.min(f0_speedups):.2f}x - {np.max(f0_speedups):.2f}x\n")

    if iaif_speedups:
        print(f"IAIF:")
        print(f"  Average speedup: {np.mean(iaif_speedups):.2f}x")
        print(f"  Range: {np.min(iaif_speedups):.2f}x - {np.max(iaif_speedups):.2f}x\n")

    print(f"{'=' * 70}\n")


if __name__ == "__main__":
    main()
