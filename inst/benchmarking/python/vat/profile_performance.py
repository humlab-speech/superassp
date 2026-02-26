#!/usr/bin/env python3
"""
Simple Performance Profile
Identifies bottlenecks in voice analysis
"""

import time
import cProfile
import pstats
import io
import numpy as np
import soundfile as sf
from voice_analysis.core import VoiceAnalyzer

def profile_analysis():
    print("=" * 70)
    print("  Voice Analysis Performance Profile")
    print("=" * 70)
    
    # Load test audio
    print("\nLoading test audio...")
    audio, fs = sf.read('../a1.wav')
    print(f"  Duration: {len(audio)/fs:.2f}s, Sample rate: {fs} Hz")
    
    # Create analyzer
    analyzer = VoiceAnalyzer(use_thesis_mode=False)
    
    # Warm-up run
    print("\nWarm-up run (JIT compilation)...")
    t0 = time.time()
    measures, F0 = analyzer.analyze(audio, fs)
    t1 = time.time()
    print(f"  Time: {t1-t0:.3f}s, Measures: {len(measures)}")
    
    # Profiled run
    print("\nProfiling analysis...")
    profiler = cProfile.Profile()
    profiler.enable()
    
    measures, F0 = analyzer.analyze(audio, fs)
    
    profiler.disable()
    
    # Print statistics
    print("\n" + "=" * 70)
    print("  Top 30 Functions by Cumulative Time")
    print("=" * 70)
    
    s = io.StringIO()
    stats = pstats.Stats(profiler, stream=s)
    stats.sort_stats('cumulative')
    stats.print_stats(30)
    
    output = s.getvalue()
    for line in output.split('\n'):
        if 'voice_analysis' in line or 'function calls' in line or 'ncalls' in line:
            print(line)
    
    print("\n" + "=" * 70)
    print("  Top 30 Functions by Total Time")
    print("=" * 70)
    
    s = io.StringIO()
    stats = pstats.Stats(profiler, stream=s)
    stats.sort_stats('tottime')
    stats.print_stats(30)
    
    output = s.getvalue()
    lines_to_print = []
    for line in output.split('\n')[:50]:  # First 50 lines usually have the important info
        if ('voice_analysis' in line or 'rpde' in line or 'compute' in line or 
            'jitter' in line or 'mfcc' in line or 'function calls' in line or 
            'ncalls' in line or 'tottime' in line):
            lines_to_print.append(line)
    
    for line in lines_to_print[:35]:
        print(line)
    
    # Component timing
    print("\n" + "=" * 70)
    print("  Individual Component Timing (3 runs each)")
    print("=" * 70)
    
    from voice_analysis.features import (
        compute_rpde, compute_dfa, compute_ppe,
        compute_jitter_shimmer_features, compute_hnr_nhr,
        compute_gne, compute_mfcc_features
    )
    from voice_analysis.f0_estimation import estimate_f0_swipe
    
    components = {}
    
    # F0
    print("\n  F0 Estimation (SWIPE)...")
    times = []
    for _ in range(3):
        t0 = time.time()
        F0_test, voicing = estimate_f0_swipe(audio, fs)
        times.append(time.time() - t0)
    components['F0 (SWIPE)'] = np.mean(times)
    print(f"    Mean: {np.mean(times):.3f}s")
    
    # RPDE
    print("  RPDE...")
    times = []
    for _ in range(3):
        t0 = time.time()
        rpde = compute_rpde(audio, fs=fs)
        times.append(time.time() - t0)
    components['RPDE'] = np.mean(times)
    print(f"    Mean: {np.mean(times):.3f}s")
    
    # DFA
    print("  DFA...")
    times = []
    for _ in range(3):
        t0 = time.time()
        dfa = compute_dfa(audio)
        times.append(time.time() - t0)
    components['DFA'] = np.mean(times)
    print(f"    Mean: {np.mean(times):.3f}s")
    
    # PPE
    print("  PPE...")
    times = []
    for _ in range(3):
        t0 = time.time()
        ppe = compute_ppe(F0, fs)
        times.append(time.time() - t0)
    components['PPE'] = np.mean(times)
    print(f"    Mean: {np.mean(times):.3f}s")
    
    # Jitter/Shimmer
    print("  Jitter/Shimmer...")
    times = []
    for _ in range(3):
        t0 = time.time()
        js = compute_jitter_shimmer_features(audio, fs, F0, voicing)
        times.append(time.time() - t0)
    components['Jitter/Shimmer'] = np.mean(times)
    print(f"    Mean: {np.mean(times):.3f}s")
    
    # HNR/NHR
    print("  HNR/NHR...")
    times = []
    for _ in range(3):
        t0 = time.time()
        hnr = compute_hnr_nhr(audio, fs, F0)
        times.append(time.time() - t0)
    components['HNR/NHR'] = np.mean(times)
    print(f"    Mean: {np.mean(times):.3f}s")
    
    # GNE
    print("  GNE...")
    times = []
    for _ in range(3):
        t0 = time.time()
        gne = compute_gne(audio, fs)
        times.append(time.time() - t0)
    components['GNE'] = np.mean(times)
    print(f"    Mean: {np.mean(times):.3f}s")
    
    # MFCC
    print("  MFCC...")
    times = []
    for _ in range(3):
        t0 = time.time()
        mfcc = compute_mfcc_features(audio, fs, F0)
        times.append(time.time() - t0)
    components['MFCC'] = np.mean(times)
    print(f"    Mean: {np.mean(times):.3f}s")
    
    # Summary
    total = sum(components.values())
    print(f"\n  Total component time: {total:.3f}s")
    print(f"\n  Breakdown:")
    for name, t in sorted(components.items(), key=lambda x: x[1], reverse=True):
        pct = 100 * t / total
        bar = '#' * int(pct / 2)
        print(f"    {name:20s}: {t:.3f}s ({pct:5.1f}%) {bar}")
    
    print("\n" + "=" * 70)

if __name__ == '__main__':
    profile_analysis()
