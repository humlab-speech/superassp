#!/usr/bin/env python3
"""
Comprehensive Performance Benchmark
Voice Analysis Toolbox Python Implementation
"""

import time
import numpy as np
import soundfile as sf
from voice_analysis.core import VoiceAnalyzer
from voice_analysis.features.rpde import CYTHON_AVAILABLE, NUMBA_AVAILABLE, KDTREE_AVAILABLE
import sys

def print_header(title):
    print("\n" + "=" * 70)
    print(f"  {title}")
    print("=" * 70)

def print_section(title):
    print(f"\n{title}")
    print("-" * 70)

def main():
    print_header("Voice Analysis Performance Benchmark")
    
    # Check environment
    print_section("Environment Check")
    print(f"Cython RPDE Available: {CYTHON_AVAILABLE}")
    print(f"Numba JIT Available:   {NUMBA_AVAILABLE}")
    print(f"KD-Tree Available:     {KDTREE_AVAILABLE}")
    
    try:
        import pywt
        print(f"PyWavelets Available:  True (version {pywt.__version__})")
    except ImportError:
        print(f"PyWavelets Available:  False")
    
    try:
        from PyEMD import EMD
        print(f"PyEMD Available:       True")
    except ImportError:
        print(f"PyEMD Available:       False")
    
    # Load test file
    print_section("Loading Test Data")
    try:
        audio, fs = sf.read('../a1.wav')
        duration = len(audio) / fs
        print(f"File: a1.wav")
        print(f"Sample rate: {fs} Hz")
        print(f"Duration: {duration:.2f} seconds")
        print(f"Samples: {len(audio)}")
    except Exception as e:
        print(f"ERROR: Could not load test file: {e}")
        sys.exit(1)
    
    # Test 1: Warm-up run (JIT compilation)
    print_section("Test 1: Warm-up Run (JIT Compilation)")
    analyzer = VoiceAnalyzer(use_thesis_mode=False)
    print("Performing warm-up analysis...")
    t0 = time.time()
    measures, F0 = analyzer.analyze(audio, fs)
    t1 = time.time()
    warmup_time = t1 - t0
    print(f"Warm-up time: {warmup_time:.3f}s")
    print(f"Measures computed: {len(measures)}")
    
    # Test 2: Cached runs (true performance)
    print_section("Test 2: Cached Performance (5 runs)")
    times = []
    for i in range(5):
        t0 = time.time()
        measures, F0 = analyzer.analyze(audio, fs)
        t1 = time.time()
        runtime = t1 - t0
        times.append(runtime)
        print(f"  Run {i+1}: {runtime:.3f}s")
    
    mean_time = np.mean(times)
    std_time = np.std(times)
    median_time = np.median(times)
    min_time = np.min(times)
    max_time = np.max(times)
    
    print(f"\nStatistics:")
    print(f"  Mean:   {mean_time:.3f}s ± {std_time:.3f}s")
    print(f"  Median: {median_time:.3f}s")
    print(f"  Min:    {min_time:.3f}s")
    print(f"  Max:    {max_time:.3f}s")
    print(f"  Analysis rate: {duration/mean_time:.2f}x realtime")
    
    # Test 3: Component-level profiling
    print_section("Test 3: Component Profiling")
    
    # Import specific functions
    from voice_analysis.utils.f0_extraction import extract_f0
    from voice_analysis.features.rpde import compute_rpde
    from voice_analysis.features.dfa import compute_dfa
    from voice_analysis.features.ppe import compute_ppe
    from voice_analysis.features.jitter_shimmer import compute_jitter_shimmer_features
    from voice_analysis.features.hnr import compute_hnr_nhr
    from voice_analysis.features.gne import compute_gne
    from voice_analysis.features.mfcc import compute_mfcc_features
    
    component_times = {}
    
    # F0 Extraction
    print("  Testing F0 extraction...")
    t0 = time.time()
    F0, voicing = extract_f0(audio, fs, algorithm='SWIPE')
    t1 = time.time()
    component_times['F0 Extraction'] = t1 - t0
    print(f"    F0 extraction: {t1-t0:.3f}s")
    
    # RPDE
    print("  Testing RPDE...")
    t0 = time.time()
    rpde_val = compute_rpde(audio, fs=fs)
    t1 = time.time()
    component_times['RPDE'] = t1 - t0
    print(f"    RPDE: {t1-t0:.3f}s")
    
    # DFA
    print("  Testing DFA...")
    t0 = time.time()
    dfa_val = compute_dfa(audio)
    t1 = time.time()
    component_times['DFA'] = t1 - t0
    print(f"    DFA: {t1-t0:.3f}s")
    
    # PPE
    print("  Testing PPE...")
    t0 = time.time()
    ppe_val = compute_ppe(F0, fs)
    t1 = time.time()
    component_times['PPE'] = t1 - t0
    print(f"    PPE: {t1-t0:.3f}s")
    
    # Jitter/Shimmer
    print("  Testing Jitter/Shimmer...")
    t0 = time.time()
    jit_shim = compute_jitter_shimmer_features(audio, fs, F0, voicing)
    t1 = time.time()
    component_times['Jitter/Shimmer'] = t1 - t0
    print(f"    Jitter/Shimmer: {t1-t0:.3f}s")
    
    # HNR/NHR
    print("  Testing HNR/NHR...")
    t0 = time.time()
    hnr_nhr = compute_hnr_nhr(audio, fs, F0)
    t1 = time.time()
    component_times['HNR/NHR'] = t1 - t0
    print(f"    HNR/NHR: {t1-t0:.3f}s")
    
    # GNE
    print("  Testing GNE...")
    t0 = time.time()
    gne = compute_gne(audio, fs)
    t1 = time.time()
    component_times['GNE'] = t1 - t0
    print(f"    GNE: {t1-t0:.3f}s")
    
    # MFCC
    print("  Testing MFCC...")
    t0 = time.time()
    mfcc = compute_mfcc_features(audio, fs, F0)
    t1 = time.time()
    component_times['MFCC'] = t1 - t0
    print(f"    MFCC: {t1-t0:.3f}s")
    
    # Calculate percentages
    total_component_time = sum(component_times.values())
    print(f"\n  Total component time: {total_component_time:.3f}s")
    print(f"\n  Component breakdown:")
    sorted_components = sorted(component_times.items(), key=lambda x: x[1], reverse=True)
    for name, t in sorted_components:
        pct = 100 * t / total_component_time
        print(f"    {name:20s}: {t:.3f}s ({pct:.1f}%)")
    
    # Test 4: Test Cython vs Numba for RPDE
    print_section("Test 4: RPDE Implementation Comparison")
    
    if CYTHON_AVAILABLE:
        print("  Testing Cython RPDE (3 runs)...")
        cython_times = []
        for i in range(3):
            t0 = time.time()
            rpde_val = compute_rpde(audio, fs=fs, use_cython=True, use_kdtree=False)
            t1 = time.time()
            cython_times.append(t1 - t0)
        cython_mean = np.mean(cython_times)
        print(f"    Cython mean: {cython_mean:.3f}s")
    else:
        print("  Cython not available")
        cython_mean = None
    
    if NUMBA_AVAILABLE:
        print("  Testing Numba RPDE (3 runs)...")
        numba_times = []
        for i in range(3):
            t0 = time.time()
            rpde_val = compute_rpde(audio, fs=fs, use_cython=False, use_kdtree=False)
            t1 = time.time()
            numba_times.append(t1 - t0)
        numba_mean = np.mean(numba_times)
        print(f"    Numba mean: {numba_mean:.3f}s")
    else:
        print("  Numba not available")
        numba_mean = None
    
    if cython_mean and numba_mean:
        speedup = numba_mean / cython_mean
        print(f"    Cython speedup: {speedup:.2f}x vs Numba")
    
    # Summary
    print_section("Summary")
    print(f"Best performance: {min_time:.3f}s per file ({len(measures)} measures)")
    print(f"Analysis rate: {duration/min_time:.2f}x realtime")
    print(f"Throughput: {1/min_time:.2f} files/second")
    print(f"\nBatch estimates:")
    print(f"  10 files:   {10*min_time:.1f}s sequential")
    print(f"  100 files:  {100*min_time:.1f}s sequential")
    print(f"  1000 files: {1000*min_time/60:.1f}min sequential")
    
    # Performance targets
    print_section("Performance vs Targets")
    print(f"Current:        {mean_time:.2f}s")
    print(f"Previous best:  3.98s  {'✓ FASTER' if mean_time < 3.98 else '✗ SLOWER'}")
    print(f"Target:         3.00s  {'✓ ACHIEVED' if mean_time < 3.00 else '✗ NOT YET'}")
    print(f"Stretch goal:   2.00s  {'✓ ACHIEVED' if mean_time < 2.00 else '✗ NOT YET'}")
    
    # Recommendations
    print_section("Recommendations")
    
    if mean_time > 4.0:
        print("⚠️  Performance is below baseline. Check:")
        print("   - Cython extensions built and loaded")
        print("   - Numba caching enabled")
        print("   - No debug mode active")
    elif mean_time > 3.0:
        print("✓ Good performance. Minor optimizations possible:")
        print("   - Profile to identify remaining bottlenecks")
        print("   - Consider parallel feature groups")
    elif mean_time > 2.0:
        print("✓✓ Excellent performance. Near optimal:")
        print("   - Consider platform-specific SIMD")
        print("   - Memory pooling for batch processing")
    else:
        print("✓✓✓ Outstanding performance achieved!")
    
    print("\n" + "=" * 70)

if __name__ == '__main__':
    main()
