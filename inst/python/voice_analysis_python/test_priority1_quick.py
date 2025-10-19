#!/usr/bin/env python3
"""
Quick test of Priority 1 optimizations
Tests basic functionality without extensive benchmarking
"""

import numpy as np
import soundfile as sf
import time
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from voice_analysis import VoiceAnalyzer
from voice_analysis.core_parallel import VoiceAnalyzerParallel

print("Loading test audio...")
audio, fs = sf.read('../a1.wav')
if audio.ndim > 1:
    audio = np.mean(audio, axis=1)

print(f"Audio: {len(audio)} samples, {fs} Hz, {len(audio)/fs:.2f}s")

# Test 1: Baseline
print("\n1. Baseline sequential...")
analyzer_baseline = VoiceAnalyzer(f0_algorithm='SWIPE')
start = time.time()
measures_baseline, F0 = analyzer_baseline.analyze(audio, fs)
time_baseline = time.time() - start
print(f"   Time: {time_baseline:.3f}s, Features: {len(measures_baseline)}")

# Test 2: Basic parallel
print("\n2. Feature-level parallel (no optimizations)...")
analyzer_basic = VoiceAnalyzerParallel(
    f0_algorithm='SWIPE',
    max_workers=4,
    enable_within_feature_parallel=False,
    use_rpde_kdtree=False
)
start = time.time()
measures_basic, F0 = analyzer_basic.analyze(audio, fs, verbose=False)
time_basic = time.time() - start
speedup_basic = time_baseline / time_basic
print(f"   Time: {time_basic:.3f}s, Speedup: {speedup_basic:.2f}x")

# Test 3: All optimizations
print("\n3. All Priority 1 optimizations...")
analyzer_optimized = VoiceAnalyzerParallel(
    f0_algorithm='SWIPE',
    max_workers=4,
    enable_within_feature_parallel=True,
    use_rpde_kdtree=True
)
start = time.time()
measures_opt, F0 = analyzer_optimized.analyze(audio, fs, verbose=False)
time_opt = time.time() - start
speedup_opt = time_baseline / time_opt
print(f"   Time: {time_opt:.3f}s, Speedup: {speedup_opt:.2f}x")

# Summary
print("\n" + "="*60)
print("SUMMARY")
print("="*60)
print(f"Baseline:           {time_baseline:.2f}s (1.00x)")
print(f"Basic parallel:     {time_basic:.2f}s ({speedup_basic:.2f}x)")
print(f"All optimizations:  {time_opt:.2f}s ({speedup_opt:.2f}x)")
print(f"\nImprovement: {(time_baseline - time_opt):.2f}s saved ({100*(1-time_opt/time_baseline):.1f}% faster)")

# Verify correctness
print("\n" + "="*60)
print("CORRECTNESS CHECK")
print("="*60)
n_match = sum(1 for k in measures_baseline if k in measures_opt and 
               abs(measures_baseline.get(k, 0) - measures_opt.get(k, 0)) < 1e-4)
print(f"Matching features: {n_match}/{len(measures_baseline)}")
if n_match >= len(measures_baseline) * 0.95:
    print("✓ Results are consistent (>95% match)")
else:
    print("⚠ Warning: Some results differ")
