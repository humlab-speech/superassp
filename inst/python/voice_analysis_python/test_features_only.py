#!/usr/bin/env python3
"""
Test individual optimized features
"""

import numpy as np
import soundfile as sf
import time
import sys
from pathlib import Path
from scipy import signal as scipy_signal

sys.path.insert(0, str(Path(__file__).parent))

from voice_analysis.features import compute_rpde, compute_hnr_nhr, compute_gne

print("Loading test audio...")
audio, fs = sf.read('../a1.wav')
if audio.ndim > 1:
    audio = np.mean(audio, axis=1)
print(f"Audio: {len(audio)} samples, {fs} Hz")

# Test RPDE optimization
print("\n" + "="*60)
print("Testing RPDE Optimization")
print("="*60)
audio_25k = scipy_signal.resample(audio, int(25000 * len(audio) / fs))

print("\n1. RPDE without KD-tree...")
start = time.time()
rpde_baseline = compute_rpde(audio_25k, use_kdtree=False)
t_baseline = time.time() - start
print(f"   Result: {rpde_baseline:.6f}, Time: {t_baseline:.3f}s")

print("\n2. RPDE with KD-tree...")
start = time.time()
rpde_kdtree = compute_rpde(audio_25k, use_kdtree=True)
t_kdtree = time.time() - start
print(f"   Result: {rpde_kdtree:.6f}, Time: {t_kdtree:.3f}s")

if t_baseline > 0:
    speedup = t_baseline / t_kdtree
    print(f"\n   Speedup: {speedup:.2f}x")
    match = abs(rpde_baseline - rpde_kdtree) < 0.01
    print(f"   Results match: {match} (diff: {abs(rpde_baseline - rpde_kdtree):.6f})")

# Test HNR parallelization (shorter)
print("\n" + "="*60)
print("Testing HNR/NHR Parallelization")
print("="*60)

# Use shorter audio segment for faster testing
audio_short = audio[:fs*2]  # 2 seconds

print("\n1. HNR sequential...")
start = time.time()
hnr_seq = compute_hnr_nhr(audio_short, fs, parallel=False)
t_seq = time.time() - start
print(f"   HNR_mean: {hnr_seq['HNR_mean']:.2f}, Time: {t_seq:.3f}s")

print("\n2. HNR parallel...")
start = time.time()
hnr_par = compute_hnr_nhr(audio_short, fs, parallel=True, max_workers=4)
t_par = time.time() - start
print(f"   HNR_mean: {hnr_par['HNR_mean']:.2f}, Time: {t_par:.3f}s")

if t_seq > 0:
    speedup = t_seq / t_par
    print(f"\n   Speedup: {speedup:.2f}x")
    match = abs(hnr_seq['HNR_mean'] - hnr_par['HNR_mean']) < 0.1
    print(f"   Results match: {match} (diff: {abs(hnr_seq['HNR_mean'] - hnr_par['HNR_mean']):.3f})")

# Test GNE parallelization (shorter)
print("\n" + "="*60)
print("Testing GNE Parallelization")
print("="*60)

print("\n1. GNE sequential...")
start = time.time()
gne_seq = compute_gne(audio_short, fs, parallel=False)
t_seq = time.time() - start
print(f"   GNE_mean: {gne_seq['GNE_mean']:.4f}, Time: {t_seq:.3f}s")

print("\n2. GNE parallel...")
start = time.time()
gne_par = compute_gne(audio_short, fs, parallel=True, max_workers=4)
t_par = time.time() - start
print(f"   GNE_mean: {gne_par['GNE_mean']:.4f}, Time: {t_par:.3f}s")

if t_seq > 0:
    speedup = t_seq / t_par
    print(f"\n   Speedup: {speedup:.2f}x")
    match = abs(gne_seq['GNE_mean'] - gne_par['GNE_mean']) < 0.001
    print(f"   Results match: {match} (diff: {abs(gne_seq['GNE_mean'] - gne_par['GNE_mean']):.6f})")

print("\n" + "="*60)
print("All tests completed!")
print("="*60)
