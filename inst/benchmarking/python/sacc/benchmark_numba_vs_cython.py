#!/usr/bin/env python3
"""
Direct comparison: Numba vs Cython autocorrelation
"""

import numpy as np
import time

print("="*70)
print("Numba vs Cython: Autocorrelation Performance")
print("="*70)

# Test signal (10 seconds at 8kHz)
sig = np.random.randn(80000)
frm_len = 80
nfrms = 1000
maxlags = 200
win_len = 200

print(f"\nTest parameters:")
print(f"  Signal length: {len(sig)} samples (10 seconds at 8kHz)")
print(f"  Frames: {nfrms}")
print(f"  Lags: {maxlags}")
print(f"  Window: {win_len}")
print()

# Import implementations
print("Loading implementations...")

try:
    from autocorr_numba import autocorr_numba
    print("✓ Numba loaded")
    # Warm up
    print("  Warming up Numba JIT...")
    _ = autocorr_numba(sig[:1000], frm_len, 10, maxlags, win_len, True)
    has_numba = True
except ImportError as e:
    print(f"✗ Numba not available: {e}")
    has_numba = False

try:
    from autocorr_cython import autocorr_cython
    print("✓ Cython loaded")
    has_cython = True
except ImportError as e:
    print(f"✗ Cython not available: {e}")
    has_cython = False

if not has_numba and not has_cython:
    print("\n❌ Neither implementation available!")
    exit(1)

print(f"\n{'='*70}")
print("Running benchmarks...")
print(f"{'='*70}\n")

results = {}

# Benchmark Numba
if has_numba:
    print("Benchmarking Numba...")
    times = []
    for run in range(5):
        start = time.time()
        cor_numba = autocorr_numba(sig, frm_len, nfrms, maxlags, win_len, True)
        elapsed = time.time() - start
        times.append(elapsed)
        print(f"  Run {run+1}: {elapsed:.4f}s")

    avg_numba = np.mean(times[1:])  # Skip first (warmup)
    std_numba = np.std(times[1:])
    results['Numba'] = (avg_numba, std_numba, cor_numba)
    print(f"  Average: {avg_numba:.4f}s ± {std_numba:.4f}s\n")

# Benchmark Cython
if has_cython:
    print("Benchmarking Cython...")
    times = []
    for run in range(5):
        start = time.time()
        cor_cython = autocorr_cython(sig, frm_len, nfrms, maxlags, win_len, True)
        elapsed = time.time() - start
        times.append(elapsed)
        print(f"  Run {run+1}: {elapsed:.4f}s")

    avg_cython = np.mean(times[1:])
    std_cython = np.std(times[1:])
    results['Cython'] = (avg_cython, std_cython, cor_cython)
    print(f"  Average: {avg_cython:.4f}s ± {std_cython:.4f}s\n")

# Results
print(f"{'='*70}")
print("RESULTS")
print(f"{'='*70}\n")

if 'Numba' in results:
    avg, std, _ = results['Numba']
    print(f"Numba:  {avg:.4f}s ± {std:.4f}s")

if 'Cython' in results:
    avg, std, _ = results['Cython']
    print(f"Cython: {avg:.4f}s ± {std:.4f}s")

if 'Numba' in results and 'Cython' in results:
    print()
    if avg_numba < avg_cython:
        speedup = avg_cython / avg_numba
        print(f"🏆 Numba is {speedup:.2f}x faster than Cython")
    else:
        speedup = avg_numba / avg_cython
        print(f"🏆 Cython is {speedup:.2f}x faster than Numba")

    # Numerical accuracy
    max_diff = np.max(np.abs(cor_numba - cor_cython))
    rms_diff = np.sqrt(np.mean((cor_numba - cor_cython)**2))
    print(f"\nNumerical accuracy:")
    print(f"  Max difference: {max_diff:.2e}")
    print(f"  RMS difference: {rms_diff:.2e}")

    if max_diff < 1e-10:
        print(f"  Status: ✓ Perfect agreement")
    elif max_diff < 1e-6:
        print(f"  Status: ✓ Excellent agreement")
    else:
        print(f"  Status: ⚠ Some differences")

# Deployment comparison
print(f"\n{'='*70}")
print("DEPLOYMENT COMPARISON")
print(f"{'='*70}\n")

print(f"{'Feature':<25} {'Numba':<20} {'Cython':<20}")
print("-"*70)
print(f"{'Installation':<25} {'pip install numba':<20} {'C compiler needed':<20}")
print(f"{'Cross-platform':<25} {'Yes':<20} {'Needs recompile':<20}")
print(f"{'Pure Python':<25} {'Yes':<20} {'No (generates C)':<20}")
print(f"{'First-time overhead':<25} {'JIT compile (~1s)':<20} {'None':<20}")
print(f"{'Maintenance':<25} {'Easy':<20} {'Medium':<20}")

# Final recommendation
print(f"\n{'='*70}")
print("RECOMMENDATION")
print(f"{'='*70}\n")

if 'Numba' in results and 'Cython' in results:
    perf_diff = abs(avg_numba - avg_cython) / min(avg_numba, avg_cython)

    if perf_diff < 0.20:  # Less than 20% difference
        print("💡 **Use Numba** - Similar performance, much easier deployment")
        print(f"   - Performance: {avg_numba:.4f}s vs {avg_cython:.4f}s ({perf_diff*100:.1f}% difference)")
        print("   - Installation: pip install numba")
        print("   - No compilation needed")
        print("   - Works on any platform")
    elif avg_cython < avg_numba:
        speedup = avg_numba / avg_cython
        print(f"🚀 **Use Cython** - {speedup:.2f}x faster ({perf_diff*100:.0f}% improvement)")
        print("   - Maximum performance")
        print("   - Trade-off: Requires C compiler, platform-specific builds")
        print("   - Best for: Production where speed is critical")
    else:
        speedup = avg_cython / avg_numba
        print(f"🏆 **Use Numba** - {speedup:.2f}x faster and easier to deploy")
        print("   - Better performance AND easier deployment")
        print("   - Clear winner!")

elif 'Numba' in results:
    print("✓ Numba is available and working")
    print(f"  Performance: {avg_numba:.4f}s")

elif 'Cython' in results:
    print("✓ Cython is available and working")
    print(f"  Performance: {avg_cython:.4f}s")

print()
