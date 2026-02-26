#!/usr/bin/env python3
"""
Comprehensive benchmark: Pure Python vs Numba vs Cython

Compares all three implementations of autocorrelation
"""

import numpy as np
import time

print("="*70)
print("SAcC Autocorrelation Performance Comparison")
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
implementations = {}

# 1. Pure Python (from original sbpca.py)
try:
    import sbpca
    # Use original autocorr function
    from sbpca import autocorr as autocorr_python
    implementations['Pure Python'] = autocorr_python
    print("✓ Pure Python implementation loaded")
except ImportError as e:
    print(f"✗ Pure Python not available: {e}")

# 2. Numba
try:
    from autocorr_numba import autocorr_numba, is_numba_available
    if is_numba_available():
        implementations['Numba'] = autocorr_numba
        print("✓ Numba implementation loaded")
        # Warm up JIT
        print("  (Warming up Numba JIT...)")
        _ = autocorr_numba(sig[:1000], frm_len, 10, maxlags, win_len, True)
    else:
        print("✗ Numba not available (pip install numba)")
except ImportError as e:
    print(f"✗ Numba not available: {e}")

# 3. Cython
try:
    from autocorr_cython import autocorr_cython
    implementations['Cython'] = autocorr_cython
    print("✓ Cython implementation loaded")
except ImportError as e:
    print(f"✗ Cython not available: {e}")
    print("  (Run: python setup_cython.py build_ext --inplace)")

if len(implementations) == 0:
    print("\n❌ No implementations available to benchmark!")
    exit(1)

print(f"\n{'='*70}")
print("Running benchmarks...")
print(f"{'='*70}\n")

# Benchmark each implementation
results = {}

for name, func in implementations.items():
    print(f"Benchmarking {name}...")

    times = []
    for run in range(5):
        start = time.time()
        cor = func(sig, frm_len, nfrms, maxlags, win_len, True)
        elapsed = time.time() - start
        times.append(elapsed)
        print(f"  Run {run+1}: {elapsed:.4f}s")

    avg_time = np.mean(times[1:])  # Skip first run
    std_time = np.std(times[1:])
    results[name] = {
        'avg': avg_time,
        'std': std_time,
        'cor': cor
    }

    print(f"  Average: {avg_time:.4f}s ± {std_time:.4f}s\n")

# Calculate speedups
if 'Pure Python' in results:
    baseline = results['Pure Python']['avg']
else:
    baseline = max(r['avg'] for r in results.values())

print(f"{'='*70}")
print("RESULTS SUMMARY")
print(f"{'='*70}\n")

print(f"{'Implementation':<20} {'Time (s)':<12} {'Speedup':<10} {'Status'}")
print("-"*70)

for name in ['Pure Python', 'Numba', 'Cython']:
    if name in results:
        avg = results[name]['avg']
        speedup = baseline / avg
        status = "✓"
        if speedup >= 3.0:
            status = "🚀 Fastest"
        elif speedup >= 2.0:
            status = "⚡ Very Fast"
        elif speedup >= 1.5:
            status = "✓ Fast"
        print(f"{name:<20} {avg:>6.4f}s      {speedup:>5.2f}x     {status}")

# Verify numerical accuracy
print(f"\n{'='*70}")
print("NUMERICAL ACCURACY CHECK")
print(f"{'='*70}\n")

if len(results) > 1:
    # Compare all implementations
    names = list(results.keys())
    reference = results[names[0]]['cor']

    for name in names[1:]:
        test = results[name]['cor']
        max_diff = np.max(np.abs(reference - test))
        rms_diff = np.sqrt(np.mean((reference - test)**2))

        print(f"{names[0]} vs {name}:")
        print(f"  Max difference: {max_diff:.2e}")
        print(f"  RMS difference: {rms_diff:.2e}")

        if max_diff < 1e-6:
            print(f"  Status: ✓ Excellent agreement")
        elif max_diff < 1e-3:
            print(f"  Status: ✓ Good agreement")
        else:
            print(f"  Status: ⚠ Differences detected")
        print()
else:
    print("Only one implementation available - cannot compare")

# Final recommendation
print(f"{'='*70}")
print("RECOMMENDATION")
print(f"{'='*70}\n")

if 'Cython' in results and 'Numba' in results:
    cython_time = results['Cython']['avg']
    numba_time = results['Numba']['avg']

    if cython_time < numba_time * 0.8:  # Cython significantly faster
        print("🏆 **Use Cython** for maximum performance")
        print(f"   - {baseline/cython_time:.1f}x faster than baseline")
        print(f"   - {numba_time/cython_time:.1f}x faster than Numba")
        print("   - Requires: C compiler, one-time compilation")
        print("   - Best for: Production deployments prioritizing speed")
    elif numba_time < cython_time * 0.8:  # Numba significantly faster
        print("🏆 **Use Numba** for best balance")
        print(f"   - {baseline/numba_time:.1f}x faster than baseline")
        print(f"   - {cython_time/numba_time:.1f}x faster than Cython")
        print("   - Requires: pip install numba (no C compiler)")
        print("   - Best for: Easy deployment, good performance")
    else:  # Similar performance
        print("💡 **Use Numba** for ease of deployment")
        print(f"   - Similar performance to Cython ({numba_time:.4f}s vs {cython_time:.4f}s)")
        print("   - Much easier to deploy (pip install, no compilation)")
        print("   - Recommended unless you need absolute maximum speed")

elif 'Cython' in results:
    cython_time = results['Cython']['avg']
    print(f"🏆 **Cython**: {baseline/cython_time:.1f}x speedup")
    print("   Best performance, requires C compiler")

elif 'Numba' in results:
    numba_time = results['Numba']['avg']
    print(f"🏆 **Numba**: {baseline/numba_time:.1f}x speedup")
    print("   Great performance, easy to install")

print()
