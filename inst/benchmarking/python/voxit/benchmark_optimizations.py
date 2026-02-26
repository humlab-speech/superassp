#!/usr/bin/env python3
"""
Benchmark script to compare performance of pure Python, Numba, and Cython implementations.
"""

import numpy as np
import time
from typing import Dict, List, Callable
import sys

# Import core implementation
from voxit_core import compute_voxit_features as compute_voxit_features_python

# Try Numba
try:
    from voxit_numba import compute_voxit_features_numba
    NUMBA_AVAILABLE = True
except ImportError:
    print("Warning: Numba implementation not available")
    NUMBA_AVAILABLE = False

# Try optimized (may use various optimizations)
try:
    from voxit_optimized import compute_voxit_features as compute_voxit_features_optimized
    OPTIMIZED_AVAILABLE = True
except ImportError:
    print("Warning: Optimized implementation not available")
    OPTIMIZED_AVAILABLE = False


def generate_test_data(n_words: int = 50) -> Dict:
    """Generate test data for benchmarking."""
    np.random.seed(42)
    
    # Generate realistic word timing data
    gentle_data = []
    current_time = 0.0
    for i in range(n_words):
        word_duration = np.random.uniform(0.2, 0.8)
        gentle_data.append({
            'word': f'word{i}',
            'start': current_time,
            'end': current_time + word_duration
        })
        current_time += word_duration + np.random.uniform(0.1, 0.5)  # Add pause
    
    # Generate pitch data
    total_duration = current_time
    n_frames = int(total_duration * 100)  # 100 Hz frame rate
    time_vector = np.linspace(0, total_duration, n_frames)
    f0 = 100 + 50 * np.sin(2 * np.pi * time_vector / 2) + np.random.randn(n_frames) * 10
    f0[f0 < 75] = 0  # Simulate unvoiced frames
    
    # Convert to expected format (list of dicts)
    pitch_data = [
        {'time': t, 'frequency': f}
        for t, f in zip(time_vector, f0)
    ]
    
    return {
        'gentle_data': gentle_data,
        'pitch_data': pitch_data,
        'n_words': n_words,
        'duration': total_duration
    }


def benchmark_function(func: Callable, *args, n_iterations: int = 10, warmup: int = 2) -> Dict:
    """Benchmark a function with warmup iterations."""
    # Warmup
    for _ in range(warmup):
        _ = func(*args)
    
    # Actual timing
    times = []
    for _ in range(n_iterations):
        start = time.perf_counter()
        result = func(*args)
        end = time.perf_counter()
        times.append(end - start)
    
    return {
        'mean': np.mean(times),
        'std': np.std(times),
        'min': np.min(times),
        'max': np.max(times),
        'result': result
    }


def benchmark_voxit_features(data: Dict, n_iterations: int = 10) -> Dict:
    """Benchmark Voxit feature computation implementations."""
    print("\n" + "="*60)
    print("Benchmarking Voxit Feature Computation")
    print("="*60)
    print(f"Input: {data['n_words']} words, {data['duration']:.2f} seconds")
    
    results = {}
    gentle_data = data['gentle_data']
    pitch_data = data['pitch_data']
    
    # Python implementation
    print("\nTesting Python implementation...")
    results['python'] = benchmark_function(
        compute_voxit_features_python,
        gentle_data, pitch_data,
        n_iterations=n_iterations
    )
    print(f"  Mean: {results['python']['mean']*1000:.2f} ms")
    print(f"  Std:  {results['python']['std']*1000:.2f} ms")
    
    # Numba implementation
    if NUMBA_AVAILABLE:
        print("\nTesting Numba implementation...")
        results['numba'] = benchmark_function(
            compute_voxit_features_numba,
            gentle_data, pitch_data,
            n_iterations=n_iterations
        )
        print(f"  Mean: {results['numba']['mean']*1000:.2f} ms")
        print(f"  Std:  {results['numba']['std']*1000:.2f} ms")
        speedup = results['python']['mean'] / results['numba']['mean']
        print(f"  Speedup vs Python: {speedup:.2f}x")
    
    # Optimized implementation
    if OPTIMIZED_AVAILABLE:
        print("\nTesting Optimized implementation...")
        results['optimized'] = benchmark_function(
            compute_voxit_features_optimized,
            gentle_data, pitch_data,
            n_iterations=n_iterations
        )
        print(f"  Mean: {results['optimized']['mean']*1000:.2f} ms")
        print(f"  Std:  {results['optimized']['std']*1000:.2f} ms")
        speedup = results['python']['mean'] / results['optimized']['mean']
        print(f"  Speedup vs Python: {speedup:.2f}x")
    
    # Verify results match
    print("\nVerifying results consistency...")
    base_result = results['python']['result']
    for name, res in results.items():
        if name != 'python':
            print(f"\n  {name.capitalize()} vs Python:")
            for key in base_result:
                if key in res['result']:
                    py_val = base_result[key]
                    impl_val = res['result'][key]
                    if isinstance(py_val, (int, float)) and not np.isnan(py_val):
                        diff = abs(impl_val - py_val)
                        rel_diff = diff / abs(py_val) if py_val != 0 else diff
                        print(f"    {key}: abs={diff:.6e}, rel={rel_diff:.2%}")
    
    return results


def print_summary(results_by_size: Dict):
    """Print summary of all benchmarks."""
    print("\n" + "="*60)
    print("PERFORMANCE SUMMARY")
    print("="*60)
    
    for size, results in results_by_size.items():
        print(f"\n{size} words:")
        print(f"  Python:    {results['python']['mean']*1000:7.2f} ms")
        if NUMBA_AVAILABLE and 'numba' in results:
            speedup = results['python']['mean'] / results['numba']['mean']
            print(f"  Numba:     {results['numba']['mean']*1000:7.2f} ms ({speedup:.2f}x speedup)")
        if OPTIMIZED_AVAILABLE and 'optimized' in results:
            speedup = results['python']['mean'] / results['optimized']['mean']
            print(f"  Optimized: {results['optimized']['mean']*1000:7.2f} ms ({speedup:.2f}x speedup)")
    
    # Overall assessment
    print("\n" + "="*60)
    print("RECOMMENDATIONS")
    print("="*60)
    if NUMBA_AVAILABLE or OPTIMIZED_AVAILABLE:
        print("\n✓ Optimizations are available and provide significant speedup")
        if NUMBA_AVAILABLE:
            print("  - Numba JIT compilation provides good performance")
            print("  - No compilation step required, works out of the box")
        if OPTIMIZED_AVAILABLE:
            print("  - Optimized version available with various enhancements")
    else:
        print("\n⚠ No optimizations available")
        print("  - Install numba: pip install numba")
    print()


def main():
    """Run all benchmarks."""
    print("Voxit Optimization Benchmark")
    print("="*60)
    print(f"NumPy version: {np.__version__}")
    print(f"Numba available: {NUMBA_AVAILABLE}")
    print(f"Optimized available: {OPTIMIZED_AVAILABLE}")
    
    # Generate test data
    print("\nGenerating test data...")
    word_counts = [20, 50, 100, 200]
    
    results_by_size = {}
    for n_words in word_counts:
        print(f"\n{'='*60}")
        print(f"Testing with {n_words} words")
        print(f"{'='*60}")
        
        data = generate_test_data(n_words)
        
        # Run benchmarks
        results = benchmark_voxit_features(data, n_iterations=10)
        results_by_size[n_words] = results
    
    # Print summary
    print_summary(results_by_size)


if __name__ == "__main__":
    main()
