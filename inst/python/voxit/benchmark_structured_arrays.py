#!/usr/bin/env python3
"""
Benchmark script to compare performance with and without structured arrays.
"""

import numpy as np
import time
from typing import Dict, List

# Import implementations
from voxit_optimized import compute_voxit_features_optimized, _convert_gentle_to_structured, _convert_pitch_to_structured

try:
    from voxit_numba import compute_voxit_features_numba, _convert_gentle_to_arrays, _convert_pitch_to_arrays
    NUMBA_AVAILABLE = True
except ImportError:
    print("Warning: Numba implementation not available")
    NUMBA_AVAILABLE = False


def generate_test_data_dicts(n_words: int = 100) -> Dict:
    """Generate test data as list of dicts (original format)."""
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
    
    # Realistic pitch with some variation
    base_freq = 150.0
    pitch_data = []
    for t in time_vector:
        freq = base_freq * (1 + 0.2 * np.sin(2 * np.pi * t / 2.0))
        pitch_data.append({
            'time': t,
            'frequency': freq
        })
    
    return {
        'gentle_data': gentle_data,
        'pitch_data': pitch_data,
        'duration': total_duration
    }


def generate_test_data_structured(n_words: int = 100) -> Dict:
    """Generate test data as structured arrays (optimized format)."""
    test_data = generate_test_data_dicts(n_words)
    
    # Convert to structured arrays
    n_gentle = len(test_data['gentle_data'])
    gentle_dtype = np.dtype([('start', 'f8'), ('end', 'f8'), ('is_noise', 'bool')])
    gentle_arr = np.zeros(n_gentle, dtype=gentle_dtype)
    
    for i, item in enumerate(test_data['gentle_data']):
        gentle_arr['start'][i] = item['start']
        gentle_arr['end'][i] = item['end']
        gentle_arr['is_noise'][i] = False
    
    n_pitch = len(test_data['pitch_data'])
    pitch_dtype = np.dtype([('time', 'f8'), ('frequency', 'f8')])
    pitch_arr = np.zeros(n_pitch, dtype=pitch_dtype)
    
    for i, item in enumerate(test_data['pitch_data']):
        pitch_arr['time'][i] = item['time']
        pitch_arr['frequency'][i] = item['frequency']
    
    return {
        'gentle_data': gentle_arr,
        'pitch_data': pitch_arr,
        'duration': test_data['duration']
    }


def benchmark_conversion_overhead(n_iterations: int = 100, n_words: int = 100):
    """Benchmark the overhead of converting to structured arrays."""
    print(f"\n{'='*60}")
    print(f"Benchmarking Conversion Overhead ({n_words} words, {n_iterations} iterations)")
    print(f"{'='*60}")
    
    test_data = generate_test_data_dicts(n_words)
    gentle_data = test_data['gentle_data']
    pitch_data = test_data['pitch_data']
    
    # Benchmark gentle conversion
    start = time.perf_counter()
    for _ in range(n_iterations):
        _ = _convert_gentle_to_structured(gentle_data)
    gentle_time = (time.perf_counter() - start) / n_iterations
    
    # Benchmark pitch conversion
    start = time.perf_counter()
    for _ in range(n_iterations):
        _ = _convert_pitch_to_structured(pitch_data)
    pitch_time = (time.perf_counter() - start) / n_iterations
    
    print(f"Gentle data conversion: {gentle_time*1000:.3f} ms")
    print(f"Pitch data conversion:  {pitch_time*1000:.3f} ms")
    print(f"Total conversion:       {(gentle_time + pitch_time)*1000:.3f} ms")
    
    return gentle_time + pitch_time


def benchmark_with_without_structured(n_iterations: int = 50, n_words: int = 100):
    """Compare performance with dict input vs structured array input."""
    print(f"\n{'='*60}")
    print(f"Benchmarking Dict vs Structured Array Input ({n_words} words, {n_iterations} iterations)")
    print(f"{'='*60}")
    
    # Generate both formats
    dict_data = generate_test_data_dicts(n_words)
    struct_data = generate_test_data_structured(n_words)
    
    results = {}
    
    # Test optimized version with dict input
    print("\nOptimized implementation:")
    start = time.perf_counter()
    for _ in range(n_iterations):
        result = compute_voxit_features_optimized(
            dict_data['gentle_data'],
            dict_data['pitch_data']
        )
    dict_time = (time.perf_counter() - start) / n_iterations
    print(f"  Dict input:       {dict_time*1000:.3f} ms")
    results['optimized_dict'] = dict_time
    
    # Test optimized version with structured array input
    start = time.perf_counter()
    for _ in range(n_iterations):
        result = compute_voxit_features_optimized(
            struct_data['gentle_data'],
            struct_data['pitch_data']
        )
    struct_time = (time.perf_counter() - start) / n_iterations
    print(f"  Structured input: {struct_time*1000:.3f} ms")
    print(f"  Speedup:          {dict_time/struct_time:.2f}x")
    results['optimized_struct'] = struct_time
    
    # Test numba version if available
    if NUMBA_AVAILABLE:
        print("\nNumba implementation:")
        
        # Warmup
        _ = compute_voxit_features_numba(dict_data['gentle_data'], dict_data['pitch_data'])
        
        # Test with dict input
        start = time.perf_counter()
        for _ in range(n_iterations):
            result = compute_voxit_features_numba(
                dict_data['gentle_data'],
                dict_data['pitch_data']
            )
        dict_time_numba = (time.perf_counter() - start) / n_iterations
        print(f"  Dict input:       {dict_time_numba*1000:.3f} ms")
        results['numba_dict'] = dict_time_numba
        
        # Test with structured array input
        start = time.perf_counter()
        for _ in range(n_iterations):
            result = compute_voxit_features_numba(
                struct_data['gentle_data'],
                struct_data['pitch_data']
            )
        struct_time_numba = (time.perf_counter() - start) / n_iterations
        print(f"  Structured input: {struct_time_numba*1000:.3f} ms")
        print(f"  Speedup:          {dict_time_numba/struct_time_numba:.2f}x")
        results['numba_struct'] = struct_time_numba
    
    return results


def benchmark_scaling(word_counts: List[int] = [50, 100, 200, 500]):
    """Test how performance scales with input size."""
    print(f"\n{'='*60}")
    print(f"Benchmarking Performance Scaling")
    print(f"{'='*60}")
    
    print(f"\n{'Words':<10} {'Dict (ms)':<15} {'Struct (ms)':<15} {'Speedup':<10}")
    print("-" * 60)
    
    for n_words in word_counts:
        dict_data = generate_test_data_dicts(n_words)
        struct_data = generate_test_data_structured(n_words)
        
        n_iter = max(10, 100 // (n_words // 50))  # Fewer iterations for larger datasets
        
        # Dict input
        start = time.perf_counter()
        for _ in range(n_iter):
            _ = compute_voxit_features_optimized(
                dict_data['gentle_data'],
                dict_data['pitch_data']
            )
        dict_time = (time.perf_counter() - start) / n_iter
        
        # Structured input
        start = time.perf_counter()
        for _ in range(n_iter):
            _ = compute_voxit_features_optimized(
                struct_data['gentle_data'],
                struct_data['pitch_data']
            )
        struct_time = (time.perf_counter() - start) / n_iter
        
        speedup = dict_time / struct_time
        print(f"{n_words:<10} {dict_time*1000:<15.3f} {struct_time*1000:<15.3f} {speedup:<10.2f}x")


def verify_correctness():
    """Verify that structured array version produces same results as dict version."""
    print(f"\n{'='*60}")
    print(f"Verifying Correctness")
    print(f"{'='*60}")
    
    # Generate test data
    dict_data = generate_test_data_dicts(50)
    struct_data = generate_test_data_structured(50)
    
    # Compute with both inputs
    result_dict = compute_voxit_features_optimized(
        dict_data['gentle_data'],
        dict_data['pitch_data']
    )
    
    result_struct = compute_voxit_features_optimized(
        struct_data['gentle_data'],
        struct_data['pitch_data']
    )
    
    # Compare results
    print("\nComparing results:")
    all_close = True
    for key in result_dict:
        val_dict = result_dict[key]
        val_struct = result_struct[key]
        
        if np.isnan(val_dict) and np.isnan(val_struct):
            match = True
        else:
            match = np.allclose(val_dict, val_struct, rtol=1e-10)
        
        status = "✓" if match else "✗"
        print(f"  {status} {key:<40} Dict: {val_dict:12.6f}  Struct: {val_struct:12.6f}")
        
        if not match:
            all_close = False
    
    if all_close:
        print("\n✓ All results match!")
    else:
        print("\n✗ Some results differ!")
    
    return all_close


if __name__ == "__main__":
    print("Voxit Structured Array Performance Benchmark")
    print("=" * 60)
    
    # Verify correctness first
    if not verify_correctness():
        print("\nWarning: Results differ between dict and structured array inputs!")
        print("Continuing with benchmarks anyway...")
    
    # Benchmark conversion overhead
    conversion_time = benchmark_conversion_overhead(n_iterations=1000, n_words=100)
    
    # Benchmark with/without structured arrays
    results = benchmark_with_without_structured(n_iterations=50, n_words=100)
    
    # Scaling benchmark
    benchmark_scaling([50, 100, 200, 500])
    
    # Summary
    print(f"\n{'='*60}")
    print("Summary")
    print(f"{'='*60}")
    print(f"Conversion overhead: {conversion_time*1000:.3f} ms")
    
    if 'optimized_dict' in results and 'optimized_struct' in results:
        speedup = results['optimized_dict'] / results['optimized_struct']
        saved_time = results['optimized_dict'] - results['optimized_struct']
        print(f"\nOptimized implementation:")
        print(f"  Speedup with structured arrays: {speedup:.2f}x")
        print(f"  Time saved per call: {saved_time*1000:.3f} ms")
        
        if saved_time > conversion_time:
            print(f"  ✓ Benefit exceeds conversion overhead after 1 call")
        else:
            calls_needed = int(np.ceil(conversion_time / saved_time))
            print(f"  → Break-even after ~{calls_needed} calls")
    
    if NUMBA_AVAILABLE and 'numba_dict' in results and 'numba_struct' in results:
        speedup = results['numba_dict'] / results['numba_struct']
        saved_time = results['numba_dict'] - results['numba_struct']
        print(f"\nNumba implementation:")
        print(f"  Speedup with structured arrays: {speedup:.2f}x")
        print(f"  Time saved per call: {saved_time*1000:.3f} ms")
        
        if saved_time > conversion_time:
            print(f"  ✓ Benefit exceeds conversion overhead after 1 call")
        else:
            calls_needed = int(np.ceil(conversion_time / saved_time))
            print(f"  → Break-even after ~{calls_needed} calls")
