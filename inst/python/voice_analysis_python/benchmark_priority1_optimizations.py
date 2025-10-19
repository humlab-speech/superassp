#!/usr/bin/env python3
"""
Benchmark Priority 1 Parallelization Optimizations

Tests the impact of:
1. RPDE KD-tree optimization
2. HNR/NHR frame-level parallelization
3. GNE frame-level parallelization
4. Combined optimizations

Compares against baseline implementation.
"""

import numpy as np
import soundfile as sf
import time
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent))

from voice_analysis import VoiceAnalyzer
from voice_analysis.core_parallel import VoiceAnalyzerParallel
from voice_analysis.features import compute_rpde, compute_hnr_nhr, compute_gne


def load_test_audio(filepath='../a1.wav'):
    """Load test audio file"""
    audio, fs = sf.read(filepath)
    if audio.ndim > 1:
        audio = np.mean(audio, axis=1)
    return audio, fs


def benchmark_rpde_optimization(audio, fs):
    """Test RPDE optimization with KD-tree"""
    print("\n" + "="*60)
    print("RPDE OPTIMIZATION BENCHMARK")
    print("="*60)
    
    # Resample for RPDE
    from scipy import signal
    audio_resampled = signal.resample(audio, int(25000 * len(audio) / fs))
    
    # Test without KD-tree (baseline)
    print("\n1. RPDE without KD-tree (baseline Numba)...")
    times_baseline = []
    for i in range(5):
        start = time.time()
        result_baseline = compute_rpde(audio_resampled, m=4, tau=50, epsilon=0.12, 
                                      T_max=1000, use_kdtree=False)
        elapsed = time.time() - start
        times_baseline.append(elapsed)
        print(f"   Run {i+1}: {elapsed:.3f}s, RPDE = {result_baseline:.6f}")
    
    avg_baseline = np.mean(times_baseline)
    std_baseline = np.std(times_baseline)
    print(f"   Average: {avg_baseline:.3f}s ± {std_baseline:.3f}s")
    
    # Test with KD-tree
    print("\n2. RPDE with KD-tree optimization...")
    times_kdtree = []
    for i in range(5):
        start = time.time()
        result_kdtree = compute_rpde(audio_resampled, m=4, tau=50, epsilon=0.12, 
                                    T_max=1000, use_kdtree=True)
        elapsed = time.time() - start
        times_kdtree.append(elapsed)
        print(f"   Run {i+1}: {elapsed:.3f}s, RPDE = {result_kdtree:.6f}")
    
    avg_kdtree = np.mean(times_kdtree)
    std_kdtree = np.std(times_kdtree)
    print(f"   Average: {avg_kdtree:.3f}s ± {std_kdtree:.3f}s")
    
    speedup = avg_baseline / avg_kdtree
    print(f"\n   SPEEDUP: {speedup:.2f}x")
    print(f"   Time saved: {avg_baseline - avg_kdtree:.3f}s ({100*(1-1/speedup):.1f}% reduction)")
    
    return {'baseline': avg_baseline, 'kdtree': avg_kdtree, 'speedup': speedup}


def benchmark_hnr_parallelization(audio, fs):
    """Test HNR frame-level parallelization"""
    print("\n" + "="*60)
    print("HNR/NHR PARALLELIZATION BENCHMARK")
    print("="*60)
    
    # Test sequential
    print("\n1. HNR/NHR sequential (baseline)...")
    times_seq = []
    for i in range(5):
        start = time.time()
        result_seq = compute_hnr_nhr(audio, fs, parallel=False)
        elapsed = time.time() - start
        times_seq.append(elapsed)
        print(f"   Run {i+1}: {elapsed:.3f}s, HNR = {result_seq['HNR_mean']:.2f}")
    
    avg_seq = np.mean(times_seq)
    std_seq = np.std(times_seq)
    print(f"   Average: {avg_seq:.3f}s ± {std_seq:.3f}s")
    
    # Test parallel
    print("\n2. HNR/NHR with frame-level parallelization...")
    times_par = []
    for i in range(5):
        start = time.time()
        result_par = compute_hnr_nhr(audio, fs, parallel=True, max_workers=4)
        elapsed = time.time() - start
        times_par.append(elapsed)
        print(f"   Run {i+1}: {elapsed:.3f}s, HNR = {result_par['HNR_mean']:.2f}")
    
    avg_par = np.mean(times_par)
    std_par = np.std(times_par)
    print(f"   Average: {avg_par:.3f}s ± {std_par:.3f}s")
    
    speedup = avg_seq / avg_par
    print(f"\n   SPEEDUP: {speedup:.2f}x")
    print(f"   Time saved: {avg_seq - avg_par:.3f}s ({100*(1-1/speedup):.1f}% reduction)")
    
    return {'sequential': avg_seq, 'parallel': avg_par, 'speedup': speedup}


def benchmark_gne_parallelization(audio, fs):
    """Test GNE frame-level parallelization"""
    print("\n" + "="*60)
    print("GNE PARALLELIZATION BENCHMARK")
    print("="*60)
    
    # Test sequential
    print("\n1. GNE sequential (baseline)...")
    times_seq = []
    for i in range(3):  # Fewer runs as GNE is slower
        start = time.time()
        result_seq = compute_gne(audio, fs, parallel=False)
        elapsed = time.time() - start
        times_seq.append(elapsed)
        print(f"   Run {i+1}: {elapsed:.3f}s, GNE_mean = {result_seq['GNE_mean']:.4f}")
    
    avg_seq = np.mean(times_seq)
    std_seq = np.std(times_seq)
    print(f"   Average: {avg_seq:.3f}s ± {std_seq:.3f}s")
    
    # Test parallel
    print("\n2. GNE with frame-level parallelization...")
    times_par = []
    for i in range(3):
        start = time.time()
        result_par = compute_gne(audio, fs, parallel=True, max_workers=4)
        elapsed = time.time() - start
        times_par.append(elapsed)
        print(f"   Run {i+1}: {elapsed:.3f}s, GNE_mean = {result_par['GNE_mean']:.4f}")
    
    avg_par = np.mean(times_par)
    std_par = np.std(times_par)
    print(f"   Average: {avg_par:.3f}s ± {std_par:.3f}s")
    
    speedup = avg_seq / avg_par
    print(f"\n   SPEEDUP: {speedup:.2f}x")
    print(f"   Time saved: {avg_seq - avg_par:.3f}s ({100*(1-1/speedup):.1f}% reduction)")
    
    return {'sequential': avg_seq, 'parallel': avg_par, 'speedup': speedup}


def benchmark_full_analysis(audio, fs):
    """Test full voice analysis with all optimizations"""
    print("\n" + "="*60)
    print("FULL ANALYSIS BENCHMARK")
    print("="*60)
    
    # Baseline (original sequential)
    print("\n1. Baseline sequential analysis...")
    analyzer_baseline = VoiceAnalyzer(f0_algorithm='SWIPE')
    times_baseline = []
    for i in range(3):
        start = time.time()
        measures, F0 = analyzer_baseline.analyze(audio, fs)
        elapsed = time.time() - start
        times_baseline.append(elapsed)
        print(f"   Run {i+1}: {elapsed:.3f}s, {len(measures)} features")
    
    avg_baseline = np.mean(times_baseline)
    print(f"   Average: {avg_baseline:.3f}s")
    
    # Parallel without within-feature optimization
    print("\n2. Feature-level parallel (no within-feature optimization)...")
    analyzer_parallel_basic = VoiceAnalyzerParallel(
        f0_algorithm='SWIPE',
        max_workers=4,
        enable_within_feature_parallel=False,
        use_rpde_kdtree=False
    )
    times_parallel_basic = []
    for i in range(3):
        start = time.time()
        measures, F0 = analyzer_parallel_basic.analyze(audio, fs, verbose=False)
        elapsed = time.time() - start
        times_parallel_basic.append(elapsed)
        print(f"   Run {i+1}: {elapsed:.3f}s, {len(measures)} features")
    
    avg_parallel_basic = np.mean(times_parallel_basic)
    speedup_basic = avg_baseline / avg_parallel_basic
    print(f"   Average: {avg_parallel_basic:.3f}s (speedup: {speedup_basic:.2f}x)")
    
    # Parallel with all optimizations
    print("\n3. Feature-level parallel + ALL Priority 1 optimizations...")
    analyzer_optimized = VoiceAnalyzerParallel(
        f0_algorithm='SWIPE',
        max_workers=4,
        enable_within_feature_parallel=True,
        use_rpde_kdtree=True
    )
    times_optimized = []
    for i in range(3):
        start = time.time()
        measures, F0 = analyzer_optimized.analyze(audio, fs, verbose=False)
        elapsed = time.time() - start
        times_optimized.append(elapsed)
        print(f"   Run {i+1}: {elapsed:.3f}s, {len(measures)} features")
    
    avg_optimized = np.mean(times_optimized)
    speedup_optimized = avg_baseline / avg_optimized
    print(f"   Average: {avg_optimized:.3f}s (speedup: {speedup_optimized:.2f}x)")
    
    improvement = (avg_parallel_basic - avg_optimized) / avg_parallel_basic * 100
    print(f"\n   Improvement over basic parallel: {improvement:.1f}%")
    print(f"   Total time saved vs baseline: {avg_baseline - avg_optimized:.3f}s")
    
    return {
        'baseline': avg_baseline,
        'parallel_basic': avg_parallel_basic,
        'optimized': avg_optimized,
        'speedup_basic': speedup_basic,
        'speedup_optimized': speedup_optimized
    }


def main():
    """Run all benchmarks"""
    print("\n" + "="*70)
    print(" PRIORITY 1 PARALLELIZATION OPTIMIZATIONS - COMPREHENSIVE BENCHMARK")
    print("="*70)
    
    # Load audio
    print("\nLoading test audio...")
    audio, fs = load_test_audio()
    print(f"Audio: {len(audio)} samples, {fs} Hz, {len(audio)/fs:.2f}s duration")
    
    results = {}
    
    # Run individual benchmarks
    results['rpde'] = benchmark_rpde_optimization(audio, fs)
    results['hnr'] = benchmark_hnr_parallelization(audio, fs)
    results['gne'] = benchmark_gne_parallelization(audio, fs)
    results['full'] = benchmark_full_analysis(audio, fs)
    
    # Summary
    print("\n" + "="*70)
    print(" SUMMARY OF RESULTS")
    print("="*70)
    
    print("\n1. RPDE Optimization (KD-tree):")
    print(f"   Baseline: {results['rpde']['baseline']:.3f}s")
    print(f"   Optimized: {results['rpde']['kdtree']:.3f}s")
    print(f"   Speedup: {results['rpde']['speedup']:.2f}x")
    
    print("\n2. HNR/NHR Parallelization:")
    print(f"   Sequential: {results['hnr']['sequential']:.3f}s")
    print(f"   Parallel: {results['hnr']['parallel']:.3f}s")
    print(f"   Speedup: {results['hnr']['speedup']:.2f}x")
    
    print("\n3. GNE Parallelization:")
    print(f"   Sequential: {results['gne']['sequential']:.3f}s")
    print(f"   Parallel: {results['gne']['parallel']:.3f}s")
    print(f"   Speedup: {results['gne']['speedup']:.2f}x")
    
    print("\n4. Full Analysis:")
    print(f"   Baseline (sequential): {results['full']['baseline']:.3f}s")
    print(f"   Feature-level parallel: {results['full']['parallel_basic']:.3f}s ({results['full']['speedup_basic']:.2f}x)")
    print(f"   All optimizations: {results['full']['optimized']:.3f}s ({results['full']['speedup_optimized']:.2f}x)")
    
    total_improvement = (1 - results['full']['optimized'] / results['full']['baseline']) * 100
    print(f"\n   TOTAL IMPROVEMENT: {total_improvement:.1f}% faster")
    print(f"   Time per file: {results['full']['baseline']:.2f}s → {results['full']['optimized']:.2f}s")
    
    # Extrapolate to batch processing
    n_files = 100
    batch_baseline = results['full']['baseline'] * n_files
    batch_optimized = results['full']['optimized'] * n_files / 8  # 8-core batch parallelization
    
    print(f"\n5. Expected Batch Performance ({n_files} files, 8 cores):")
    print(f"   Baseline: {batch_baseline/60:.1f} minutes")
    print(f"   Optimized: {batch_optimized/60:.1f} minutes")
    print(f"   Batch speedup: {batch_baseline/batch_optimized:.1f}x")
    
    print("\n" + "="*70)


if __name__ == '__main__':
    main()
