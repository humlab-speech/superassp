"""
Parallelization Analysis for Voice Analysis Toolbox

Analyzes current implementation for parallelization opportunities
"""

import time
import numpy as np
import soundfile as sf
from voice_analysis import VoiceAnalyzer
import cProfile
import pstats
from io import StringIO

def profile_analysis(audio_file='../a1.wav'):
    """Profile the voice analysis to identify bottlenecks"""
    
    # Load audio
    audio, fs = sf.read(audio_file)
    if audio.ndim > 1:
        audio = np.mean(audio, axis=1)
    
    # Create profiler
    profiler = cProfile.Profile()
    
    # Run analysis with profiling
    print("Running profiled analysis...")
    analyzer = VoiceAnalyzer(f0_algorithm='SWIPE')
    
    profiler.enable()
    measures, F0 = analyzer.analyze(audio, fs)
    profiler.disable()
    
    # Print results
    print(f"\nTotal measures computed: {len(measures)}")
    
    # Sort by cumulative time
    stats = pstats.Stats(profiler, stream=StringIO())
    stats.strip_dirs()
    stats.sort_stats('cumulative')
    
    print("\n" + "="*80)
    print("TOP 30 FUNCTIONS BY CUMULATIVE TIME")
    print("="*80)
    stats.print_stats(30)
    
    print("\n" + "="*80)
    print("TOP 30 FUNCTIONS BY TOTAL TIME")
    print("="*80)
    stats.sort_stats('tottime')
    stats.print_stats(30)
    
    return stats

def time_individual_features(audio_file='../a1.wav'):
    """Time each feature computation individually"""
    from voice_analysis.f0_estimation import estimate_f0_swipe
    from voice_analysis.features import (
        compute_jitter_shimmer_features,
        compute_hnr_nhr,
        compute_mfcc_features,
        compute_wavelet_features,
        compute_ppe,
        compute_dfa,
        compute_rpde,
        compute_gne,
        compute_emd_features,
        compute_glottal_quotient,
        compute_vfer,
    )
    from scipy import signal as scipy_signal
    
    # Load audio
    audio, fs = sf.read(audio_file)
    if audio.ndim > 1:
        audio = np.mean(audio, axis=1)
    
    # Preprocess
    audio = audio - np.mean(audio)
    audio = audio / (np.max(np.abs(audio)) + 1e-10)
    
    timings = {}
    
    # F0 estimation
    print("Timing F0 estimation...")
    start = time.time()
    F0 = estimate_f0_swipe(audio, fs, 50, 500)
    timings['F0_estimation'] = time.time() - start
    print(f"  F0 estimation: {timings['F0_estimation']:.3f}s")
    
    # Amplitude contour (simple version)
    frame_shift = 0.01
    frame_len = int(frame_shift * fs)
    n_frames = len(audio) // frame_len
    A0 = np.zeros(n_frames)
    for i in range(n_frames):
        start_idx = i * frame_len
        end_idx = min(start_idx + frame_len, len(audio))
        A0[i] = np.max(np.abs(audio[start_idx:end_idx]))
    
    # Individual features
    features_to_time = [
        ('Jitter', lambda: compute_jitter_shimmer_features(F0, 'jitter'), 22),
        ('Shimmer', lambda: compute_jitter_shimmer_features(A0, 'shimmer'), 22),
        ('HNR/NHR', lambda: compute_hnr_nhr(audio, fs, 50, 500), 4),
        ('DFA', lambda: compute_dfa(audio, scales=np.arange(50, 201, 20)), 1),
        ('RPDE', lambda: compute_rpde(scipy_signal.resample(audio, int(25000 * len(audio) / fs)), 
                                       m=4, tau=50, epsilon=0.12, T_max=1000), 1),
        ('PPE', lambda: compute_ppe(F0, fs, f0_mean_healthy=120), 1),
        ('GNE', lambda: compute_gne(audio, fs), 6),
        ('MFCC', lambda: compute_mfcc_features(audio, fs), 78),
        ('Wavelet', lambda: compute_wavelet_features(F0), 50),
        ('Glottal_Quotient', lambda: compute_glottal_quotient(audio, fs, 50, 500), 3),
        ('VFER', lambda: compute_vfer(audio, fs), 7),
        ('EMD', lambda: compute_emd_features(audio), 6),
    ]
    
    total_time = 0
    print("\nTiming individual features:")
    for name, func, count in features_to_time:
        print(f"  {name}...", end=' ', flush=True)
        start = time.time()
        try:
            result = func()
            elapsed = time.time() - start
            timings[name] = elapsed
            total_time += elapsed
            print(f"{elapsed:.3f}s ({count} features, {elapsed/count*1000:.1f}ms per feature)")
        except Exception as e:
            print(f"FAILED: {e}")
            timings[name] = 0
    
    timings['Total'] = total_time
    print(f"\nTotal feature computation time: {total_time:.3f}s")
    print(f"F0 + Features total: {timings['F0_estimation'] + total_time:.3f}s")
    
    # Sort by time
    sorted_timings = sorted([(k, v) for k, v in timings.items() if k != 'Total'], 
                           key=lambda x: x[1], reverse=True)
    
    print("\n" + "="*80)
    print("RANKED BY COMPUTATION TIME")
    print("="*80)
    for name, elapsed in sorted_timings:
        pct = 100 * elapsed / total_time
        print(f"  {name:20s}: {elapsed:6.3f}s ({pct:5.1f}%)")
    
    return timings

def analyze_parallelization_opportunities():
    """Analyze opportunities for parallelization"""
    
    print("\n" + "="*80)
    print("PARALLELIZATION OPPORTUNITIES ANALYSIS")
    print("="*80)
    
    print("\n1. FEATURE-LEVEL PARALLELIZATION (Independent Groups)")
    print("-" * 80)
    
    independent_groups = {
        'Group A - Time series analysis': [
            'Jitter (22 features)',
            'Shimmer (22 features)',
            'PPE (1 feature)',
        ],
        'Group B - Frequency domain analysis': [
            'HNR/NHR (4 features)',
            'GNE (6 features)',
        ],
        'Group C - Nonlinear dynamics': [
            'DFA (1 feature)',
            'RPDE (1 feature)',
        ],
        'Group D - Spectral analysis': [
            'MFCC (78 features)',
            'Wavelet (50 features)',
        ],
        'Group E - Complex analysis': [
            'Glottal Quotient (3 features)',
            'VFER (7 features)',
            'EMD (6 features)',
        ],
    }
    
    print("\nFeatures can be computed in parallel once F0 and A0 are available:")
    for group, features in independent_groups.items():
        print(f"\n  {group}:")
        for feature in features:
            print(f"    - {feature}")
    
    print("\n2. WITHIN-FEATURE PARALLELIZATION")
    print("-" * 80)
    
    parallelizable_features = {
        'Jitter/Shimmer': {
            'opportunity': 'K=3, K=5, K=11 perturbation quotients can be computed in parallel',
            'speedup': '1.5-2x',
            'complexity': 'Low',
        },
        'HNR/NHR': {
            'opportunity': 'Frame-by-frame analysis can be parallelized',
            'speedup': '2-4x on multi-core',
            'complexity': 'Medium',
        },
        'MFCC': {
            'opportunity': 'Librosa internally uses parallel FFT, already optimized',
            'speedup': 'Already optimized',
            'complexity': 'N/A',
        },
        'Wavelet': {
            'opportunity': 'Feature extraction from each detail coefficient independent',
            'speedup': '2-3x',
            'complexity': 'Low',
        },
        'DFA': {
            'opportunity': 'Scale-by-scale computation can be parallelized',
            'speedup': '2-4x',
            'complexity': 'Medium (already Numba optimized)',
        },
        'RPDE': {
            'opportunity': 'Distance computations can be vectorized/parallelized',
            'speedup': '2-5x',
            'complexity': 'High (already Numba optimized)',
        },
        'GNE': {
            'opportunity': 'Multiple frequency band analysis can be parallelized',
            'speedup': '3-6x',
            'complexity': 'Medium',
        },
        'EMD': {
            'opportunity': 'Feature extraction from IMFs can be parallelized',
            'speedup': '1.5-2x',
            'complexity': 'Low',
        },
    }
    
    print("\nWithin-feature parallelization opportunities:")
    for feature, info in parallelizable_features.items():
        print(f"\n  {feature}:")
        print(f"    Opportunity: {info['opportunity']}")
        print(f"    Potential speedup: {info['speedup']}")
        print(f"    Implementation complexity: {info['complexity']}")
    
    print("\n3. RECOMMENDED PARALLELIZATION STRATEGY")
    print("-" * 80)
    
    recommendations = [
        ("Priority 1: Feature-level parallelization", 
         "Use concurrent.futures.ProcessPoolExecutor or ThreadPoolExecutor to run independent feature groups in parallel",
         "Expected speedup: 3-5x on 4-8 core systems",
         "Effort: Low-Medium, just need proper dependency management"),
        
        ("Priority 2: Within-feature parallelization for heavy features",
         "Focus on HNR/NHR, GNE (frame-based), and DFA/RPDE (if not using Numba)",
         "Expected speedup: Additional 1.5-2x for these specific features",
         "Effort: Medium, requires careful refactoring"),
        
        ("Priority 3: Batch processing optimization",
         "When processing multiple files, parallelize at file level",
         "Expected speedup: Linear with number of cores",
         "Effort: Low, straightforward parallelization"),
    ]
    
    for i, (title, desc, speedup, effort) in enumerate(recommendations, 1):
        print(f"\n  {title}:")
        print(f"    Description: {desc}")
        print(f"    {speedup}")
        print(f"    {effort}")
    
    print("\n4. IMPLEMENTATION CONSTRAINTS")
    print("-" * 80)
    print("""
  Threading vs Multiprocessing:
    - ThreadPoolExecutor: Good for I/O and numpy operations (releases GIL)
    - ProcessPoolExecutor: Better for pure Python code, but overhead from pickling
    - Recommendation: Start with ThreadPoolExecutor for feature-level parallelization
  
  Numba considerations:
    - Already optimized DFA and RPDE with @jit decorators
    - Numba releases GIL, works well with threading
    - Don't parallelize already-Numba-optimized inner loops
  
  Memory considerations:
    - Each parallel worker needs its own copy of audio data
    - Feature-level parallelization has minimal memory overhead
    - Frame-level parallelization could increase memory usage significantly
    """)
    
    print("\n5. EXPECTED OVERALL SPEEDUP")
    print("-" * 80)
    print("""
  Current state:
    - Sequential execution of all features
    - Some Numba optimization (DFA, RPDE)
    - Total time: ~5-15 seconds per file (depending on length)
  
  With feature-level parallelization (4 cores):
    - Expected speedup: 3-4x
    - Total time: ~2-5 seconds per file
  
  With additional within-feature parallelization:
    - Expected speedup: 4-6x total
    - Total time: ~1-3 seconds per file
  
  Batch processing (100 files, 8 cores):
    - Expected speedup: 6-7x per file
    - Total time: ~20-50 seconds for 100 files
    """)

if __name__ == '__main__':
    print("Voice Analysis Toolbox - Parallelization Analysis")
    print("=" * 80)
    
    # Time individual features
    print("\n" + "="*80)
    print("PHASE 1: INDIVIDUAL FEATURE TIMING")
    print("="*80)
    timings = time_individual_features()
    
    # Analyze parallelization opportunities
    analyze_parallelization_opportunities()
    
    print("\n" + "="*80)
    print("ANALYSIS COMPLETE")
    print("="*80)
