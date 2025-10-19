# Voice Analysis Toolbox - Files Overview

## Documentation Files

### Main Documentation
- **README.md** - Main project documentation with installation and usage
- **PYTHON_IMPLEMENTATION_README.md** - Original implementation notes

### Parallelization Documentation (NEW!)
- **PARALLELIZATION_QUICKSTART.md** - Quick start guide for parallel processing
- **PARALLELIZATION_ANALYSIS.md** - Detailed technical analysis of parallelization opportunities
- **PARALLELIZATION_SUMMARY.txt** - Quick reference summary with performance metrics

### Historical Documentation
- **IMPLEMENTATION_COMPLETE.md** - Original implementation completion notes
- **NUMBA_OPTIMIZATION_ANALYSIS.md** - Numba optimization analysis
- **OPTIMIZATION_GUIDE.md** - General optimization guide
- **PERFORMANCE_SUMMARY.md** - Performance summary from initial implementation
- **COMPLETE_UPDATE.md** - Update notes

### Algorithm Documentation
- **Naylor2007.pdf** - DYPSA algorithm reference
- **d2l-d04.pdf** - DYPSA technical report
- **The_DYPSA_algorithm_for_estimation_of_glottal_clos.pdf** - DYPSA detailed description

## Code Structure

### Main Package (`voice_analysis/`)

```
voice_analysis/
├── __init__.py              # Package exports (VoiceAnalyzer, VoiceAnalyzerParallel, etc.)
├── core.py                  # Main sequential VoiceAnalyzer class
├── core_parallel.py         # NEW! Parallel VoiceAnalyzerParallel class
├── cli.py                   # Command-line interface
│
├── f0_estimation/           # F0 (pitch) estimation
│   ├── __init__.py
│   ├── swipe.py            # SWIPE algorithm (using pySPTK)
│   └── praat.py            # Praat-style autocorrelation
│
├── features/                # Feature computation modules
│   ├── __init__.py
│   ├── jitter_shimmer.py   # Perturbation measures (44 features)
│   ├── hnr.py              # Harmonics-to-Noise Ratio (4 features)
│   ├── mfcc.py             # MFCCs with deltas (78-84 features)
│   ├── wavelet.py          # Wavelet features (~50 features)
│   ├── gne.py              # Glottal-to-Noise Excitation (6 features)
│   ├── ppe.py              # Pitch Period Entropy (1 feature)
│   ├── dfa.py              # Detrended Fluctuation Analysis (1 feature)
│   ├── rpde.py             # Recurrence Period Density Entropy (1 feature)
│   ├── gq.py               # Glottal Quotient (3 features)
│   ├── vfer.py             # VFER measures (7 features)
│   └── emd.py              # Empirical Mode Decomposition (6 features)
│
└── utils/                   # Utility functions
    ├── __init__.py
    ├── dypsa.py            # DYPSA glottal closure detection
    ├── tkeo.py             # Teager-Kaiser Energy Operator
    ├── perturbation.py     # Perturbation quotient computations
    └── entropy.py          # Entropy calculations
```

## Key Files for Users

### For Basic Usage
1. **README.md** - Start here for installation and basic usage
2. **examples/basic_usage.py** - Simple usage example

### For Parallel Processing
1. **PARALLELIZATION_QUICKSTART.md** - How to use parallel features
2. **voice_analysis/core_parallel.py** - Parallel implementation

### For Performance Analysis
1. **PARALLELIZATION_SUMMARY.txt** - Performance overview
2. **benchmark_parallel.py** - Run performance benchmarks
3. **analyze_parallelization.py** - Detailed profiling

### For Development
1. **voice_analysis/core.py** - Main analysis implementation
2. **tests/** - Test suite
3. **setup.py** - Package configuration

## Scripts and Benchmarks

### Analysis Scripts
- **analyze_parallelization.py** - Analyzes parallelization opportunities
- **benchmark_parallel.py** - Benchmarks sequential vs parallel performance
- **benchmark_numba.py** - Benchmarks Numba optimization

### Test Scripts
- **tests/test_voice_analysis.py** - Main test suite
- **tests/test_new_features.py** - Feature-specific tests
- **tests/benchmark_optimization.py** - Optimization benchmarks

### Example Scripts
- **examples/basic_usage.py** - Basic usage example

## Installation and Setup
- **requirements.txt** - Python dependencies
- **setup.py** - Package installation script
- **install.sh** - Installation helper script

## Performance Characteristics

### Computation Time Distribution (4-second audio)
```
RPDE:              2.8s (58.7%) ← BOTTLENECK
MFCC:              0.5s (9.8%)
Glottal Quotient:  0.5s (9.5%)
VFER:              0.4s (9.2%)
F0 Estimation:     0.3s (6.6%)
Jitter:            0.3s (6.0%)
GNE:               0.2s (4.5%)
Others:            0.2s (4.3%)
────────────────────────────────
Total:             4.7s
```

### Parallel Performance
- Sequential: 4.3s per file
- Parallel (4 workers): 3.9s per file (1.08x speedup)
- Batch (8 workers, 100 files): 60s total (7.2x speedup)

## Key Algorithms

### F0 Estimation
- **SWIPE** (default): Sawtooth wave pitch estimator
- **Praat**: Autocorrelation-based

### Time-Domain Features
- **Jitter**: F0 period perturbations (RAP, PQ3, PQ5, PQ11, etc.)
- **Shimmer**: Amplitude perturbations (same variants as jitter)

### Frequency-Domain Features
- **HNR/NHR**: Harmonics-to-noise ratio
- **MFCCs**: Mel-frequency cepstral coefficients with deltas

### Nonlinear Dynamics Features
- **DFA**: Quantifies fractal scaling (Numba-optimized)
- **RPDE**: Recurrence period density entropy (Numba-optimized, but still slow)
- **PPE**: Pitch period entropy

### Voice Quality Features
- **GNE**: Glottal-to-noise excitation
- **Glottal Quotient**: Based on DYPSA algorithm
- **VFER**: Vocal fold excitation ratio

## Optimization Status

### ✅ Implemented Optimizations
- Numba JIT compilation for DFA and RPDE
- Vectorized operations throughout
- Parallel processing at feature and batch level
- Efficient NumPy/SciPy usage

### 🔴 Priority Optimizations Needed
1. **RPDE** - Major bottleneck (58.7% of time)
   - Consider KD-tree for nearest neighbors
   - Cython rewrite of critical loops
   - Expected gain: 2-5x on RPDE → 1.7-2.1x overall

2. **Frame-based parallelization** for HNR/NHR and GNE
   - Expected gain: Additional 1.2-1.5x

### ⚪ Future Considerations
- Julia port for 2.6-3.9x speedup
- GPU acceleration for spectral features
- Advanced caching strategies

## Quick Reference: What to Read When

### "I just want to use the toolbox"
→ README.md → examples/basic_usage.py

### "I need to process many files"
→ PARALLELIZATION_QUICKSTART.md

### "I want to understand performance"
→ PARALLELIZATION_SUMMARY.txt

### "I need to optimize further"
→ PARALLELIZATION_ANALYSIS.md

### "I'm contributing to development"
→ All documentation + source code in voice_analysis/

### "I need algorithm details"
→ PDF files (Naylor2007.pdf, etc.) + source code comments
