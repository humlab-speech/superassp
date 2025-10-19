# Voice Analysis Toolbox - Python Implementation

A faithful Python reimplementation of the MATLAB Voice Analysis Toolbox by Athanasios Tsanas (2014).

## Features

Computes **132 dysphonia measures** from sustained vowel recordings:

- **Jitter** (22 measures): F0 period perturbations
- **Shimmer** (22 measures): Amplitude perturbations
- **HNR/NHR** (4 measures): Harmonics-to-Noise ratios
- **MFCCs** (84 measures): Mel-frequency cepstral coefficients with deltas
- **Wavelet features** (~50 measures): Multi-resolution analysis
- **GNE** (6 measures): Glottal-to-Noise Excitation
- **PPE** (1 measure): Pitch Period Entropy
- **DFA** (1 measure): Detrended Fluctuation Analysis
- **RPDE** (1 measure): Recurrence Period Density Entropy

## Installation

```bash
pip install -r requirements.txt
```

## Quick Start

### Python API

```python
from voice_analysis import analyze_voice_file

# Analyze a voice recording
measures, F0 = analyze_voice_file('audio.wav')

print(f"Computed {len(measures)} measures")
print(f"Jitter RAP: {measures['jitter_RAP']:.6f}")
print(f"HNR mean: {measures['HNR_mean']:.2f} dB")
```

### Parallel Processing (NEW!)

For improved performance when processing multiple files:

```python
from voice_analysis import analyze_batch_parallel

# Process multiple files in parallel (7-8x faster)
file_list = ['file1.wav', 'file2.wav', 'file3.wav', ...]
results = analyze_batch_parallel(file_list, max_workers=8)

for filename, result in results.items():
    if result is not None:
        measures, F0 = result
        print(f"{filename}: {len(measures)} features")
```

See [PARALLELIZATION_QUICKSTART.md](PARALLELIZATION_QUICKSTART.md) for details.

### Command Line

```bash
# Basic usage
python -m voice_analysis audio.wav

# With custom F0 range
python -m voice_analysis audio.wav --f0-min 75 --f0-max 300

# Save results to JSON
python -m voice_analysis audio.wav --output results.json

# Use PRAAT F0 algorithm instead of SWIPE
python -m voice_analysis audio.wav --f0-algorithm PRAAT
```

## Dependencies

- numpy >= 1.21.0
- scipy >= 1.7.0
- soundfile >= 0.10.0
- librosa >= 0.9.0
- PyWavelets >= 1.1.1
- pySPTK >= 0.1.0 (for SWIPE F0 algorithm)
- nolds >= 0.5.0 (for DFA/RPDE)
- PyEMD >= 1.3.0 (optional, for EMD features)
- numba >= 0.54.0 (optional, for performance)

## Citation

**Required citations:**

1. **Original toolbox:**
   > Tsanas, A., Little, M., McSharry, P., & Ramig, L. (2011). Nonlinear speech analysis algorithms mapped to a standard metric achieve clinically useful quantification of average Parkinson's disease symptom severity. *Journal of the Royal Society Interface*, 8(59), 842-855.

2. **Methodology:**
   > Tsanas, A. (2012). Accurate telemonitoring of Parkinson's disease symptom severity using nonlinear speech signal processing and statistical machine learning. D.Phil. thesis, University of Oxford.

## License

GPL-3.0 (matching the original MATLAB implementation)

## Original MATLAB Code

Copyright (c) Athanasios Tsanas, 2014

## Python Implementation

Python port: 2025

## Validation

This implementation has been validated against the original MATLAB version. Key features:

- F0 estimation using SWIPE (pySPTK) or Praat-style autocorrelation
- All jitter/shimmer variants including Perturbation Quotients
- HNR/NHR computation matching MATLAB output
- MFCCs with deltas and delta-deltas
- Wavelet decomposition features
- Nonlinear measures (DFA, RPDE, PPE)
- GNE measures

## Examples

See `examples/` directory for detailed usage examples.

## Performance and Parallelization

The toolbox includes parallel processing support for improved performance:

- **Single file**: ~4-5 seconds per file (sequential), ~3.9 seconds (parallel with 4 workers)
- **Batch processing**: 7-8x speedup when processing multiple files in parallel

For detailed performance analysis and optimization strategies, see:
- [PARALLELIZATION_QUICKSTART.md](PARALLELIZATION_QUICKSTART.md) - Quick start guide
- [PARALLELIZATION_ANALYSIS.md](PARALLELIZATION_ANALYSIS.md) - Detailed technical analysis
- [PARALLELIZATION_SUMMARY.txt](PARALLELIZATION_SUMMARY.txt) - Performance summary

## Testing

```bash
pytest tests/
```

## Development

```bash
# Install in development mode
pip install -e .

# Run tests with coverage
pytest --cov=voice_analysis tests/

# Format code
black voice_analysis/

# Lint
flake8 voice_analysis/
```

## Notes

- Recommended input: Sustained vowel /a/ recordings (3-5 seconds)
- Sampling rate: 44.1 kHz or higher recommended
- F0 range: Adjust f0_min/f0_max based on speaker (male/female/child)
- DYPSA-dependent features (GQ, VFER) use amplitude envelope fallback

## Support

For questions about the original algorithms, contact the original author:
- Email: tsanasthanasis@gmail.com

For questions about this Python implementation, please open an issue on the repository.
