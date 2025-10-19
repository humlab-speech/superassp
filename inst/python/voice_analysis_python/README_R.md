# Voice Analysis Toolbox - R Integration

This directory contains the Python implementation of the Voice Analysis Toolbox, a faithful reimplementation of the MATLAB Voice Analysis Toolbox by Athanasios Tsanas for computing 132 dysphonia measures from sustained vowel recordings.

## Installation from R

### Automatic Installation

The easiest way to install the module is through R:

```r
library(superassp)

# Install with automatic method selection (tries Cython, falls back to pure Python)
install_voice_analysis()

# Or force Cython build for maximum performance (requires C compiler)
install_voice_analysis(method = "cython")

# Or install pure Python/Numba version (no compilation needed)
install_voice_analysis(method = "pure")
```

### Manual Installation

Alternatively, you can install manually:

```bash
# Navigate to this directory
cd inst/python/voice_analysis_python

# Install dependencies
pip install -r requirements.txt

# Install with Cython optimizations (recommended)
python setup_cython.py install

# Or install without Cython
python setup.py install
```

## Usage from R

### Basic Usage

```r
library(superassp)

# Analyze a sustained vowel recording
result <- lst_vat("sustained_vowel.wav")

# View all 132 measures
names(result$measures)

# Access specific measures
result$measures$DFA    # Detrended Fluctuation Analysis
result$measures$RPDE   # Recurrence Period Density Entropy
result$measures$PPE    # Pitch Period Entropy
```

### Advanced Options

```r
# Custom F0 range (e.g., for male voice)
result <- lst_vat("male_vowel.wav", f0_min = 75, f0_max = 300)

# Time windowing - analyze 1.0 to 3.0 seconds
result <- lst_vat("recording.wav", beginTime = 1.0, endTime = 3.0)

# Parallel processing for speed
result <- lst_vat("vowel.wav", n_cores = 8)

# Different F0 algorithm
result <- lst_vat("vowel.wav", f0_algorithm = "PRAAT")  # Default is "SWIPE"

# Return F0 contour for inspection
result <- lst_vat("vowel.wav", return_f0 = TRUE)
plot(result$f0, type = "l", main = "F0 Contour")

# Thesis mode (for replication of Tsanas 2012)
result <- lst_vat("vowel.wav", use_thesis_mode = TRUE)
```

### Batch Processing

```r
# Process multiple files
files <- c("vowel1.wav", "vowel2.wav", "vowel3.wav")
results <- lst_vat(files)

# Extract specific measure from all files
dfa_values <- sapply(results, function(r) r$measures$DFA)
```

### Non-WAV Formats

The function automatically handles any media format supported by FFmpeg:

```r
# Process video file (extracts audio automatically)
result <- lst_vat("interview.mp4", beginTime = 10, endTime = 13)

# Process MP3
result <- lst_vat("audio.mp3")
```

## Feature Categories

The toolbox computes 132 measures across 9 categories:

1. **Jitter (22-25 measures)**: Fundamental frequency perturbation
   - Local jitter, RAP, PPQ5, DDP
   - Absolute and relative variants
   - Optional: AR-based jitter (thesis mode)

2. **Shimmer (22-25 measures)**: Amplitude perturbation
   - Local shimmer, APQ3, APQ5, APQ11, DDA
   - dB and percentage variants
   - Optional: NMSP (thesis mode)

3. **Harmonic-to-Noise (4 measures)**:
   - HNR (0-500 Hz, 0-1500 Hz)
   - NHR (0-500 Hz, 0-1500 Hz)

4. **Nonlinear Dynamics (3 measures)**:
   - DFA (Detrended Fluctuation Analysis)
   - RPDE (Recurrence Period Density Entropy)
   - PPE (Pitch Period Entropy)

5. **Glottal Measures (9 measures)**:
   - GNE (Glottal-to-Noise Excitation, 6 variants)
   - Glottal Quotient (3 variants)

6. **VFER (7 measures)**: Vocal Fold Excitation Ratio

7. **MFCCs (84 measures)**: Mel-Frequency Cepstral Coefficients
   - 12 MFCCs with delta and delta-delta (36 total)
   - Statistics: mean, std, min, max, range, skewness, kurtosis

8. **Wavelet (~50 measures)**: Multi-scale wavelet decomposition
   - Energy and entropy at multiple scales

9. **EMD (6 measures)**: Empirical Mode Decomposition
   - Energy and frequency of intrinsic mode functions

## Performance

For a typical 3-second sustained vowel:
- Sequential processing: 5-15 seconds
- Parallel processing (8 cores): 2-5 seconds
- With Cython optimizations: 2-3x faster

### Optimization Tips

```r
# Check if Cython extensions are available
info <- voice_analysis_info()
print(info$cython_available)  # Should be TRUE for best performance

# Use recommended number of workers
print(info$recommended_workers)  # Optimal for your system

# Apply parallelization
result <- lst_vat("vowel.wav", n_cores = info$recommended_workers)
```

## System Requirements

### R Packages
- superassp (this package)
- reticulate (for Python integration)
- av (for media file loading)

### Python Dependencies
- numpy >= 1.21.0
- scipy >= 1.7.0
- soundfile >= 0.10.0
- librosa >= 0.9.0
- pywt >= 1.1.1
- pysptk >= 0.1.0
- nolds >= 0.5.0
- EMD-signal >= 1.3.0
- numba >= 0.54.0 (for performance)
- joblib >= 1.0.0 (for parallelization)

### Optional (for Cython)
- cython >= 0.29.0
- C compiler (gcc, clang, or MSVC)

## Troubleshooting

### Module Not Available

```r
# Check if module is installed
voice_analysis_available()

# If FALSE, install it
install_voice_analysis()
```

### Cython Compilation Failed

If Cython compilation fails, the installer will automatically fall back to pure Python:

```r
# Or explicitly install pure Python version
install_voice_analysis(method = "pure")
```

### Import Errors

```r
# Reinstall with dependencies
install_voice_analysis(force_reinstall = TRUE)
```

### Performance Issues

```r
# Check system configuration
info <- voice_analysis_info()
print(info)

# Use parallel processing
result <- lst_vat("vowel.wav", n_cores = info$recommended_workers)
```

## Citation

If you use this toolbox, please cite:

> Tsanas, A., Little, M., McSharry, P., & Ramig, L. (2011).
> Nonlinear speech analysis algorithms mapped to a standard metric achieve
> clinically useful quantification of average Parkinson's disease symptom severity.
> Journal of the Royal Society Interface, 8(59), 842-855.

> Tsanas, A. (2012). Accurate telemonitoring of Parkinson's disease symptom severity
> using nonlinear speech signal processing and statistical machine learning.
> D.Phil. thesis, University of Oxford.

## License

GPL-3.0

## Authors

- Original MATLAB code: Athanasios Tsanas (2014)
- Python implementation: 2025
- R integration: superassp package team
