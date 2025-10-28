# Voxit: Voice and Articulation Complexity Measures

Python implementation of the Voxit toolbox for prosodic complexity analysis.

## Features

Voxit computes 11 prosodic and rhythmic features from speech audio:

### Temporal Features
- **WPM**: Words per minute (speaking rate)
- **pause_count**: Number of pauses (100-3000ms)
- **long_pause_count**: Number of pauses > 3s
- **average_pause_length**: Mean pause duration (seconds)
- **average_pause_rate**: Pauses per second
- **rhythmic_complexity_of_pauses**: Normalized Lempel-Ziv complexity (%)

### Pitch Features
- **average_pitch**: Mean F0 (Hz)
- **pitch_range**: F0 range (octaves)
- **pitch_speed**: F0 velocity (octaves/second)
- **pitch_acceleration**: F0 acceleration (octaves/second²)
- **pitch_entropy**: F0 distribution entropy (bits)

## Installation

### Basic Installation

```bash
pip install numpy scipy lempel_ziv_complexity
```

### With Numba (2-3x speedup)

```bash
pip install numba
```

No compilation required - instant speedup through JIT compilation.

### With Cython (3-5x speedup)

Requires C compiler (gcc/clang/MSVC):

```bash
pip install Cython
python setup_cython.py build_ext --inplace
```

Maximum performance through compiled C extensions.

## Usage

### Python

```python
from voxit import compute_features

# Prepare input data
gentle_data = [
    {'word': 'hello', 'case': 'success', 'start': 0.0, 'end': 0.5},
    {'word': 'world', 'case': 'success', 'start': 0.7, 'end': 1.2},
]

pitch_data = [
    {'time': 0.0, 'frequency': 120.0},
    {'time': 0.01, 'frequency': 121.5},
    # ... more pitch points
]

# Compute features
features = compute_features(
    gentle_data=gentle_data,
    pitch_data=pitch_data,
    start_time=None,  # Optional time windowing
    end_time=None
)

print(features['WPM'])  # Speaking rate
print(features['average_pitch'])  # Mean F0
print(features['pitch_range'])  # F0 range in octaves
```

### R (via superassp)

```r
library(superassp)

# Install voxit
install_voxit(install_numba = TRUE)

# Extract features
features <- lst_voxit(
  "audio.wav",
  alignmentFiles = "alignments.csv"
)

print(features$WPM)
print(features$average_pitch)
```

## Input Formats

### Word Alignments (gentle_data)

List of dictionaries with:
- `word`: Word text
- `case`: Category (e.g., "success", "[noise]")
- `start`: Start time in seconds
- `end`: End time in seconds

### Pitch Track (pitch_data)

List of dictionaries with:
- `time`: Time in seconds
- `frequency`: F0 in Hz (0 for unvoiced)

## Performance

### Benchmarks (on typical 5-second audio)

- **Standard Python**: ~200ms
- **With Numba**: ~80ms (2.5x faster)
- **With Cython**: ~60ms (3.3x faster)

### Optimization Tips

1. **Use Numba**: Fastest without compilation
   ```python
   features = compute_features(..., use_numba=True)
   ```

2. **Use Cython**: Maximum performance
   ```python
   features = compute_features(..., use_cython=True)
   ```

3. **Automatic selection**: Picks best available
   ```python
   features = compute_features(...)  # Auto-detects optimizations
   ```

## Algorithm Details

### Rhythmic Complexity
- Samples speech/pause pattern at 100 Hz
- Binary sequence: 1 (voiced) / 0 (pause, 100-3000ms)
- Normalized Lempel-Ziv complexity

### Pitch Dynamics
- Pitch converted to log₂ scale (octaves)
- Savitzky-Golay smoothing (order=2, window=7) before derivatives
- Velocity: First derivative of smoothed pitch contour
- Acceleration: Second derivative of smoothed pitch contour
- Signed directionless statistics: `mean(abs(x)) * sign(mean(x))`

### Pitch Entropy
- 25-bin histogram over ±1 octave from mean
- Shannon entropy: `-Σ p(i) * log₂(p(i))`

## Dependencies

- **numpy**: Numerical operations
- **scipy**: Savitzky-Golay filter (`scipy.signal.savgol_filter`)
- **lempel_ziv_complexity**: Rhythmic complexity calculation

Optional:
- **numba**: JIT compilation for 2-3x speedup
- **Cython**: C extensions for 3-5x speedup

## References

- Original MATLAB implementation: Voxit toolbox
- Python reimplementation with optimizations (2024)
- SAcC pitch tracker: Ellis & Weiss (2010)

## License

See LICENSE file for details.

## Contributing

Contributions welcome! Please ensure:
- Code follows PEP 8 style
- All tests pass
- Optimized versions match reference implementation
- Documentation is updated

## Testing

```bash
python -m pytest tests/
```

## Version History

- **1.0.0**: Initial release
  - Core functionality
  - Numba optimization
  - Cython compilation support
  - Full compatibility with MATLAB version
