# Praat-to-Parselmouth Python Implementations Summary

This document summarizes the Python/Parselmouth implementations created for Praat signal processing functions in the superassp package.

## Overview

Based on analysis of `/Users/frkkan96/Documents/src/superassp/R/praat_ssff.R`, the following R functions call Praat scripts for signal processing:

1. `praat_formant_burg` - Formant analysis using Burg algorithm
2. `praat_formantpath_burg` - Formant tracking using FormantPath
3. `praat_sauce` - PraatSauce voice quality measures
4. `praat_intensity` - Intensity (loudness) contour
5. `praat_moments` - Spectral moments analysis
6. `praat_pitch` - Multi-method pitch tracking

## Implementation Status

### ✅ Already Implemented

1. **praat_formant_burg.py**
   - Location: `/Users/frkkan96/Documents/src/superassp/inst/python/praat_formant_burg.py`
   - Praat script: `inst/praat/formant_burg.praat`
   - Status: **ALREADY EXISTS** (skipped)
   - Function: Basic Burg formant analysis

### ✅ Newly Implemented

2. **praat_pitch.py**
   - Location: `/Users/frkkan96/Documents/src/superassp/inst/python/praat_pitch.py`
   - Praat script: `inst/praat/praat_pitch.praat`
   - Status: **CREATED**
   - Function: Multi-method pitch tracking
   - Methods implemented:
     - Cross-correlation (cc)
     - Autocorrelation (ac)
     - SPINET (optional)
     - Subharmonic-to-harmonic ratio (SHS, optional)
   - Key features:
     - Supports time windowing with `extract_part()`
     - Returns DataFrame with pitch tracks from multiple algorithms
     - Boolean flag `only_correlation_methods` to enable/disable SPINET and SHS

3. **praat_intensity.py**
   - Location: `/Users/frkkan96/Documents/src/superassp/inst/python/praat_intensity.py`
   - Praat script: `inst/praat/intensity.praat`
   - Status: **CREATED**
   - Function: Intensity (loudness) contour extraction
   - Key features:
     - Computes dB intensity over time
     - Supports mean subtraction for normalization
     - Configurable minimum F0 for time resolution
     - Returns DataFrame with Time(s) and Intensity(dB) columns

4. **praat_spectral_moments.py**
   - Location: `/Users/frkkan96/Documents/src/superassp/inst/python/praat_spectral_moments.py`
   - Praat script: `inst/praat/praat_spectral_moments.praat`
   - Status: **CREATED**
   - Function: Spectral shape analysis using moments 1-4
   - Measures computed:
     - Center of Gravity (Moment 1): Mean frequency
     - Standard Deviation (Moment 2): Frequency spread
     - Skewness (Moment 3): Spectral asymmetry
     - Kurtosis (Moment 4): Spectral peakedness
   - Key features:
     - Frame-by-frame spectral analysis
     - Configurable window length and frequency resolution
     - Automatic Nyquist frequency detection
     - Returns DataFrame with Time, CenterOfGravity, SD, Skewness, Kurtosis

5. **praat_formantpath_burg.py**
   - Location: `/Users/frkkan96/Documents/src/superassp/inst/python/praat_formantpath_burg.py`
   - Praat script: `inst/praat/formantpath_burg.praat`
   - Status: **CREATED**
   - Function: Advanced formant tracking with automatic ceiling optimization
   - Key features:
     - Uses FormantPath for robust formant tracking
     - Automatically explores multiple formant ceilings
     - Optional formant tracking with cost functions
     - Extracts formant frequencies (F1-F5), bandwidths (B1-B5), and amplitudes (L1-L5)
     - Computes formant amplitudes from spectrogram
     - Returns comprehensive DataFrame with all formant parameters

### ⚠️ Not Implemented

6. **praat_sauce** (praatsauce.praat)
   - Status: **NOT IMPLEMENTED**
   - Reason: This is a complex voice quality analysis script that:
     - Requires multiple dependent Praat scripts
     - Includes formant analysis, harmonics-to-noise ratio, spectral tilt, etc.
     - Would require significant additional work to implement
   - Recommendation: Could be implemented in a future update if needed

## Implementation Patterns

All implementations follow a consistent pattern:

### Standard Function Signature
```python
def function_name(
    soundFile,
    beginTime=0.0,
    endTime=0.0,
    # ... specific parameters ...
    windowShape=pm.WindowShape.GAUSSIAN1,
    relativeWidth=1.0):
```

### Common Features

1. **Time Windowing**
   ```python
   snd = pm.Sound(soundFile)
   dur = snd.get_total_duration()
   if beginTime > 0.0 or endTime > 0.0:
       if beginTime >= 0.0 and endTime <= dur:
           snd = snd.extract_part(beginTime, endTime, windowShape, relativeWidth, True)
   ```

2. **Praat Object Manipulation**
   - Use `pm.praat.call()` for all Praat commands
   - Follow Praat's exact parameter order and naming

3. **Table Conversion**
   ```python
   result = pd.read_table(io.StringIO(pm.praat.call(table, "List", True)))
   return result
   ```

4. **Error Handling**
   - Undefined values from Praat appear as `"--undefined--"` strings
   - These are automatically converted to NaN by pandas

## Key Implementation Notes

### praat_pitch.py
- Implements four different pitch tracking algorithms
- The `only_correlation_methods` parameter controls whether SPINET and SHS are computed
- Default is True (only cc and ac) for faster computation
- Matrix merging follows Praat's "append rows" then "transpose" pattern

### praat_intensity.py
- Simple and direct implementation
- The `time_step=0.0` parameter triggers automatic calculation in Praat
- Uses IntensityTier → TableOfReal → Table conversion pipeline

### praat_spectral_moments.py
- Most complex implementation due to frame-by-frame processing
- Creates and populates a Praat Table object iteratively
- Each frame requires creating a spectrum slice, computing moments, and cleanup
- Automatic Nyquist frequency detection when `maximum_frequency=0.0`

### praat_formantpath_burg.py
- Uses the newer FormantPath algorithm (more robust than basic Burg)
- Formant tracking is optional but recommended for clean tracks
- Amplitude extraction requires a separate spectrogram analysis
- Handles undefined formant values gracefully

## Testing Recommendations

For each implementation, test with:

1. **Basic usage**: Default parameters with a simple speech file
2. **Time windowing**: Non-zero beginTime and endTime
3. **Parameter variations**: Different window sizes, frequency ranges, etc.
4. **Edge cases**: Very short files, silence, extreme parameter values
5. **Output validation**: Compare with original Praat script outputs

## Dependencies

All implementations require:
- `parselmouth` - Praat integration
- `pandas` - DataFrame output
- `numpy` - Numerical operations
- `io` - String I/O for table conversion
- `math` - Mathematical functions (formantpath only)

## Integration with R Package

These Python functions are called from R via the `reticulate` package. The R wrapper functions in `praat_ssff.R`:

1. Load the Python module
2. Convert R parameters to Python types
3. Call the Python function
4. Convert DataFrame results back to R data structures
5. Optionally save as SSFF files

## Future Work

Potential enhancements:

1. **Implement praat_sauce**: Full voice quality analysis suite
2. **Add more pitch algorithms**: RAPT, SWIPE, etc. (some may already exist)
3. **Optimize performance**: Vectorized operations where possible
4. **Add visualization**: Optional plot generation
5. **Enhanced error handling**: More descriptive error messages
6. **Batch processing**: Multi-file processing optimization

## File Locations

### Python Implementations
- `/Users/frkkan96/Documents/src/superassp/inst/python/praat_pitch.py`
- `/Users/frkkan96/Documents/src/superassp/inst/python/praat_intensity.py`
- `/Users/frkkan96/Documents/src/superassp/inst/python/praat_spectral_moments.py`
- `/Users/frkkan96/Documents/src/superassp/inst/python/praat_formantpath_burg.py`
- `/Users/frkkan96/Documents/src/superassp/inst/python/praat_formant_burg.py` (existing)

### R Wrapper Functions
- `/Users/frkkan96/Documents/src/superassp/R/praat_ssff.R`

### Praat Scripts (Reference)
- `/Users/frkkan96/Documents/src/superassp/inst/praat/praat_pitch.praat`
- `/Users/frkkan96/Documents/src/superassp/inst/praat/intensity.praat`
- `/Users/frkkan96/Documents/src/superassp/inst/praat/praat_spectral_moments.praat`
- `/Users/frkkan96/Documents/src/superassp/inst/praat/formantpath_burg.praat`
- `/Users/frkkan96/Documents/src/superassp/inst/praat/formant_burg.praat`

---

**Document Created**: 2025-10-14
**Author**: Claude Code Analysis
**Package**: superassp
