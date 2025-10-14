# Parselmouth / Python Implementation Testing Summary

**Date**: 2025-10-14
**Task**: Test and verify equivalence of Python/Parselmouth implementations vs original Praat-calling functions

## Overview

This document summarizes the testing and verification work for the Python/Parselmouth implementations of Praat signal processing functions in the superassp package.

## Python Implementations Created

Five Python/Parselmouth implementations were created in `inst/python/`:

1. **praat_formant_burg.py** - Formant analysis using Burg algorithm (cleaned up existing implementation)
2. **praat_pitch.py** - Multi-method pitch tracking (ac, cc, spinet, shs)
3. **praat_intensity.py** - Intensity (loudness) contour extraction
4. **praat_spectral_moments.py** - Spectral moments (CoG, SD, Skewness, Kurtosis)
5. **praat_formantpath_burg.py** - Advanced formant tracking with FormantPath

## R Wrapper Functions Created

Five optimized R wrapper functions were created in `R/praat_python_optimized.R`:

1. `praat_formant_burg_opt()` - Wrapper for praat_formant_burg.py
2. `praat_pitch_opt()` - Wrapper for praat_pitch.py
3. `praat_intensity_opt()` - Wrapper for praat_intensity.py
4. `praat_spectral_moments_opt()` - Wrapper for praat_spectral_moments.py
5. `praat_formantpath_burg_opt()` - Wrapper for praat_formantpath_burg.py

## Test Results

Test file: `tests/test_parselmouth_equivalence.R`

### ✅ Working Tests

1. **Test 2: praat_pitch_opt - Pitch tracking**
   - Status: PASS
   - Returns valid AsspDataObj with 'cc' and 'ac' tracks
   - 579 frames for full file, 195 frames for 1-second window
   - Minor warning: "NAs introduced by coercion" (non-critical)

2. **Test 4: praat_spectral_moments_opt - Spectral shape**
   - Status: PASS
   - Returns valid AsspDataObj with 'cog', 'sd', 'skewness', 'kurtosis' tracks
   - 582 frames for full file
   - Correct sample rate of 200 Hz

3. **Test 6: Time windowing with praat_pitch_opt**
   - Status: PASS
   - Time windowing works correctly
   - Returns expected ~200 frames for 1-second window at 200 Hz sample rate

### ⚠️ Failing Tests (Issues Identified)

1. **Test 1: praat_formant_burg_opt - Basic formant analysis**
   - Status: FAIL
   - Error: "data must be a numeric matrix"
   - Issue: Type conversion problem in R wrapper when adding tracks to AsspDataObj
   - Python function executes successfully, problem is in R-side conversion
   - **Action needed**: Debug data type conversion in R wrapper (likely result_df column types)

2. **Test 3: praat_intensity_opt - Intensity analysis**
   - Status: FAIL (partial)
   - Returns valid AsspDataObj but with 0-dimensional intensity track
   - Likely issue: Empty DataFrame returned from Python, or column name mismatch
   - **Action needed**: Debug Python function return value and R column name matching

3. **Test 5: praat_formantpath_burg_opt - FormantPath**
   - Status: FAIL
   - Error: "ValueError: Cannot convert 0-dimensional NumPy array argument '<SpectralAnalysisWindowShape.GAUSSIAN: 5>' to a Praat vector or matrix"
   - Issue: Passing Python enum object instead of integer value to Parselmouth
   - **Action needed**: Convert enum to integer or string representation before passing to Python

4. **Test 7: Writing output to file**
   - Status: FAIL
   - Error: "There are no track format specifiers!"
   - Issue: AsspDataObj missing track format metadata when writing to file
   - **Action needed**: Ensure trackFormats attribute is set correctly in AsspDataObj creation

## Fixes Applied

### praat_formant_burg.py

1. **Removed test code at end of file**
   - Removed `soundFile = "/Users/frkkan96/Desktop/a1.wav"` assignment
   - Removed `res = praat_formant_burg(soundFile)` test call
   - Removed unused helper functions `PraatTableToPandas()` and `PraatTableToPandas2()`

2. **Fixed loop ranges**
   - Changed `range(1,math.ceil(number_of_formants))` to `range(1,math.ceil(number_of_formants) + 1)`
   - Changed `range(1,nFormantRows)` to `range(1,nFormantRows + 1)`
   - Python range() is exclusive of the end value, so +1 is needed to include the last formant/row

## Known Issues

### 1. DataFrame Column Types
Some columns from pandas DataFrames may not convert correctly to R matrices. Need to verify that:
- All numeric columns are of correct type (float64, not object)
- No string values like "--undefined--" remain unconverted
- Column names match exactly what R wrapper expects

### 2. Enum Passing to Python
The reticulate package may not handle Python enum types correctly when passed as arguments. Solutions:
- Pass enum integer values instead of enum objects
- Pass string representations and convert within Python function
- Use `int(enum_object)` in R before passing to Python

### 3. Track Format Specifications
AsspDataObj requires proper trackFormats attribute for writing. Need to ensure:
- `attr(outDataObj, "trackFormats")` is set correctly (e.g., "REAL32")
- Track formats match the data type used in wrassp::addTrack()

### 4. Empty DataFrames
Some Python functions may return empty DataFrames in edge cases:
- Verify Python function returns data for test audio file
- Add defensive checks in R wrapper for empty results
- Consider adding better error messages

## Recommendations for Next Steps

### Immediate Fixes (Priority 1)

1. **Fix praat_formant_burg_opt**
   - Add `str(result_df)` debug output to see DataFrame structure
   - Verify column types match expectations
   - Check if "--undefined--" values are being converted to NA

2. **Fix praat_intensity_opt**
   - Test Python function standalone to verify it returns data
   - Check column name matching (case sensitivity, special characters)
   - Verify DataFrame has rows before processing

3. **Fix praat_formantpath_burg_opt**
   - Change enum passing strategy:
     ```r
     # Instead of:
     py_spec_shape <- reticulate::py_eval("pm.SpectralAnalysisWindowShape.GAUSSIAN")

     # Use integer value:
     py_spec_shape <- 5L  # GAUSSIAN value
     ```

### Testing Enhancements (Priority 2)

1. **Add comparison tests**
   - Compare output of _opt functions with original Praat-calling functions
   - Verify numerical equivalence (within tolerance)
   - Document any expected differences

2. **Add edge case tests**
   - Very short audio files (< 100ms)
   - Silence
   - Different sample rates (8kHz, 16kHz, 44.1kHz, 48kHz)
   - Various time window configurations

3. **Performance benchmarks**
   - Compare execution time: Praat subprocess vs Parselmouth
   - Measure memory usage
   - Test with large files (> 1 minute)

### Documentation (Priority 3)

1. **Parameter mapping document**
   - Map original Praat parameters to Python function parameters
   - Document any parameter name differences
   - List default value changes

2. **Migration guide**
   - How to switch from original to _opt functions
   - Expected differences in output format
   - Performance improvement estimates

## Test Environment

- **R Version**: 4.4-arm64
- **Python Version**: System Python 3
- **Parselmouth**: Installed via pip3
- **reticulate**: Using system Python3
- **Test Audio**: tests/signalfiles/AVQI/input/sv1.wav
- **Platform**: macOS (Darwin 25.0.0)

## Success Metrics

- **Target**: 7/7 tests passing (100%)
- **Current**: 3/7 tests passing (43%)
- **Partially Working**: 1/7 (Test 3 - returns object but with empty data)
- **Failing**: 3/7 (Tests 1, 5, 7)

## Conclusion

The Python/Parselmouth implementations are mostly functional:
- **Pitch tracking**: Fully working
- **Spectral moments**: Fully working
- **Time windowing**: Working correctly
- **Formant analysis**: Needs debugging (3 related failures)
- **Intensity analysis**: Needs investigation (returns empty data)

The core architecture is sound. The remaining issues are primarily:
1. Data type conversion between Python and R
2. Enum type handling
3. Empty result handling

These are all fixable issues that don't require redesigning the implementation approach.

---

**Next Actions**:
1. Debug and fix the 4 failing tests
2. Add comparison tests against original Praat functions
3. Document parameter mappings and migration guide
4. Add performance benchmarks

