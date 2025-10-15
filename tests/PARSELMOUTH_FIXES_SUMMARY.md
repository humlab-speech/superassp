# Python/Parselmouth Implementation - Fixes Applied

## Test Results Summary

**Before fixes:** 3/7 tests passing (43%)  
**After fixes:** 7/7 tests passing (100%) ✓

All Python/Parselmouth-based optimized functions are now working correctly!

---

## Issues Fixed

### Issue 1: Data Type Conversion Error (praat_formant_burg_opt)
**Error:** "data must be a numeric matrix"

**Root Cause:** 
- pandas DataFrames returned columns F4(Hz), B4(Hz), F5(Hz), B5(Hz), L4(dB), L5(dB) as character type instead of numeric
- This occurred when pandas encountered "--undefined--" or "?" values in the data
- R's as.matrix() failed when trying to convert character vectors to numeric matrices

**Fix Applied:**
- Added `as.numeric()` conversion before `as.matrix()` in all track-adding loops
- Applied to both praat_formant_burg_opt and praat_formantpath_burg_opt
- NAs are replaced with 0 after conversion

**Files Modified:**
- `R/praat_python_optimized.R` lines 189-225 (praat_formant_burg_opt)
- `R/praat_python_optimized.R` lines 985-1021 (praat_formantpath_burg_opt)

---

### Issue 2: Empty Data (praat_intensity_opt)
**Error:** Intensity track returned 0 dimensions

**Root Cause:**
- Python function returned column names with spaces: "Time (s)" and "Intensity (dB)"
- R code used backtick notation `result_df$\`Intensity(dB)\`` looking for column without spaces
- Column access failed silently, returning NULL
- Result: empty intensity track with 0 dimensions

**Fix Applied:**
- Changed from backtick `$` operator to double bracket `[[]]` notation
- Updated column names to match Python output exactly:
  - `result_df$\`Time(s)\`` → `result_df[["Time (s)"]]`
  - `result_df$\`Intensity(dB)\`` → `result_df[["Intensity (dB)"]]`

**Files Modified:**
- `R/praat_python_optimized.R` lines 561-584 (praat_intensity_opt)

---

### Issue 3: Enum Type Error (praat_formantpath_burg_opt)
**Error:** "ValueError: Cannot convert 0-dimensional NumPy array argument '<SpectralAnalysisWindowShape.GAUSSIAN: 5>'"

**Root Cause:**
- R was passing Python enum object through reticulate to Python function
- When Python function tried to pass this enum to `pm.praat.call()`, it didn't convert properly
- Parselmouth couldn't handle the 0-dimensional NumPy array representation of the enum

**Fix Applied:**
1. Changed Python function signature to accept string instead of enum:
   ```python
   spectrogram_window_shape="Gaussian"  # String instead of enum
   ```

2. Added string handling in Python function:
   ```python
   spec_window = spectrogram_window_shape if isinstance(spectrogram_window_shape, str) else "Gaussian"
   spectrogram = pm.praat.call(snd, "To Spectrogram", ..., spec_window)
   ```

3. Updated R wrapper to pass string directly instead of enum object

**Files Modified:**
- `inst/python/praat_formantpath_burg.py` lines 31, 145-147
- `R/praat_python_optimized.R` line 369 (removed enum conversion)

---

### Issue 4: Missing Track Format Specifiers (File Writing)
**Error:** "There are no track format specifiers!" → "Not enough format specifiers for the data tracks"

**Root Cause:**
- Discovered bug in wrassp::addTrack() function:
  ```r
  else append(attr(dobj, "trackFormats"), format)  # Bug: result not assigned!
  ```
- The append() result wasn't being assigned back to the attribute
- trackFormats stayed empty even after adding tracks
- wrassp::write.AsspDataObj() requires trackFormats to specify data type for each track

**Fix Applied:**
- Manually append to trackFormats after each wrassp::addTrack() call:
  ```r
  outDataObj <- wrassp::addTrack(outDataObj, "intensity", as.matrix(intensity_data), "REAL32")
  # Manually fix trackFormats due to wrassp::addTrack bug
  attr(outDataObj, "trackFormats") <- c(attr(outDataObj, "trackFormats"), "REAL32")
  ```
- Applied to all 5 optimized functions after every addTrack() call

**Files Modified:**
- `R/praat_python_optimized.R`:
  - praat_formant_burg_opt: lines 195, 210, 225
  - praat_pitch_opt: lines 420, 430, 440, 449
  - praat_intensity_opt: line 589
  - praat_spectral_moments_opt: lines 757, 766, 775, 784
  - praat_formantpath_burg_opt: lines 991, 1006, 1021

---

## Technical Details

### pandas Mixed Type Columns
When pandas DataFrames contain mixed types (numeric and string like "--undefined--"), pandas converts the entire column to character type to maintain consistency. This requires explicit `as.numeric()` conversion in R.

### R Column Access Methods
- **Backtick operator `$`**: Doesn't work reliably with column names containing spaces
- **Double bracket `[[]]`**: Works correctly with any column name, including those with spaces

### reticulate and Python Enums
Python enum objects don't convert well through reticulate when passed to Parselmouth's `pm.praat.call()`. Using string representations and converting inside Python functions provides better compatibility.

### wrassp::addTrack() Bug
The function has a bug where `append()` is called without assignment:
```r
else append(attr(dobj, "trackFormats"), format)  # Bug!
```
Should be:
```r
else attr(dobj, "trackFormats") <- append(attr(dobj, "trackFormats"), format)
```
Until this is fixed in wrassp, manual trackFormats management is required.

---

## Verification

All tests now pass successfully:

1. ✓ praat_formant_burg_opt - Basic formant analysis (574 frames, 15 tracks)
2. ✓ praat_pitch_opt - Pitch tracking (579 frames, 2 tracks)
3. ✓ praat_intensity_opt - Intensity analysis (175 frames, 1 track)
4. ✓ praat_spectral_moments_opt - Spectral shape (582 frames, 4 tracks)
5. ✓ praat_formantpath_burg_opt - FormantPath (574 frames, 9 tracks)
6. ✓ Time windowing test (195 frames as expected)
7. ✓ File writing test (844 bytes written successfully)

---

## Notes

### Warnings
Several tests show "NAs introduced by coercion" warnings. These are expected and harmless - they occur when converting strings like "--undefined--" or "?" to numeric, which produces NA values that are immediately replaced with 0.

### Performance
The Python/Parselmouth implementations provide significant performance improvements over calling external Praat processes while maintaining output compatibility.

---

## Files Changed Summary

- `R/praat_python_optimized.R` - All 5 _opt functions updated
- `inst/python/praat_formantpath_burg.py` - Changed spectrogram_window_shape to string
- Debug scripts created (can be removed):
  - `tests/debug_formant_burg.R`
  - `tests/debug_intensity.R`
  - `tests/debug_file_writing.R`

---

Date: 2025-10-14
Status: All issues resolved ✓
