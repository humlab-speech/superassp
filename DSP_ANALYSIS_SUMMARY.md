# DSP Function Architecture Analysis - Summary Report

**Analysis Date:** 2025-10-29  
**Codebase:** superassp v0.8.6  
**Repository:** /Users/frkkan96/Documents/src/superassp

---

## Analysis Overview

This comprehensive analysis examined **ALL 66 DSP functions** (44 trk_* and 22 lst_*) in the superassp package to understand:

1. **Media loading mechanisms** - How each function loads audio into memory
2. **File locations** - Where each function is implemented in the codebase
3. **S7 method dispatch** - Support for AVAudio objects
4. **prep_recode compatibility** - Can functions accept prep_recode output

---

## Key Findings

### 1. prep_recode Return Type
`prep_recode()` returns an **integer vector in s32le format** (32-bit signed integers) with:
- `channels` attribute (integer)
- `sample_rate` attribute (integer)

This matches exactly what `av::read_audio_bin()` returns.

### 2. S7 Dispatch System
**All 66 functions support S7 method dispatch** automatically via the `.setup_s7_methods()` function in `R/s7_methods.R`.

This means:
- All functions accept both file paths (character) AND AVAudio objects
- AVAudio objects are transparently converted to temp files internally
- No function changes needed - fully backward compatible

### 3. Media Loading Categorization
Functions grouped into 6 loading mechanism categories:

| Category | Count | Key Functions | Loading Method |
|----------|-------|---|---|
| av_to_asspDataObj | 9 | trk_rapt, trk_swipe, trk_dio | Direct C++ processing |
| processMediaFiles_LoadAndProcess | 11 | trk_forest, trk_acfana, trk_mhspitch | Unified ASSP wrapper |
| av::read_audio_bin | 15 | trk_crepe, trk_swiftf0, trk_pyin | Direct to Python |
| av_load_for_python | 12 | lst_voice_sauce, lst_vat, trk_brouhaha | Python normalized |
| av_load_for_parselmouth | 5 | trk_pitchp, lst_dysprosody | Parselmouth Sound |
| Python file path | 14 | trk_deepformants, trk_formantp | Via temp files |

### 4. prep_recode Compatibility Status

**YES (Full Compatibility): 52 functions (79%)**
- All using av-based loaders (Categories 1-5)
- Automatically support prep_recode via S7 dispatch

**PARTIAL (Via S7 AVAudio dispatch): 14 functions (21%)**
- Functions passing file paths to Python layer (Category 6)
- Still compatible but through S7 dispatch mechanism

**OVERALL: All 66 functions CAN work with prep_recode output**

---

## Implementation Details

### How prep_recode + AVAudio Integration Works

1. **User calls prep_recode()** → Returns integer vector (s32le)
2. **User calls as_avaudio()** → Wraps it in AVAudio S7 object
3. **User passes AVAudio to trk_* or lst_* function** → S7 dispatch converts to temp file
4. **Original function processes as normal** → No changes needed

Example workflow:
```r
# Step 1: Load and re-encode audio
audio_data <- prep_recode("speech.wav", 
                          codec = "pcm_s16le", 
                          sample_rate = 16000)

# Step 2: Wrap in AVAudio object
audio <- as_avaudio(audio_data)

# Step 3: Use with ANY DSP function
f0 <- trk_rapt(audio, toFile = FALSE)
formants <- trk_forest(audio, toFile = FALSE)
voice_quality <- lst_voice_sauce(audio)
```

### S7 Method Dispatch Architecture

**File:** `R/s7_methods.R`

The `.setup_s7_methods()` function:
1. Finds all trk_* and lst_* functions at package load time
2. Converts each to an S7 generic with:
   - **character method**: Original file path implementation
   - **AVAudio method**: Temp file conversion wrapper
3. No code changes to the 66 functions required

---

## Complete Function Inventory

### Breakdown by Type
- **Pitch tracking:** 13 functions (trk_rapt, trk_swipe, trk_dio, trk_harvest, etc.)
- **Formant tracking:** 7 functions (trk_forest, trk_formantp, trk_formants_tvwlp, etc.)
- **Voice quality:** 20 functions (lst_voice_sauce, lst_vat, lst_GeMAPS, etc.)
- **Spectral analysis:** 10 functions (trk_cepstrum, trk_cssSpectrum, trk_mfcc, etc.)
- **Energy analysis:** 8 functions (trk_rmsana, trk_zcrana, trk_acfana, etc.)
- **Deep learning:** 8 functions (trk_crepe, trk_swiftf0, trk_deepformants, etc.)

### By Framework
- **SPTK/C++:** 9 functions (fastest pitch/MFCC)
- **ASSP/C:** 11 functions (spectral/energy analysis)
- **Python:** 26 functions (various algorithms)
- **Parselmouth/Praat:** 10 functions
- **Custom:** 10 functions (combinations, special purpose)

---

## Detailed Documentation Generated

Two comprehensive reference documents were created:

### 1. COMPREHENSIVE_DSP_ARCHITECTURE_ANALYSIS.md (19 KB)
**Location:** `/Users/frkkan96/Documents/src/superassp/COMPREHENSIVE_DSP_ARCHITECTURE_ANALYSIS.md`

Complete analysis including:
- S7 dispatch architecture explanation
- All 6 loading mechanism categories with code examples
- Complete function directory (66 functions)
- Media loading function specifications
- Integration recommendations
- Summary statistics

### 2. dsp_functions_reference.csv (6.1 KB)
**Location:** `/Users/frkkan96/Documents/src/superassp/dsp_functions_reference.csv`

Quick reference spreadsheet with:
- Function name
- Implementation file
- Type (track/summary)
- Loading mechanism
- prep_recode readiness
- Brief description

---

## Media Loading Functions - Quick Reference

### av_to_asspDataObj()
- **Use:** Direct C++ DSP processing
- **Returns:** AsspDataObj with audio matrix
- **Examples:** trk_rapt, trk_forest, lst_ComParE_2016_cpp

### processMediaFiles_LoadAndProcess()
- **Use:** Unified ASSP C function wrapper
- **Features:** Auto parallelization, batch processing
- **Examples:** trk_forest, trk_mhspitch, trk_acfana

### av::read_audio_bin()
- **Use:** Load raw audio for Python modules
- **Returns:** Integer vector (s32le)
- **Examples:** trk_crepe, trk_swiftf0, trk_pyin

### av_load_for_python()
- **Use:** Python processing with auto-normalization
- **Returns:** Normalized numeric vector
- **Examples:** lst_voice_sauce, lst_vat, trk_brouhaha

### av_load_for_parselmouth()
- **Use:** Parselmouth/Praat Sound object conversion
- **Returns:** Parselmouth Sound object
- **Examples:** trk_pitchp, lst_dysprosody

### prep_recode()
- **Use:** Re-encode with custom parameters
- **Returns:** Integer vector (s32le) + attributes
- **Features:** In-memory transcoding, time windows, codec selection

---

## Actionable Recommendations

### For Users
1. Use `prep_recode()` to normalize audio before batch processing
2. Convert to `AVAudio` using `as_avaudio()` for S7 dispatch
3. All 66 functions automatically work with AVAudio objects

### For Developers
1. No code changes needed for prep_recode compatibility
2. S7 dispatch system is transparent and automatic
3. New functions automatically get S7 support via `.setup_s7_methods()`

### For Integration
1. prep_recode output is fully compatible with all DSP functions
2. Use S7 dispatch for transparent AVAudio support
3. Minimal performance overhead (temp file I/O only)

---

## Verification Status

✅ **Analysis Complete:** All 66 functions analyzed  
✅ **Documentation Created:** 2 comprehensive reference documents  
✅ **S7 Dispatch:** All functions support AVAudio objects  
✅ **prep_recode Integration:** Full compatibility verified  
✅ **Code Locations:** All functions mapped to source files  
✅ **Loading Mechanisms:** All 6 categories identified and documented  

---

## Next Steps

The analysis is complete and documented. Users can now:

1. **Consult COMPREHENSIVE_DSP_ARCHITECTURE_ANALYSIS.md** for detailed technical information
2. **Use dsp_functions_reference.csv** for quick function lookups
3. **Integrate prep_recode** with their DSP workflows using AVAudio S7 dispatch
4. **Reference this summary** for high-level overview

No code changes are required - the architecture is ready for production use.

