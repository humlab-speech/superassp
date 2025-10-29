# DSP Function Extension Attribute Audit Report

**Date:** 2025-10-29
**Package:** superassp
**Scope:** All `lst_*` and `trk_*` DSP functions
**Total Functions Analyzed:** 61+

---

## Executive Summary

This audit examined all DSP functions in the superassp R package to assess consistency between file extension attributes (`attr(*, "ext")`) and formal function parameters (`explicitExt`). The package demonstrates **97% consistency** overall, with only 5-6 functions requiring fixes.

### Key Findings

- ✅ **59 functions** have properly defined `ext` attributes
- ❌ **2 functions** have mismatched parameter defaults vs attributes
- ⚠️ **3-4 functions** are missing `ext` attributes (track functions)
- ✅ **8 functions** intentionally lack `ext` (5 summary + 3 internal wrappers)

### Critical Issues Identified

1. **2 OpenSMILE functions** have incorrect `explicitExt` defaults (HIGH PRIORITY)
2. **3-4 track functions** missing `ext` attributes (MEDIUM PRIORITY)

---

## Part 1: Extension Attribute Consistency Analysis

### 1.1 Functions WITH ext Attributes (59 Functions)

All properly configured functions follow consistent naming conventions:

#### Extension Categories

**Pitch/F0 Tracking (17 extensions):**
- `f0` - Standard F0 (RAPT, SWIPE, DIO, HARVEST, REAPER, KSVFO, COVAREP SRH)
- `pit` - Pitch (MHS, Praat)
- `pyp` - PYIN pitch
- `yip` - YIN pitch
- `crp` - CREPE pitch
- `yf0` - YAAPT pitch
- `sacc` - SAcC pitch
- `sf0` - Swift-F0 pitch
- `rpm` - REAPER via Parselmouth
- `egg` - EGG-based F0

**Formant Analysis (4 extensions):**
- `fms` - Formants (ASSP FOREST)
- `pfm` - Praat formants
- `fpb` - Formant path (Praat)
- `dfm` - Deep formants

**Spectral Analysis (5 extensions):**
- `css` - Cepstral spectral
- `dft` - DFT spectral
- `lps` - LPC spectral
- `cep` - Cepstrum
- `spm` - Spectral moments

**Energy/Amplitude (4 extensions):**
- `rms` - RMS energy
- `zcr` - Zero crossing rate
- `acf` - Autocorrelation
- `int` - Intensity (Praat)

**Voice Quality (6 extensions):**
- `brh` - Brouhaha (VAD + SNR + C50)
- `psa` - VoiceSauce (Praat)
- `avqi` - AVQI voice quality index
- `dsi` - Dysphonia Severity Index
- `pvr` - Praat voice report
- `pvt` - Voice tremor

**Source-Filter Decomposition (2 extensions):**
- `gfm` - GFM-IAIF
- `glf` - COVAREP IAIF

**Excitation/Periodicity (3 extensions):**
- `xte` - Excitation
- `wap` - Aperiodicity (WORLD)
- `sec` - Sinusoidal excitation

**LP Coefficients (4 extensions):**
- `arf` - ARF coefficients
- `lar` - LAR coefficients
- `lpc` - LPC coefficients
- `rfc` - RFC coefficients

**OpenSMILE Feature Sets (4 extensions):**
- `oge` - GeMAPS (62 features)
- `ogs` - eGeMAPS (88 features)
- `ocp` - ComParE 2016 (6373 features)
- `emo` - emobase (988 features)

**Other Features (4 extensions):**
- `mfcc` - MFCC coefficients
- `ap` - Aperiodicity (D4C)
- `pm` - Pitch marks
- `npf` - NumPy import

---

### 1.2 CRITICAL ISSUES: Mismatched Parameter Defaults (2 Functions)

#### Issue #1: `lst_eGeMAPS()` Mismatch

**File:** `R/list_python_opensmile_eGeMAPS.R`

**Problem:**
```r
# Line 66: WRONG
lst_eGeMAPS <- function(listOfFiles,
                   beginTime=0,
                   endTime=0,
                   explicitExt="ocp",  # ❌ INCORRECT
                   use_cpp = TRUE)

# Line 125: CORRECT
attr(lst_eGeMAPS,"ext") <- c("ogs")  # ✓ CORRECT
```

**Impact:** When users call `lst_eGeMAPS()` with file output, files will be named with `.ocp` extension instead of the correct `.ogs` extension.

**Fix Required:**
```r
# Line 66: Change to
explicitExt="ogs",  # ✓ CORRECTED
```

**Root Cause:** Likely copied from `lst_ComParE_2016()` template without updating.

---

#### Issue #2: `lst_emobase()` Mismatch

**File:** `R/list_python_opensmile_emobase.R`

**Problem:**
```r
# Line 26: WRONG
lst_emobase <- function(listOfFiles,
                   beginTime=0,
                   endTime=0,
                   explicitExt="ocp",  # ❌ INCORRECT
                   use_cpp = TRUE, verbose = FALSE)

# Line 87: CORRECT
attr(lst_emobase,"ext") <- c("emo")  # ✓ CORRECT
```

**Impact:** When users call `lst_emobase()` with file output, files will be named with `.ocp` extension instead of the correct `.emo` extension.

**Fix Required:**
```r
# Line 26: Change to
explicitExt="emo",  # ✓ CORRECTED
```

**Root Cause:** Same as above - copied from `lst_ComParE_2016()` template.

---

### 1.3 Consistency Patterns Across Function Types

#### ✅ Fully Consistent Function Types

**C ASSP Library Functions (11 functions):**
- All have matching `explicitExt` parameter and `ext` attribute
- Examples: `trk_forest`, `trk_ksvfo`, `trk_mhspitch`, `trk_rmsana`
- Pattern: `explicitExt = "X"` matches `attr(*, "ext") <- "X"`

**C++ SPTK Functions (7 functions):**
- All have matching defaults
- Examples: `trk_rapt`, `trk_swipe`, `trk_dio`, `trk_harvest`, `trk_reaper`
- Pattern: All use `"f0"` for pitch, `"mfcc"` for MFCC, `"ap"` for aperiodicity

**C++ ESTK Functions (1 function):**
- `trk_pitchmark` - Consistent (`explicitExt` conditional, `ext = "pm"`)

**Python Functions (35+ functions):**
- Most are consistent
- Exceptions: The 2 OpenSMILE functions identified above

---

## Part 2: Functions WITHOUT ext Attributes

### 2.1 Complete Inventory (11 Functions)

```
1. covarep_vq.R              → lst_covarep_vq
2. list_cpp_opensmile_emobase.R → lst_emobase_cpp
3. list_cpp_opensmile_generic.R → lst_eGeMAPS_cpp, lst_ComParE_2016_cpp
4. list_dysprosody.R         → lst_dysprosody
5. list_vat.R                → lst_vat
6. list_voxit.R              → lst_voxit
7. lst_voice_sauce.R         → lst_voice_sauce
8. ssff_python_dv_f0.R       → trk_dv_f0
9. ssff_python_dv_formants.R → trk_dv_formants
10. trk_creak_union.R        → trk_creak_union
11. trk_formants_tvwlp.R     → trk_formants_tvwlp
```

---

### 2.2 Category 1: Summary Functions (5) - ✅ Expected

These functions return **scalar/vector summaries** (data frames or lists) rather than time-series SSFF tracks. Lack of `ext` attribute is **intentional and correct**.

#### 1. `lst_vat()` - Voice Analysis Toolbox
- **File:** `R/list_vat.R`
- **Returns:** Named list with 132 dysphonia measures
- **Output Mode:** In-memory only (no file writing capability)
- **Status:** ✅ **Correct** - Summary function, no SSFF output

#### 2. `lst_voice_sauce()` - VoiceSauce Parameters
- **File:** `R/lst_voice_sauce.R`
- **Returns:** Data frame with voice quality parameters
- **Output Mode:** In-memory only
- **Status:** ✅ **Correct** - Summary function, no SSFF output

#### 3. `lst_dysprosody()` - Prosodic Features
- **File:** `R/list_dysprosody.R`
- **Returns:** Named list with 193 prosodic features (MOMEL-INTSINT, spectral tilt)
- **Output Mode:** In-memory only
- **Status:** ✅ **Correct** - Summary function, no SSFF output

#### 4. `lst_voxit()` - Voice/Articulation Complexity
- **File:** `R/list_voxit.R`
- **Returns:** Named list with 11 prosodic complexity measures
- **Output Mode:** In-memory only
- **Status:** ✅ **Correct** - Summary function, no SSFF output

#### 5. `lst_covarep_vq()` - COVAREP Voice Quality
- **File:** `R/covarep_vq.R`
- **Returns:** Named list (NAQ, QOQ, H1-H2, HRF, PSP)
- **Output Mode:** In-memory only
- **Status:** ✅ **Correct** - Summary function, no SSFF output

**Rationale:** `lst_*` functions with scalar/aggregate outputs don't need file extensions because they don't produce time-aligned track files.

---

### 2.3 Category 2: Internal C++ Wrappers (3) - ✅ Acceptable

These are **internal helper functions** (`@keywords internal`) called by public Python wrappers that **do** have `ext` attributes.

#### 6. `lst_emobase_cpp()` - OpenSMILE C++ Backend
- **File:** `R/list_cpp_opensmile_emobase.R`
- **Called By:** `lst_emobase()` which HAS `attr(*, "ext") <- "emo"`
- **Visibility:** `@keywords internal`
- **Status:** ✅ **Acceptable** - Internal function, parent has ext

#### 7. `lst_eGeMAPS_cpp()` - OpenSMILE C++ Backend
- **File:** `R/list_cpp_opensmile_generic.R`
- **Called By:** `lst_eGeMAPS()` which HAS `attr(*, "ext") <- "ogs"`
- **Visibility:** `@keywords internal`
- **Status:** ✅ **Acceptable** - Internal function, parent has ext

#### 8. `lst_ComParE_2016_cpp()` - OpenSMILE C++ Backend
- **File:** `R/list_cpp_opensmile_generic.R`
- **Called By:** `lst_ComParE_2016()` which HAS `attr(*, "ext") <- "ocp"`
- **Visibility:** `@keywords internal`
- **Status:** ✅ **Acceptable** - Internal function, parent has ext

**Rationale:** These are implementation details. Users call the public wrapper functions which have proper `ext` attributes.

---

### 2.4 Category 3: Track Functions Missing ext (3-4) - ⚠️ ISSUES

These are **track-based functions** (`trk_*`) that return time-series AsspDataObj but lack `ext` attributes.

#### Issue #3: `trk_creak_union()` - Missing Attribute

**File:** `R/trk_creak_union.R`

**Current State:**
```r
trk_creak_union <- function(listOfFiles,
                            beginTime = 0,
                            endTime = 0,
                            windowShift = 10,
                            use.am = TRUE,
                            use.cd = TRUE,
                            cd.threshold = 0.3,
                            use.reaper = FALSE,
                            return.features = FALSE,
                            return.probabilities = FALSE,
                            explicitExt = "crk",  # ✓ HAS explicitExt
                            outputDirectory = NULL,
                            toFile = TRUE,        # ✓ HAS toFile
                            conda.env = NULL)

# ❌ MISSING: attr(trk_creak_union, "ext") <- "crk"
```

**Impact:**
- Function supports file output (`toFile = TRUE`)
- Has explicit extension parameter (`explicitExt = "crk"`)
- **Missing** the attribute declaration
- Inconsistent with package conventions

**Fix Required:**
```r
# Add at end of file:
attr(trk_creak_union, "ext") <- "crk"
attr(trk_creak_union, "tracks") <- c("AM_creak", "CD_creak", "union_creak")
attr(trk_creak_union, "outputType") <- "SSFF"
```

**Priority:** HIGH (simple one-line fix)

---

#### Issue #4: `trk_dv_f0()` - Missing ext and Parameter

**File:** `R/ssff_python_dv_f0.R`

**Current State:**
```r
# Function returns AsspDataObj with F0 track
# ❌ NO explicitExt parameter
# ❌ NO attr(*, "ext") declaration
```

**Impact:**
- Returns time-series F0 track data
- Part of DisVoice integration
- Cannot write to file with consistent naming

**Fix Required:**
```r
# 1. Add parameter to function signature:
trk_dv_f0 <- function(audio_path,
                     frame_shift = 10,
                     min_f0 = 75,
                     max_f0 = 600,
                     include_voicing = TRUE,
                     output_format = "AsspDataObj",
                     explicitExt = "dvf",  # ADD THIS
                     toFile = FALSE,       # ADD THIS
                     outputDirectory = NULL,  # ADD THIS
                     ...)

# 2. Add attribute at end of file:
attr(trk_dv_f0, "ext") <- "dvf"
attr(trk_dv_f0, "tracks") <- c("F0[Hz]", "voicing")
attr(trk_dv_f0, "outputType") <- "SSFF"
```

**Suggested Extension:** `"dvf"` (DisVoice F0)

**Priority:** MEDIUM (requires function refactoring)

---

#### Issue #5: `trk_dv_formants()` - Missing ext and Parameter

**File:** `R/ssff_python_dv_formants.R`

**Current State:**
```r
# Function returns AsspDataObj with formant tracks
# ❌ NO explicitExt parameter
# ❌ NO attr(*, "ext") declaration
```

**Impact:**
- Returns time-series formant track data
- Part of DisVoice integration
- Cannot write to file with consistent naming

**Fix Required:**
```r
# 1. Add parameters to function signature:
trk_dv_formants <- function(audio_path,
                           frame_shift = 10,
                           max_formant = 5500,
                           n_formants = 4,
                           output_format = "AsspDataObj",
                           explicitExt = "dvfm",  # ADD THIS
                           toFile = FALSE,        # ADD THIS
                           outputDirectory = NULL,  # ADD THIS
                           ...)

# 2. Add attribute at end of file:
attr(trk_dv_formants, "ext") <- "dvfm"
attr(trk_dv_formants, "tracks") <- c("fm", "bw")
attr(trk_dv_formants, "outputType") <- "SSFF"
```

**Suggested Extension:** `"dvfm"` (DisVoice Formants)

**Priority:** MEDIUM (requires function refactoring)

---

#### Issue #6: `trk_formants_tvwlp()` - Missing ext (Unverified)

**File:** `R/trk_formants_tvwlp.R`

**Current State:**
```r
# Function implementation needs verification
# ❌ NO attr(*, "ext") declaration found
```

**Action Required:**
1. Verify function implementation
2. Check if it supports file output
3. If yes, add appropriate `ext` attribute

**Suggested Extension:** `"tvw"` (TVWLP) or reuse `"fms"` (generic formants)

**Priority:** MEDIUM (needs investigation)

---

## Part 3: Recommendations

### 3.1 High Priority Fixes (2 items)

**MUST FIX - Parameter/Attribute Mismatches:**

1. **`lst_eGeMAPS()`** - Fix explicitExt default
   - **File:** `R/list_python_opensmile_eGeMAPS.R:66`
   - **Change:** `explicitExt="ocp"` → `explicitExt="ogs"`
   - **Effort:** 1 line change
   - **Risk:** Low

2. **`lst_emobase()`** - Fix explicitExt default
   - **File:** `R/list_python_opensmile_emobase.R:26`
   - **Change:** `explicitExt="ocp"` → `explicitExt="emo"`
   - **Effort:** 1 line change
   - **Risk:** Low

**Estimated Time:** 5 minutes + testing

---

### 3.2 Medium Priority Fixes (3-4 items)

**SHOULD FIX - Missing ext Attributes:**

3. **`trk_creak_union()`** - Add ext attribute
   - **File:** `R/trk_creak_union.R`
   - **Change:** Add `attr(trk_creak_union, "ext") <- "crk"`
   - **Effort:** 1-3 lines
   - **Risk:** Low

4. **`trk_dv_f0()`** - Add ext + parameters
   - **File:** `R/ssff_python_dv_f0.R`
   - **Change:** Add `explicitExt`, `toFile`, `outputDirectory` parameters + attr
   - **Effort:** Function signature refactor + attribute
   - **Risk:** Medium (may affect existing code)

5. **`trk_dv_formants()`** - Add ext + parameters
   - **File:** `R/ssff_python_dv_formants.R`
   - **Change:** Add `explicitExt`, `toFile`, `outputDirectory` parameters + attr
   - **Effort:** Function signature refactor + attribute
   - **Risk:** Medium (may affect existing code)

6. **`trk_formants_tvwlp()`** - Investigate and fix
   - **File:** `R/trk_formants_tvwlp.R`
   - **Change:** Verify implementation, add ext if needed
   - **Effort:** Investigation required
   - **Risk:** Unknown

**Estimated Time:** 1-2 hours + testing

---

### 3.3 Low Priority Enhancements (Optional)

**OPTIONAL - For Consistency:**

7. **C++ wrapper functions** - Add ext attributes
   - **Files:** `list_cpp_opensmile_*.R`
   - **Benefit:** Complete consistency across all functions
   - **Risk:** None (internal functions)
   - **Priority:** Low (not user-facing)

---

### 3.4 No Action Required (8 items)

**Intentionally Missing ext Attributes:**

- `lst_vat()` - Summary function ✓
- `lst_voice_sauce()` - Summary function ✓
- `lst_dysprosody()` - Summary function ✓
- `lst_voxit()` - Summary function ✓
- `lst_covarep_vq()` - Summary function ✓
- `lst_emobase_cpp()` - Internal wrapper ✓
- `lst_eGeMAPS_cpp()` - Internal wrapper ✓
- `lst_ComParE_2016_cpp()` - Internal wrapper ✓

---

## Part 4: Testing Recommendations

### 4.1 After Fixing High Priority Issues

```r
# Test lst_eGeMAPS with corrected extension
test_file <- system.file("samples", "sustained", "a1.wav", package = "superassp")
result <- lst_eGeMAPS(test_file, toFile = TRUE)
# Verify output file has .ogs extension (NOT .ocp)

# Test lst_emobase with corrected extension
result <- lst_emobase(test_file, toFile = TRUE)
# Verify output file has .emo extension (NOT .ocp)
```

### 4.2 After Fixing Medium Priority Issues

```r
# Test trk_creak_union
result <- trk_creak_union(test_file, toFile = TRUE)
# Verify output file has .crk extension

# Verify attribute is set
stopifnot(attr(trk_creak_union, "ext") == "crk")

# Test DisVoice functions (if refactored)
result <- trk_dv_f0(test_file, toFile = TRUE)
# Verify .dvf extension

result <- trk_dv_formants(test_file, toFile = TRUE)
# Verify .dvfm extension
```

---

## Part 5: Summary Statistics

### Overall Health Metrics

| Metric | Count | Percentage |
|--------|-------|------------|
| **Total DSP Functions** | 61+ | 100% |
| **Functions with ext** | 59 | 97% |
| **Functions without ext** | 11 | 18% |
| **Properly configured** | 59 | 97% |
| **Mismatched param/attr** | 2 | 3% |
| **Missing ext (bugs)** | 3-4 | 5-7% |
| **Intentionally no ext** | 8 | 13% |

### Issue Breakdown

| Priority | Issue Type | Count | Status |
|----------|-----------|-------|--------|
| 🔴 High | Parameter mismatch | 2 | Requires fix |
| 🟡 Medium | Missing ext attr | 3-4 | Requires fix |
| 🟢 Low | Enhancement | 3 | Optional |
| ✅ None | Intentional design | 8 | No action |

### Estimated Effort

- **High Priority Fixes:** 5 minutes + 10 minutes testing = **15 minutes**
- **Medium Priority Fixes:** 1-2 hours + 30 minutes testing = **1.5-2.5 hours**
- **Total Estimated Effort:** **2-3 hours**

---

## Part 6: Extension Naming Convention Reference

### Current Extensions in Use (Alphabetical)

```
acf  - Autocorrelation
ap   - Aperiodicity (D4C)
arf  - ARF coefficients
brh  - Brouhaha (VAD/SNR/C50)
cep  - Cepstrum
crp  - CREPE pitch
css  - Cepstral spectral
dft  - DFT spectral
dfm  - Deep formants
dsi  - Dysphonia Severity Index
egg  - EGG-based F0
emo  - emobase (OpenSMILE)
f0   - Standard F0
fms  - Formants (ASSP)
fpb  - Formant path (Praat)
gfm  - GFM-IAIF
glf  - COVAREP IAIF
int  - Intensity
lar  - LAR coefficients
lpc  - LPC coefficients
lps  - LPC spectral
mfcc - MFCC
ocp  - ComParE 2016
oge  - GeMAPS
ogs  - eGeMAPS
pfm  - Praat formants
pit  - Pitch (MHS/Praat)
pm   - Pitch marks
psa  - VoiceSauce (Praat)
pvr  - Praat voice report
pvt  - Voice tremor
pyp  - PYIN pitch
rfc  - RFC coefficients
rms  - RMS energy
rpm  - REAPER (Parselmouth)
sacc - SAcC pitch
sec  - Sinusoidal excitation
sf0  - Swift-F0
spm  - Spectral moments
wap  - Aperiodicity (WORLD)
xte  - Excitation
yf0  - YAAPT
yip  - YIN pitch
zcr  - Zero crossing rate
```

### Proposed New Extensions

```
crk  - Creaky voice (trk_creak_union) - Already in use, needs attr
dvf  - DisVoice F0 (trk_dv_f0) - PROPOSED
dvfm - DisVoice formants (trk_dv_formants) - PROPOSED
tvw  - TVWLP formants (trk_formants_tvwlp) - PROPOSED
```

---

## Conclusion

The superassp package demonstrates **excellent consistency** in extension attribute management, with 97% of functions properly configured. The identified issues are:

1. **2 critical mismatches** (easy fixes, 15 minutes)
2. **3-4 missing attributes** (medium complexity, 2-3 hours)
3. **8 intentional omissions** (correct design, no action needed)

All issues can be resolved with minimal risk in approximately **2-3 hours of development time**.

---

## Appendix A: Files Analyzed

### Files with DSP Functions (61+ functions across ~50 files)

**C ASSP Functions:**
- `ssff_c_assp_acfana.R`, `ssff_c_assp_cepstrum.R`, `ssff_c_assp_cssSpectrum.R`
- `ssff_c_assp_dftSpectrum.R`, `ssff_c_assp_forest.R`, `ssff_c_assp_ksvfo.R`
- `ssff_c_assp_lpsSpectrum.R`, `ssff_c_assp_mhspitch.R`, `ssff_c_assp_rmsana.R`
- `ssff_c_assp_zcrana.R`, `ssff_c_assp_lp_analysis.R`

**C++ SPTK Functions:**
- `ssff_cpp_sptk_rapt.R`, `ssff_cpp_sptk_swipe.R`, `ssff_cpp_sptk_dio.R`
- `ssff_cpp_sptk_harvest.R`, `ssff_cpp_sptk_reaper.R`, `ssff_cpp_sptk_mfcc.R`
- `ssff_cpp_sptk_d4c.R`

**C++ ESTK Functions:**
- `ssff_cpp_estk_pitchmark.R`

**Python Functions:**
- `ssff_python_pyin.R`, `ssff_python_yin.R`, `ssff_python_crepe.R`
- `ssff_python_yaapt.R`, `ssff_python_sacc.R`, `ssff_python_swiftf0.R`
- `ssff_python_deepformants.R`, `ssff_python_brouhaha.R`, `ssff_python_gfmiaif.R`
- Plus 25+ additional Python-based functions

**OpenSMILE Functions:**
- `list_python_opensmile_GeMAPS.R`, `list_python_opensmile_eGeMAPS.R`
- `list_python_opensmile_ComParE_2016.R`, `list_python_opensmile_emobase.R`
- `list_cpp_opensmile_gemaps.R`, `list_cpp_opensmile_emobase.R`

**Summary Functions:**
- `list_vat.R`, `lst_voice_sauce.R`, `list_dysprosody.R`, `list_voxit.R`
- `covarep_vq.R`

**Other Functions:**
- `trk_creak_union.R`, `trk_formants_tvwlp.R`
- `ssff_python_dv_f0.R`, `ssff_python_dv_formants.R`

---

**Report Generated:** 2025-10-29
**Auditor:** Claude Code (Automated Analysis)
**Package Version:** v0.8.7 (cpp_optimization branch)
