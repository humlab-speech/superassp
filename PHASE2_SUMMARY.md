# Phase 2 Implementation Summary

**Date:** 2025-10-22
**Branch:** units
**Status:** ✅ Complete

## Overview

Phase 2 implements the complete runtime infrastructure for the unified track naming system, enabling automatic template expansion, clean R-friendly column names, units assignment, and professional scientific notation with subscripts in plots.

## Commits

1. **9c731a6** - feat: Implement Phase 2 core - as.data.frame.AsspDataObj() with template expansion
2. **99f7758** - feat: Add plotmath expression support for scientific subscripts in axis labels
3. **3b48e86** - feat: Update function definitions to use template notation for multi-column tracks
4. **b6a59f4** - fix: Correct trk_forest fname and add tracks attribute setting

## New Files Created

### R/track_helpers.R (420 lines)
Internal helper functions for track name manipulation:

- `.has_placeholder(name)` - Detects template placeholders (`[A-Z]i` pattern)
- `.expand_track_template(template, n_cols)` - Expands `Fi[Hz]` → `F1[Hz]`, `F2[Hz]`, ...
- `.clean_track_names(names)` - Converts `fo[Hz]` → `fo_Hz`
- `.parse_unit_from_colname(col)` - Extracts unit suffix (`_Hz`, `_dB`, etc.)
- `.assign_track_units(df)` - Assigns R units package units to columns
- `.generate_track_labels(names)` - Creates short labels ("F1 [Hz]")
- `.generate_track_descriptions(names)` - Creates full labels ("First formant frequency [Hz]")
- `.get_track_label_mapping()` - Predefined label registry with 60+ common tracks

**Key Pattern:** Uniform placeholder detection using `[A-Z]i(\[|$)` regex.

### R/assp_dataobj.R (284 lines)
S3 methods for converting AsspDataObj to data.frame/tibble:

**Main Method:**
```r
as.data.frame.AsspDataObj(x, ...,
                          convert_units = TRUE,
                          clean_names = TRUE,
                          na.zeros = FALSE)
```

**Features:**
- Runtime template expansion based on actual matrix dimensions
- Automatic track name cleaning (brackets → underscores)
- Units assignment (Hz, dB, %, μs, s)
- Label attribute storage for plotting
- Handles `fo[Hz]`, `Fi[Hz]`, `LPCi`, etc.

**Also Includes:**
- `as_tibble.AsspDataObj()` - Tibble conversion wrapper
- `get_track_label(df, col, full = FALSE)` - Extract plotting labels

### R/track_labels_plotmath.R (230 lines)
Plotmath expression generation for scientific subscripts:

**Main Function:**
```r
get_track_label_expr(df, col, full = FALSE, use_subscripts = TRUE)
```

**Pattern Matching:**
- `fo` → `f[o]` (f subscript o, NOT f subscript zero)
- `F1` → `F[1]` (F subscript 1)
- `H1_H2c` → `H[1]-H[2c]` (harmonic differences with proper subscripts)
- `LPC12` → `LPC[12]` (LP coefficients)

**Implementation:**
- Uses R's built-in `bquote()` for expression creation
- No external dependencies (latex2exp not needed)
- Proper spacing with `~` operator: `f[o]~"[Hz]"`
- Pattern order matters: check specific patterns (H1-H2) before general ones (H1)

### R/ggtrack.R (157 lines - updated)
Auto-labeled ggplot2 plotting with subscripts:

```r
ggtrack(data, mapping = aes(), ...,
        full_labels = FALSE,
        use_subscripts = TRUE)
```

**Features:**
- Automatic axis label extraction from track metadata
- Plotmath expressions when `use_subscripts = TRUE`
- Plain text when `use_subscripts = FALSE`
- Seamless integration with ggplot2 layers

**Usage:**
```r
ggtrack(df, aes(x = frame_time, y = F1_Hz)) +
  geom_line()
# Y-axis shows: F₁ [Hz] (with subscript)
```

## Updated Function Definitions

Changed `attr(,"tracks")` to use template notation:

### Formant Analysis
```r
# trk_forest
attr(trk_forest, "tracks") <- c("Fi[Hz]", "Bi[Hz]")  # was: c("F[Hz]", "B[Hz]")

# trk_formantp (Praat/Parselmouth)
attr(trk_formantp, "tracks") <- c("fmi", "bwi", "lvi")  # was: c("fm", "bw", "lv")

# trk_formantpathp (Praat formant tracking)
attr(trk_formantpathp, "tracks") <- c("fmi", "bwi", "lvi")  # was: c("fm", "bw", "lv")
```

### Linear Prediction Analysis
```r
# LPC coefficients
attr(lpcana, "tracks") <- c("RMS[dB]", "gain[dB]", "LPCi")  # was: "LPC"

# Area function coefficients
attr(arfana, "tracks") <- c("RMS[dB]", "gain[dB]", "ARFi")  # was: "ARF"

# Log area ratio coefficients
attr(larana, "tracks") <- c("RMS[dB]", "gain[dB]", "LARi")  # was: "LAR"

# Reflection coefficients
attr(rfcana, "tracks") <- c("RMS[dB]", "gain[dB]", "RFCi")  # was: "RFC"
```

## Bug Fixes

### 1. trk_forest C Function Name
**Issue:** `fname = "trk_forest"` in processMediaFiles_LoadAndProcess call
**Fix:** Changed to `fname = "forest"` (actual ASSP C library function name)
**Impact:** trk_forest was completely broken since DSP function rename

### 2. Tracks Attribute Propagation
**Issue:** Template notation in function attributes wasn't being copied to returned AsspDataObj
**Fix:** Added code to set `attr(result, "tracks")` after simplification
**Implementation:**
```r
# After simplification for single file
if(!toFile) {
  tracks_attr <- attr(trk_forest, "tracks")
  if(n_files == 1) {
    attr(externalRes, "tracks") <- tracks_attr
  } else {
    for(i in seq_along(externalRes)) {
      attr(externalRes[[i]], "tracks") <- tracks_attr
    }
  }
}
```

### 3. Plotmath Pattern Matching Order
**Issue:** `H1_H2c_dB` matched single harmonic pattern instead of difference pattern
**Fix:** Check specific patterns (harmonic differences) before general patterns (single harmonics)
**Pattern Order:**
1. `H<n1>_H<n2>` (harmonic differences) - CHECK FIRST
2. `H<n>_A<n>` (harmonic to amplitude) - CHECK SECOND
3. `H<number>` (single harmonics) - CHECK THIRD

## Test Suite

Created `tests/testthat/test-track-naming-phase2.R` with 13 test groups:

1. ✅ Template expansion for formant tracks
2. ✅ Template expansion for LP coefficients
3. ✅ clean_names parameter functionality
4. ✅ Units assignment
5. ✅ Track label generation and storage
6. ✅ get_track_label() retrieval
7. ✅ Plotmath expression generation
8. ✅ Plotmath special cases (fo, frame_time)
9. ✅ na.zeros parameter
10. ✅ as_tibble.AsspDataObj()
11. ✅ ggtrack() with subscripts
12. ✅ Integration test: Full workflow
13. ✅ Error handling

## Three-Layer Naming Strategy

The implementation successfully enables the three-layer strategy:

### Layer 1: SSFF/AsspDataObj (Scientific Notation)
```r
attr(obj, "tracks") <- c("Fi[Hz]", "Bi[Hz]")
# Template notation with placeholders
```

### Layer 2: data.frame (R-Friendly)
```r
names(df)
# [1] "frame_time" "F1_Hz" "F2_Hz" "F3_Hz" "F4_Hz"
# [6] "B1_Hz" "B2_Hz" "B3_Hz" "B4_Hz"
# Clean, no backticks needed
```

### Layer 3: Plotting (Display Labels)
```r
# Short labels
get_track_label(df, "F1_Hz")
# [1] "F1 [Hz]"

# Full descriptive labels
get_track_label(df, "F1_Hz", full = TRUE)
# [1] "First formant frequency [Hz]"

# Plotmath expressions
get_track_label_expr(df, "F1_Hz", use_subscripts = TRUE)
# expression(F[1]~"[Hz]")  → renders as F₁ [Hz]
```

## Template Expansion Examples

### Runtime Expansion
```r
# User requests 4 formants
fms <- trk_forest("audio.wav", numFormants = 4, toFile = FALSE)
attr(fms, "tracks")  # c("Fi[Hz]", "Bi[Hz]")

# Automatic expansion during conversion
df <- as.data.frame(fms)
names(df)
# [1] "frame_time" "F1_Hz" "F2_Hz" "F3_Hz" "F4_Hz"
# [6] "B1_Hz" "B2_Hz" "B3_Hz" "B4_Hz"

# User requests 8 formants - same template, different expansion
fms8 <- trk_forest("audio.wav", numFormants = 8, toFile = FALSE)
df8 <- as.data.frame(fms8)
names(df8)
# [1] "frame_time" "F1_Hz" ... "F8_Hz" "B1_Hz" ... "B8_Hz"
```

### LP Coefficients
```r
lpc <- lpcana("audio.wav", order = 12, toFile = FALSE)
attr(lpc, "tracks")  # c("RMS[dB]", "gain[dB]", "LPCi")

df <- as.data.frame(lpc)
names(df)
# [1] "frame_time" "RMS_dB" "gain_dB" "LPC1" "LPC2" ... "LPC12"
```

## Known Issues

### devtools::load_all() Attribute Preservation
**Issue:** Function attributes set after function definition aren't preserved by `devtools::load_all()`
**Workaround:** Use installed package or source files directly
**Impact:** Development workflow only, doesn't affect users
**Status:** Known R/devtools limitation, not a package bug

## Integration with Existing Code

### Backward Compatibility
- All new parameters are optional with sensible defaults
- Existing code continues to work unchanged
- `clean_names = TRUE` by default (R-friendly names)
- `convert_units = TRUE` by default (type-safe units)
- `use_subscripts = TRUE` by default (professional plots)

### wrassp Compatibility
- Works with both wrassp and superassp AsspDataObj objects
- wrassp objects use fallback to `names(x)` when `tracks` attribute missing
- Template expansion works even without explicit templates

## Documentation

All new functions are fully documented with roxygen2:

- Parameter descriptions
- Return value specifications
- Usage examples
- Details sections
- Cross-references with `@seealso`

## Performance

- Template expansion: O(n) where n = number of columns
- Name cleaning: Regex-based, efficient
- Units assignment: Conditional, only when requested
- Label generation: Cached in attributes after first call

## Code Quality

- Consistent coding style
- Comprehensive error handling
- Input validation
- Clear function naming (`.internal_helper()` convention)
- Well-commented code
- Type safety with units package integration

## Next Steps

Phase 2 is complete. Recommended next steps:

1. **Phase 3: Documentation & Vignettes**
   - Create user guide for new naming system
   - Add plotting examples with subscripts
   - Document migration path for existing users

2. **Remaining Function Updates**
   - Update remaining ASSP C functions to set tracks attributes
   - Update Python-based functions (Praat, pysptk, etc.)
   - Update SPTK C++ functions

3. **Extended Features**
   - Add support for custom label mappings
   - Support for non-English labels
   - LaTeX output option for publications

4. **Testing & Validation**
   - Run full package test suite
   - Integration tests with emuR
   - Performance benchmarks

5. **Release Preparation**
   - Update NEWS.md
   - Version bump
   - Merge to master
   - CRAN submission preparation

---

**Implementation Time:** ~4 hours
**Files Changed:** 7 files (4 new, 3 updated)
**Lines of Code:** ~1200 new lines
**Test Coverage:** 13 test groups with 40+ assertions
