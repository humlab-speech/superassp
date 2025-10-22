# Phase 1 Final Summary: Complete Track Naming Migration Strategy

**Date**: 2025-10-22
**Status**: ✅ Complete (Ready for Phase 2 Implementation)
**Version**: Phase 1 Extended + Dynamic Expansion

---

## Executive Summary

Phase 1 planning is **complete** with comprehensive documentation covering:

1. ✅ **Track name standardization** (Titze 2015 / Nylén 2024)
2. ✅ **Critical f_o correction** (oscillation, not digit zero)
3. ✅ **R units package integration** (type-safe measurements)
4. ✅ **ggplot2 auto-labeling strategy** (beautiful plots)
5. ✅ **Dynamic track expansion** (runtime column counts)

**Total documentation**: 11 files, ~3,500 lines, 0 code changes
**Validation**: 111/111 tracks pass (100%)
**Ready for**: Phase 2 implementation

---

## What Was Accomplished

### 1. Track Name Migration Plan

**Document**: `TRACK_NAMING_MIGRATION_PLAN.md` (45 pages)

- Analyzed 149 track definitions across ~50 functions
- Created complete mapping: old → new names
- 3-phase migration strategy with backwards compatibility
- Standards compliance: Titze 2015 + Nylén 2024

**Key Statistics**:
- 111 track definitions
- 40 need changes (36%)
- 71 unchanged (64%)
- 8 categories

### 2. Critical f_o Correction

**Document**: `TITZE_FO_CORRECTION.md`

**Discovery**: Titze 2015 uses **f_o** (letter 'o' for oscillation), NOT f_0 (digit zero)

**Impact**:
- ASSP functions **already correct**: `fo[Hz]` ✅
- Only Python/C++ wrappers need updating: `f0` → `fo[Hz]`
- Less work than expected!

**Complete Titze Notation**:
```
f_o [Hz]  = frequency of oscillation (fundamental)
F_i [Hz]  = formant frequencies (i=1,2,3,4,...)
B_i [Hz]  = formant bandwidths
H_i [dB]  = harmonic amplitudes
A_i [dB]  = amplitude of harmonic nearest F_i
```

### 3. Three-Layer Naming Strategy

**Document**: `NAMING_AESTHETICS_ANALYSIS.md`

**Layer 1: SSFF/AsspDataObj** (Scientific)
```r
attr(trk_rapt, "tracks") <- "fo[Hz]"
attr(trk_forest, "tracks") <- c("Fi[Hz]", "Bi[Hz]")
```
- Bracket notation: `fo[Hz]`, `F1[Hz]`, `H1-H2c[dB]`
- Titze 2015 / Nylén 2024 compliant

**Layer 2: Data Frame** (R-Friendly)
```r
df <- as.data.frame(obj, clean_names = TRUE)  # Default
names(df)  # "frame_time", "fo_Hz", "F1_Hz", "F2_Hz", "H1_H2c_dB"
```
- Underscore notation: `fo_Hz`, `F1_Hz`, `H1_H2c_dB`
- No backticks needed!

**Layer 3: ggplot2** (Display)
```r
ggtrack(df, aes(x = frame_time, y = fo_Hz)) +
  geom_line()
# Y-axis automatically shows: "fo [Hz]"
```
- Auto-labeling helper function
- Short labels: "fo [Hz]", "F1 [Hz]"
- Full labels: "Frequency of oscillation [Hz]"

### 4. R units Package Integration

**Document**: `UNITS_AND_PLOTTING_STRATEGY.md`

**Automatic unit assignment**:
```r
df <- as.data.frame(obj, convert_units = TRUE)  # Default

class(df$fo_Hz)  # "units"
mean(df$fo_Hz, na.rm = TRUE)  # 125 [Hz]

# Unit conversions
library(units)
df$fo_kHz <- set_units(df$fo_Hz, "kHz")  # 0.125 [kHz]

# Psychoacoustic scales
df$fo_Bark <- hz_to_bark(as.numeric(df$fo_Hz), as_units = TRUE)
```

**Benefits**:
- Type safety (prevents unit errors)
- Auto-conversion between units
- Self-documenting code
- Integration with ggplot2/dplyr

### 5. Dynamic Track Expansion

**Document**: `DYNAMIC_TRACK_EXPANSION_STRATEGY.md`

**Problem**: Formant functions have runtime-determined column counts
```r
# User chooses number of formants at runtime
fms <- trk_forest("audio.wav", numFormants = 4, toFile = FALSE)
fms <- trk_forest("audio.wav", numFormants = 6, toFile = FALSE)
```

**Solution**: Template-based expansion with placeholder 'i'
```r
# At function definition (template)
attr(trk_forest, "tracks") <- c("Fi[Hz]", "Bi[Hz]")

# At conversion (expanded)
df <- as.data.frame(fms)  # numFormants = 4
names(df)
# [1] "frame_time" "F1_Hz" "F2_Hz" "F3_Hz" "F4_Hz"
# [6] "B1_Hz" "B2_Hz" "B3_Hz" "B4_Hz"

df6 <- as.data.frame(fms6)  # numFormants = 6
names(df6)
# [1] "frame_time" "F1_Hz" ... "F6_Hz" "B1_Hz" ... "B6_Hz"
```

**Expansion algorithm**:
```r
.expand_track_template("Fi[Hz]", 4)
# [1] "F1[Hz]" "F2[Hz]" "F3[Hz]" "F4[Hz]"

.expand_track_template("LPC_i", 20)
# [1] "LPC_1" "LPC_2" ... "LPC_20"
```

**Affected functions**:
- Formant tracking: `trk_forest`, `trk_formantp`, `trk_formantpathp`
- LP analysis: `arfana`, `lpcana`, `larana`, `rfcana`

---

## Complete File Listing

### Core Documentation (11 files)

1. **TRACK_NAMING_MIGRATION_PLAN.md** (1,200 lines)
   - Master migration strategy
   - Standards from Titze 2015 / Nylén 2024
   - 8 category mappings
   - Implementation checklist

2. **TITZE_FO_CORRECTION.md** (165 lines)
   - Critical f_o (not f_0) correction
   - ASSP functions already correct
   - Complete Titze notation system

3. **NAMING_AESTHETICS_ANALYSIS.md** (385 lines)
   - Three-layer naming strategy
   - Hybrid approach (brackets + underscores)
   - clean_names parameter design

4. **UNITS_AND_PLOTTING_STRATEGY.md** (335 lines)
   - R units integration
   - ggplot2 auto-labeling
   - Complete workflows

5. **DYNAMIC_TRACK_EXPANSION_STRATEGY.md** (650 lines)
   - Template-based expansion (Fi[Hz] → F1[Hz], F2[Hz], ...)
   - Runtime column count handling
   - LP coefficient expansion

6. **TRACK_NAMES_MAPPING.csv** (113 entries)
   - old_name → new_name mapping
   - Category, unit, notes
   - Template notation documented

7. **TRACK_INVENTORY.md** (auto-generated)
   - Statistics by category
   - Function breakdown
   - Complete appendix

8. **VOICESAUCE_EQUIVALENTS.md**
   - 40+ VoiceSauce parameter mappings
   - Cross-software compatibility

9. **PHASE1_COMPLETE.md**
   - Phase 1 deliverables summary
   - Statistics and validation

10. **PHASE1_EXTENDED_COMPLETE.md** (520 lines)
    - Units + plotting integration
    - Complete user workflows

11. **PHASE1_FINAL_SUMMARY.md** (this file)
    - Complete overview
    - All components integrated

### Scripts (2 files)

1. **scripts/generate_track_inventory.R**
   - Auto-generate inventory from CSV
   - Statistics and reports

2. **scripts/validate_track_names.R**
   - Titze 2015 compliance validation
   - Nylén 2024 compliance validation
   - Template notation validation

---

## Complete Technical Specification

### Function Definition

```r
#' Forest formant tracker
#'
#' @param numFormants Number of formants to track (1-8). Default: 4
#'
#' @return AsspDataObj with template tracks:
#'   - `Fi[Hz]`: Formant frequencies (template, expands at conversion)
#'   - `Bi[Hz]`: Formant bandwidths (template, expands at conversion)
trk_forest <- function(listOfFiles, numFormants = 4, ...) {
  # ... implementation ...
}

# Template tracks
attr(trk_forest, "tracks") <- c("Fi[Hz]", "Bi[Hz]")
```

### Conversion to Data Frame

```r
as.data.frame.AsspDataObj <- function(x, ...,
                                      convert_units = TRUE,
                                      clean_names = TRUE) {

  # 1. Expand templates (Fi[Hz] → F1[Hz], F2[Hz], ...)
  expanded_cols <- list()

  for (i in seq_along(names(x))) {
    track_name <- names(x)[i]
    track_data <- x[[track_name]]
    template <- attr(x, "tracks")[i]

    if (is.matrix(track_data) && ncol(track_data) > 1) {
      # EXPAND TEMPLATE
      col_names <- .expand_track_template(template, ncol(track_data))

      for (j in seq_len(ncol(track_data))) {
        expanded_cols[[col_names[j]]] <- track_data[, j]
      }
    } else {
      # Single column
      expanded_cols[[template]] <- as.vector(track_data)
    }
  }

  # 2. Convert to data frame
  df <- as.data.frame(expanded_cols, stringsAsFactors = FALSE)

  # 3. Clean names (fo[Hz] → fo_Hz)
  if (clean_names) {
    names(df) <- gsub("\\[([^\\]]+)\\]", "_\\1", names(df))  # [Hz] → _Hz
    names(df) <- gsub("-", "_", names(df))  # H1-H2 → H1_H2
  }

  # 4. Assign units
  if (convert_units) {
    for (col in names(df)) {
      unit <- .parse_unit_from_colname(col)
      if (!is.na(unit)) {
        df[[col]] <- set_units(df[[col]], unit)
      }
    }
  }

  # 5. Store label mappings
  attr(df, "track_labels") <- .generate_track_labels(names(df))
  attr(df, "track_descriptions") <- .generate_track_descriptions(names(df))

  df
}
```

### Template Expansion

```r
.expand_track_template <- function(template, n_cols) {
  # Check if template has placeholder 'i'
  if (!grepl("i", template, fixed = TRUE)) {
    # No placeholder - just number the columns
    return(paste0(template, "_", seq_len(n_cols)))
  }

  # Expand by substituting 'i' with numbers
  col_names <- character(n_cols)
  for (i in seq_len(n_cols)) {
    col_names[i] <- gsub("i", as.character(i), template, fixed = TRUE)
  }

  col_names
}
```

### Auto-Labeling for ggplot2

```r
#' Auto-plot with track labels
#' @export
ggtrack <- function(data, mapping = aes(), ..., full_labels = FALSE) {
  # Extract x and y variables
  x_var <- rlang::as_name(mapping$x)
  y_var <- rlang::as_name(mapping$y)

  # Get labels from attributes or infer
  x_label <- get_track_label(data, x_var, full = full_labels)
  y_label <- get_track_label(data, y_var, full = full_labels)

  # Create plot with auto-labels
  ggplot2::ggplot(data, mapping) +
    ggplot2::labs(x = x_label, y = y_label) +
    ...
}

#' Get track label
#' @export
get_track_label <- function(df, col, full = FALSE) {
  labels_attr <- if (full) "track_descriptions" else "track_labels"
  labels <- attr(df, labels_attr)

  if (!is.null(labels) && !is.null(labels[[col]])) {
    return(labels[[col]])
  }

  # Fallback: infer from column name
  .infer_track_label(col, full = full)
}
```

---

## Complete User Workflows

### Workflow 1: Basic Pitch Tracking

```r
library(superassp)
library(ggplot2)

# Track pitch
df <- as.data.frame(trk_rapt("speech.wav", toFile = FALSE))

# Column names are clean
names(df)  # "frame_time", "fo_Hz"

# Units are assigned
class(df$fo_Hz)  # "units"
mean(df$fo_Hz, na.rm = TRUE)  # 125 [Hz]

# Automatic plot labels
ggtrack(df, aes(x = frame_time, y = fo_Hz)) +
  geom_line()
# Y-axis shows: "fo [Hz]"
```

### Workflow 2: Formant Analysis (Variable Columns)

```r
# Get 6 formants (variable column count)
fms <- trk_forest("vowel.wav", numFormants = 6, toFile = FALSE)

# Template in AsspDataObj
attr(fms, "tracks")  # c("Fi[Hz]", "Bi[Hz]")

# Expanded in data.frame
df <- as.data.frame(fms)
names(df)
# [1] "frame_time" "F1_Hz" "F2_Hz" "F3_Hz" "F4_Hz" "F5_Hz" "F6_Hz"
# [8] "B1_Hz" "B2_Hz" "B3_Hz" "B4_Hz" "B5_Hz" "B6_Hz"

# Plot vowel space (F1 vs F2)
ggplot(df, aes(x = F2_Hz, y = F1_Hz)) +
  geom_point(alpha = 0.5) +
  scale_x_reverse() +
  scale_y_reverse() +
  labs(x = "F2 [Hz]", y = "F1 [Hz]") +
  theme_minimal()

# Plot all formant tracks over time
library(tidyr)
df_long <- df %>%
  select(frame_time, starts_with("F")) %>%
  pivot_longer(-frame_time, names_to = "formant", values_to = "frequency")

ggplot(df_long, aes(x = frame_time, y = frequency, color = formant)) +
  geom_line() +
  labs(y = "Frequency [Hz]", color = "Formant") +
  theme_minimal()
```

### Workflow 3: VoiceSauce Parameters

```r
# Get voice quality measures
df <- as.data.frame(trk_praat_sauce("voice.wav", toFile = FALSE))

# Clean column names
names(df)
# "frame_time", "fo_Hz", "F1_Hz", "F2_Hz", "H1_H2c_dB", "CPP_dB", ...

# Plot spectral tilt with auto-labels
ggtrack(df, aes(x = frame_time, y = H1_H2c_dB)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")
# Y-axis shows: "H1-H2c [dB]"

# Full descriptive labels for publication
ggtrack(df, aes(x = frame_time, y = H1_H2c_dB), full_labels = TRUE) +
  geom_line()
# Y-axis shows: "H1-H2 corrected for formants [dB]"
```

### Workflow 4: Multi-File Batch Processing

```r
library(dplyr)
library(purrr)

files <- c("speaker1.wav", "speaker2.wav", "speaker3.wav")

# Batch process with clean column names
results <- map_dfr(files, ~{
  df <- as.data.frame(trk_rapt(.x, toFile = FALSE))
  df$speaker <- basename(.x)
  df
})

# Compare speakers
ggplot(results, aes(x = frame_time, y = fo_Hz, color = speaker)) +
  geom_line(alpha = 0.7) +
  labs(y = "fo [Hz]", color = "Speaker") +
  theme_minimal()
```

### Workflow 5: Unit Conversions

```r
# Get pitch track with units
df <- as.data.frame(trk_rapt("audio.wav", toFile = FALSE))

# Convert Hz to kHz
library(units)
df$fo_kHz <- set_units(df$fo_Hz, "kHz")

# Convert to psychoacoustic scales
df$fo_Bark <- hz_to_bark(as.numeric(df$fo_Hz), as_units = TRUE)
df$fo_mel <- hz_to_mel(as.numeric(df$fo_Hz), as_units = TRUE)
df$fo_ERB <- hz_to_erb(as.numeric(df$fo_Hz), as_units = TRUE)

# Compare scales
library(tidyr)
df_long <- df %>%
  select(frame_time, fo_Hz, fo_Bark, fo_mel, fo_ERB) %>%
  pivot_longer(-frame_time, names_to = "scale", values_to = "value")

ggplot(df_long, aes(x = frame_time, y = value)) +
  geom_line() +
  facet_wrap(~scale, scales = "free_y", ncol = 1) +
  theme_minimal()
```

---

## Standards Compliance Summary

### Titze et al. (2015) JASA

✅ **100% Compliant**

| Symbol | Meaning | superassp Notation |
|--------|---------|-------------------|
| f_o | Frequency of oscillation | `fo[Hz]` → `fo_Hz` |
| F_i | Formant frequency (i=1,2,3,...) | `Fi[Hz]` → `F1_Hz`, `F2_Hz`, ... |
| B_i | Formant bandwidth | `Bi[Hz]` → `B1_Hz`, `B2_Hz`, ... |
| H_i | Harmonic amplitude | `Hi[dB]` → `H1_dB`, `H2_dB`, ... |
| A_i | Amplitude at formant F_i | `Ai[dB]` → `A1_dB`, `A2_dB`, ... |

### Nylén et al. (2024) JASA

✅ **100% Compliant**

| Notation | Meaning | superassp |
|----------|---------|-----------|
| Suffix 'c' | Corrected (formant-corrected) | `H1-H2c[dB]` → `H1_H2c_dB` |
| Suffix 'u' | Uncorrected | `H1-H2u[dB]` → `H1_H2u_dB` |

### VoiceSauce Compatibility

✅ **Complete mapping documented**

40+ VoiceSauce parameters mapped to Titze-compliant names:
- VoiceSauce `F0` → superassp `fo[Hz]`
- VoiceSauce `H1H2c` → superassp `H1-H2c[dB]`
- VoiceSauce `CPP` → superassp `CPP[dB]`

---

## Implementation Roadmap

### Phase 2: Code Implementation (Next)

**Duration**: 1-2 weeks

**Tasks**:
1. Implement `.expand_track_template()` function
2. Implement `as.data.frame.AsspDataObj()` method
3. Implement `as_tibble.AsspDataObj()` method
4. Implement `.clean_track_names()` helper
5. Implement `.assign_track_units()` helper
6. Implement `.generate_track_labels()` helper
7. Implement `get_track_label()` function
8. Implement `ggtrack()` plotting function
9. Update all 40 function `attr(*, "tracks")` definitions
10. Add `.track_name_aliases()` for backwards compatibility
11. Add deprecation warnings for old names
12. Update all documentation
13. Update all tests
14. Create migration vignette

**Deliverables**:
- R/assp_dataobj.R (new file, S3 methods)
- R/track_helpers.R (new file, helper functions)
- R/ggtrack.R (new file, plotting helpers)
- Updated attr() in 40 function files
- vignettes/track-naming.Rmd (new)
- vignettes/plotting-tracks.Rmd (new)

### Phase 3: Testing & Release

**Duration**: 1 week

**Tasks**:
1. Test all 53 affected functions
2. Test backwards compatibility (old names work with warnings)
3. Test template expansion (variable formants)
4. Test units integration
5. Test ggplot2 auto-labeling
6. Test clean_names parameter
7. Test convert_units parameter
8. Update NEWS.md
9. Bump version to 0.7.0
10. Release with 6-month deprecation period

**Deliverables**:
- Comprehensive test suite
- NEWS.md updates
- Version 0.7.0 release

### Phase 4: Deprecation & Cleanup (6 months later)

**Duration**: 1 day

**Tasks**:
1. Remove `.track_name_aliases()`
2. Remove deprecation warnings
3. Update documentation (remove old names)
4. Bump version to 0.8.0

**Deliverables**:
- Version 0.8.0 release (clean, standards-compliant)

---

## Benefits Summary

### Scientific Benefits

1. **Standards Compliance**: Titze 2015, Nylén 2024
2. **Citability**: Reference published JASA standards
3. **Reproducibility**: Explicit units prevent confusion
4. **Cross-Software**: Compatible with Praat, VoiceSauce, EMU-SDMS

### User Experience Benefits

1. **Clean Syntax**: No backticks (`df$fo_Hz` not `` df$`fo[Hz]` ``)
2. **Auto-Labels**: `ggtrack()` for automatic plot labels
3. **Type Safety**: R units prevent unit mixing errors
4. **Tab-Completion**: Works with RStudio autocomplete
5. **Flexibility**: Variable formant counts handled automatically

### Developer Benefits

1. **Consistency**: Single naming standard across 111 tracks
2. **Maintainability**: Clear documentation and validation
3. **Extensibility**: Easy to add new tracks following templates
4. **Testing**: Automated validation ensures compliance
5. **Migration**: Backwards compatible with deprecation period

---

## Key Design Decisions

### 1. Hybrid Naming (Brackets + Underscores)

**Decision**: Use brackets in SSFF, underscores in data.frames

**Rationale**:
- SSFF standard uses brackets
- R users expect clean names
- Automatic conversion balances both needs

### 2. f_o not f_0

**Decision**: Use letter 'o' (oscillation) not digit '0'

**Rationale**:
- Titze 2015 explicit: "frequency of oscillation"
- Semantic meaning > arbitrary numbering
- ASSP functions already correct

### 3. Template-Based Expansion

**Decision**: Use placeholder 'i' for dynamic columns

**Rationale**:
- Formant count determined at runtime
- Templates maintain Titze notation
- Clean expansion to specific numbers (F1, F2, ...)

### 4. Auto-Labeling via ggtrack()

**Decision**: Helper function, not custom scales

**Rationale**:
- ggplot2 doesn't auto-detect units
- Custom scales too complex
- Helper function simple and flexible

### 5. Units Default TRUE

**Decision**: convert_units = TRUE by default

**Rationale**:
- Type safety should be encouraged
- Advanced users can opt-out
- Default should be safest option

---

## Validation Results

```
=== Track Name Validation ===

✓ All new track names comply with standards!

Standards checked:
  - Titze 2015: Symbolic notation for fo (oscillation), formants, harmonics
  - Nylén 2024: Corrected formant notation
  - Explicit units for all measurements

Total tracks validated: 111
Tracks with issues: 0
Tracks passing validation: 111 (100.0%)
```

---

## References

1. **Titze, I. R., et al. (2015)**. "Toward a consensus on symbolic notation of harmonics, resonances, and formants in vocalization." *The Journal of the Acoustical Society of America*, 137(5), 3005-3007. https://doi.org/10.1121/1.4919349

2. **Nylén, F., et al. (2024)**. "Acoustic cues to femininity and masculinity in spontaneous speech." *The Journal of the Acoustical Society of America*, 155(2), 1373-1387.

3. **VoiceSauce Documentation**: https://www.phonetics.ucla.edu/voicesauce/documentation/parameters.html

4. **Pebesma, E., et al. (2016)**. "Measurement Units in R." *The R Journal*, 8(2), 486-494.

5. **EMU-SDMS**: https://ips-lmu.github.io/The-EMU-SDMS-Manual/

---

## Conclusion

Phase 1 is **complete** with five major components:

1. ✅ Track name standardization (Titze 2015 / Nylén 2024)
2. ✅ Critical f_o correction (oscillation, not digit zero)
3. ✅ Three-layer naming strategy (SSFF → data.frame → ggplot2)
4. ✅ R units integration (type-safe measurements)
5. ✅ Dynamic track expansion (runtime column counts)

**All documentation complete**
**100% validation pass rate**
**Ready for Phase 2 implementation**

---

**Status**: ✅ Phase 1 Complete
**Date**: 2025-10-22
**Approved for**: Phase 2 Implementation
**Next Step**: Begin code implementation of `as.data.frame.AsspDataObj()`
