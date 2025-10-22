# Phase 1 Extended: Track Name Migration + Units + Plotting Strategy

**Status**: ✅ Complete (Documentation)
**Date**: 2025-10-22
**Phase**: Pre-Implementation Planning

## Overview

This document summarizes the complete Phase 1 planning for track name migration, including the critical **f_o correction** and the **units/plotting integration strategy**.

## Critical Correction: f_o (not f_0)

**Discovery**: Titze et al. (2015) uses **f_o** (subscript letter 'o' for **oscillation**), NOT f_0 (subscript digit zero).

### Impact

**Good news**: ASSP functions were **already correct**!
```r
attr(trk_ksvfo, "tracks") <- "fo[Hz]"  # ✅ Already Titze-compliant
attr(fo, "tracks") <- "fo[Hz]"          # ✅ Already Titze-compliant
```

**Needs updating**: Only Python/C++ wrappers (20 functions)
```r
# Before
attr(trk_rapt, "tracks") <- "f0"        # ✗ Wrong
attr(trk_swiftf0, "tracks") <- "F0"     # ✗ Wrong

# After
attr(trk_rapt, "tracks") <- "fo[Hz]"    # ✓ Correct (Titze 2015)
attr(trk_swiftf0, "tracks") <- "fo[Hz]" # ✓ Correct (Titze 2015)
```

### Complete Titze 2015 Notation

| Symbol | Meaning | Unit | Example |
|--------|---------|------|---------|
| f_o | Frequency of oscillation (fundamental) | Hz | 100 Hz |
| F_i | Formant frequency (i=1,2,3,4...) | Hz | F₁ = 800 Hz |
| B_i | Formant bandwidth (i=1,2,3,4...) | Hz | B₁ = 80 Hz |
| H_i | Harmonic amplitude (i=1,2,3...) | dB | H₁ = 60 dB |
| A_i | Amplitude of harmonic nearest F_i | dB | A₁ = 55 dB |

**Reference**: Titze et al. (2015). JASA 137(5), 3005-3007.

## Three-Layer Naming Strategy

### Layer 1: SSFF/AsspDataObj (Scientific)

**Bracket notation** for scientific compliance:
```r
attr(trk_rapt, "tracks") <- "fo[Hz]"
attr(trk_forest, "tracks") <- c("F1[Hz]", "F2[Hz]", "F3[Hz]", "F4[Hz]",
                                "B1[Hz]", "B2[Hz]", "B3[Hz]", "B4[Hz]")
attr(trk_praat_sauce, "tracks") <- c("fo[Hz]", "F1[Hz]", "H1-H2c[dB]")
```

**Why**:
- Titze 2015 / Nylén 2024 compliant
- Standard SSFF format
- Citable in publications

### Layer 2: Data Frame (R-Friendly)

**Underscore notation** for R ergonomics:
```r
df <- as.data.frame(obj, clean_names = TRUE)  # Default
names(df)
# [1] "frame_time" "fo_Hz" "F1_Hz" "F2_Hz" "H1_H2c_dB"

# No backticks needed!
mean(df$fo_Hz, na.rm = TRUE)
```

**Automatic conversion**:
```r
as.data.frame.AsspDataObj <- function(x, ..., clean_names = TRUE) {
  # ...
  if (clean_names) {
    # fo[Hz] → fo_Hz
    names(df) <- gsub("\\[([^\\]]+)\\]", "_\\1", names(df))
    # H1-H2c[dB] → H1_H2c_dB
    names(df) <- gsub("-", "_", names(df))
  }
  # ...
}
```

### Layer 3: ggplot2 Aesthetics (Display)

**Auto-labels** for beautiful plots:
```r
library(superassp)
library(ggplot2)

df <- as.data.frame(trk_rapt("audio.wav", toFile = FALSE))

# Automatic track labels
ggtrack(df, aes(x = frame_time, y = fo_Hz)) +
  geom_line()
# Y-axis shows: "fo [Hz]"

# Full descriptive labels for publications
ggtrack(df, aes(x = frame_time, y = fo_Hz), full_labels = TRUE) +
  geom_line()
# Y-axis shows: "Frequency of oscillation [Hz]"
```

**Implementation**: Helper function with label mapping
```r
get_track_label(df, "fo_Hz")            # "fo [Hz]"
get_track_label(df, "fo_Hz", full = TRUE)  # "Frequency of oscillation [Hz]"
```

## R units Package Integration

### Automatic Unit Assignment

```r
df <- as.data.frame(obj, convert_units = TRUE)  # Default

class(df$fo_Hz)     # "units"
print(df$fo_Hz[1])  # 120 [Hz]

# Units-aware calculations
mean(df$fo_Hz, na.rm = TRUE)  # 125 [Hz]
range(df$fo_Hz, na.rm = TRUE) # 100 [Hz] 150 [Hz]
```

### Unit Conversions

```r
library(units)

# Convert Hz to kHz
df$fo_kHz <- set_units(df$fo_Hz, "kHz")
# Automatic: 0.125 [kHz]

# Convert to psychoacoustic scales
df$fo_Bark <- hz_to_bark(as.numeric(df$fo_Hz), as_units = TRUE)
df$fo_mel <- hz_to_mel(as.numeric(df$fo_Hz), as_units = TRUE)
```

### Benefits

1. **Type safety**: Prevents mixing incompatible units
2. **Auto-conversion**: Units package handles conversions
3. **Documentation**: Units are self-documenting
4. **Integration**: Works with ggplot2, dplyr, etc.

## Complete User Workflows

### Workflow 1: Basic Pitch Tracking

```r
library(superassp)
library(ggplot2)

# Track pitch
df <- as.data.frame(trk_rapt("speech.wav", toFile = FALSE))

# Quick plot with auto-labels
ggtrack(df, aes(x = frame_time, y = fo_Hz)) +
  geom_line() +
  theme_minimal()
# Y-axis automatically shows: "fo [Hz]"

# Statistical summary with units
summary(df$fo_Hz)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  80 [Hz] 110 [Hz] 120 [Hz] 125 [Hz] 140 [Hz] 180 [Hz]
```

### Workflow 2: Formant Analysis

```r
# Track formants
df <- as.data.frame(trk_forest("vowel.wav", toFile = FALSE))

# Reshape for faceted plot
library(tidyr)
df_long <- df %>%
  select(frame_time, F1_Hz, F2_Hz, F3_Hz) %>%
  pivot_longer(-frame_time, names_to = "formant", values_to = "frequency")

# Plot all formants
ggplot(df_long, aes(x = frame_time, y = frequency, color = formant)) +
  geom_line() +
  labs(y = "Frequency [Hz]", color = "Formant") +
  theme_bw()
```

### Workflow 3: VoiceSauce Parameters

```r
# Get voice quality measures
df <- as.data.frame(trk_praat_sauce("voice.wav", toFile = FALSE))

# Plot spectral tilt
ggtrack(df, aes(x = frame_time, y = H1_H2c_dB)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_minimal()
# Y-axis shows: "H1-H2c [dB]"
```

### Workflow 4: Multi-File Batch Processing

```r
library(dplyr)
library(purrr)

files <- c("speaker1.wav", "speaker2.wav", "speaker3.wav")

# Batch process
results <- map_dfr(files, ~{
  df <- as.data.frame(trk_rapt(.x, toFile = FALSE))
  df$file <- .x
  df
})

# Compare speakers
ggplot(results, aes(x = frame_time, y = fo_Hz, color = file)) +
  geom_line(alpha = 0.7) +
  labs(y = "fo [Hz]", color = "Speaker") +
  theme_minimal()
```

### Workflow 5: Publication-Quality Figure

```r
library(superassp)
library(ggplot2)
library(patchwork)

# Track multiple parameters
fo_df <- as.data.frame(trk_rapt("speech.wav", toFile = FALSE))
formant_df <- as.data.frame(trk_forest("speech.wav", toFile = FALSE))

# Pitch contour
p1 <- ggtrack(fo_df, aes(x = frame_time, y = fo_Hz), full_labels = TRUE) +
  geom_line(color = "steelblue") +
  theme_bw() +
  theme(axis.title.x = element_blank())

# Formant tracks
p2 <- ggplot(formant_df, aes(x = frame_time)) +
  geom_line(aes(y = F1_Hz, color = "F1")) +
  geom_line(aes(y = F2_Hz, color = "F2")) +
  geom_line(aes(y = F3_Hz, color = "F3")) +
  labs(y = "Formant frequency [Hz]", x = "Time [s]", color = "Formant") +
  theme_bw()

# Combine
p1 / p2 +
  plot_annotation(title = "Acoustic analysis of /a/ vowel")
```

## Phase 1 Deliverables Summary

### Documentation Files (9 files)

1. **TRACK_NAMING_MIGRATION_PLAN.md** (45 pages)
   - Complete 3-phase strategy
   - 149 track definitions analyzed
   - 8 category mappings

2. **TRACK_NAMES_MAPPING.csv** (111 entries)
   - old_name → new_name mapping
   - Category, unit, notes

3. **TRACK_INVENTORY.md**
   - Statistics by category
   - Function breakdown
   - Complete appendix

4. **VOICESAUCE_EQUIVALENTS.md**
   - 40+ parameter mappings
   - VoiceSauce → superassp equivalents

5. **NAMING_AESTHETICS_ANALYSIS.md**
   - Hybrid approach rationale
   - R usability considerations
   - clean_names implementation

6. **TITZE_FO_CORRECTION.md**
   - Critical f_o (not f_0) correction
   - ASSP functions already correct
   - Simplified migration scope

7. **UNITS_AND_PLOTTING_STRATEGY.md** (NEW)
   - Three-layer naming strategy
   - R units integration
   - ggplot2 auto-labeling
   - Complete workflows

8. **PHASE1_COMPLETE.md**
   - Summary of Phase 1
   - Statistics and validation

9. **PHASE1_EXTENDED_COMPLETE.md** (this file)
   - Complete overview
   - Units + plotting integration

### Scripts (2 files)

1. **scripts/generate_track_inventory.R**
   - Automated inventory generation
   - CSV → Markdown conversion

2. **scripts/validate_track_names.R**
   - Titze 2015 compliance checking
   - Nylén 2024 compliance checking
   - 100% validation pass rate

## Validation Results

```
✓ All new track names comply with standards!

Standards checked:
  - Titze 2015: Symbolic notation for fo (oscillation), formants, harmonics
  - Nylén 2024: Corrected formant notation
  - Explicit units for all measurements

Total tracks validated: 111
Tracks with issues: 0
Tracks passing validation: 111 (100.0%)
```

## Migration Statistics

| Category | Tracks | Need Changes | Unchanged |
|----------|--------|--------------|-----------|
| fo | 24 | 20 (83%) | 4 (17%) |
| Formant | 2 | 2 (100%) | 0 (0%) |
| Bandwidth | 2 | 2 (100%) | 0 (0%) |
| Harmonic | 9 | 9 (100%) | 0 (0%) |
| Harmonic diff | 8 | 8 (100%) | 0 (0%) |
| Voice quality | 11 | 11 (100%) | 0 (0%) |
| Spectral | 10 | 0 (0%) | 10 (100%) |
| Metadata | 16 | 0 (0%) | 16 (100%) |
| Other | 45 | 4 (9%) | 41 (91%) |
| **TOTAL** | **111** | **40 (36%)** | **71 (64%)** |

## Implementation Roadmap

### Phase 2: Code Implementation (Next)

**Goal**: Implement track name changes with backwards compatibility

**Tasks**:
1. Update all `attr(*, "tracks")` definitions (40 functions)
2. Implement `.track_name_aliases()` for old → new mapping
3. Add `clean_names` parameter to `as.data.frame.AsspDataObj()`
4. Add `convert_units` parameter for R units integration
5. Implement `get_track_label()` helper
6. Implement `ggtrack()` plotting function
7. Add deprecation warnings for old names
8. Update all function documentation
9. Update all tests
10. Create migration vignette

**Duration**: 1-2 weeks

### Phase 3: Testing & Release (Following)

**Goal**: Comprehensive testing and v0.7.0 release

**Tasks**:
1. Test all 53 affected functions
2. Test backwards compatibility (old names still work)
3. Test units integration
4. Test ggplot2 auto-labeling
5. Create plotting vignette
6. Update NEWS.md
7. Bump version to 0.7.0
8. Release with deprecation period (6 months)

**Duration**: 1 week

### Phase 4: Deprecation & Cleanup (Future)

**Goal**: Remove old names in v0.8.0 (after 6 months)

**Tasks**:
1. Remove aliases for old names
2. Remove deprecation warnings
3. Update documentation to remove old names
4. Bump version to 0.8.0

## Benefits Summary

### Scientific Benefits

1. **Standards Compliance**: Titze 2015, Nylén 2024
2. **Citability**: Reference published JASA standards
3. **Cross-Software**: Compatible with Praat, VoiceSauce
4. **Reproducibility**: Explicit units prevent confusion

### User Experience Benefits

1. **No Backticks**: Clean column names (`fo_Hz` not `` `fo[Hz]` ``)
2. **Auto-Labels**: `ggtrack()` for automatic plot labels
3. **Type Safety**: R units prevent unit errors
4. **Tab-Completion**: Works with RStudio autocomplete

### Developer Benefits

1. **Consistency**: Single naming standard across 111 tracks
2. **Maintainability**: Clear mapping documentation
3. **Extensibility**: Easy to add new tracks following standard
4. **Testing**: Validation script ensures compliance

## Key Design Decisions

### Decision 1: Hybrid Naming (Brackets + Underscores)

**Rationale**:
- SSFF files need scientific notation (`fo[Hz]`)
- R data frames need clean names (`fo_Hz`)
- Best of both worlds via `clean_names` parameter

**Alternative considered**: Only underscores everywhere
**Rejected because**: SSFF standard uses brackets

### Decision 2: f_o not f_0

**Rationale**:
- Titze 2015 uses letter 'o' for oscillation
- Semantic meaning > arbitrary zero
- ASSP functions already correct

**Alternative considered**: Change ASSP to f_0
**Rejected because**: They were already right!

### Decision 3: Auto-labeling via ggtrack()

**Rationale**:
- ggplot2 doesn't auto-detect units for labels
- Manual labels tedious and error-prone
- Helper function best balance

**Alternative considered**: Custom ggplot2 scales
**Rejected because**: Too complex, limited benefit

### Decision 4: Opt-in Units (Default TRUE)

**Rationale**:
- Units provide safety and documentation
- Advanced users can disable if needed
- Default should be safest option

**Alternative considered**: Opt-out (default FALSE)
**Rejected because**: Units should be encouraged

## References

1. **Titze, I. R., et al. (2015)**. "Toward a consensus on symbolic notation of harmonics, resonances, and formants in vocalization." *The Journal of the Acoustical Society of America*, 137(5), 3005-3007. https://doi.org/10.1121/1.4919349

2. **Nylén, F., et al. (2024)**. "Acoustic cues to femininity and masculinity in spontaneous speech." *The Journal of the Acoustical Society of America*, 155(2), 1373-1387.

3. **VoiceSauce Documentation**: https://www.phonetics.ucla.edu/voicesauce/documentation/parameters.html

4. **R units package**: Pebesma, E., Mailund, T., & Hiebert, J. (2016). Measurement Units in R. *The R Journal*, 8(2), 486-494.

## Conclusion

Phase 1 Extended is **complete** with three major components:

1. ✅ **Track name migration plan** (Titze 2015 / Nylén 2024)
2. ✅ **Critical f_o correction** (letter o, not digit zero)
3. ✅ **Units & plotting strategy** (seamless ggplot2 integration)

All documentation complete, 100% validation pass rate, ready for Phase 2 implementation.

**Next step**: Begin Phase 2 code implementation with backwards-compatible track name updates.

---

**Status**: ✅ Phase 1 Extended Complete
**Date**: 2025-10-22
**Approved for**: Phase 2 Implementation
