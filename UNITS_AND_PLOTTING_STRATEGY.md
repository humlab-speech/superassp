# Units Package Integration & ggplot2 Aesthetics Strategy

**Date**: 2025-10-22
**Status**: Phase 1.5 - Before Code Implementation
**Goal**: Seamless integration of R `units` package with track names and ggplot2 plotting

## Problem Statement

We need a strategy that:

1. **Assigns units** to AsspDataObj track data using R `units` package
2. **Preserves units** through `as.data.frame()` and `as_tibble()` conversions
3. **Enables clean ggplot2 plotting** with automatic axis labels
4. **Balances three naming needs**:
   - Scientific notation in SSFF files: `fo[Hz]`
   - R-friendly column names: `fo_Hz`
   - ggplot2 aesthetic labels: "fo [Hz]" or "Frequency of oscillation [Hz]"

## Current State Analysis

### R units Package Behavior

```r
library(units)

# Units attach to numeric vectors
x <- set_units(c(100, 120, 150), 'Hz')
# Units: [Hz]
# [1] 100 120 150

# Units survive in data frames
df <- data.frame(fo = x)
class(df$fo)  # "units"

# But ggplot2 doesn't auto-label from units
library(ggplot2)
p <- ggplot(df, aes(x = 1:3, y = fo)) + geom_line()
p$labels$y  # "fo" (not "fo [Hz]")
```

**Key Finding**: Units package preserves unit metadata in vectors/data frames, but ggplot2 doesn't automatically use units for axis labels.

### Current superassp Implementation

From `demo_unit_conversion.R:27`:
```r
df_with_units <- as.data.frame(f0_obj, convert_units = TRUE)
class(df_with_units[["fo[Hz]"]])  # "units"
```

**Issue**: Column name `fo[Hz]` requires backticks in ggplot2:
```r
ggplot(df, aes(x = frame_time, y = `fo[Hz]`))  # Awkward!
```

## Proposed Solution: Three-Layer Naming Strategy

### Layer 1: SSFF/AsspDataObj Level (Scientific)

**Track names with brackets**:
```r
attr(trk_rapt, "tracks") <- "fo[Hz]"
attr(trk_forest, "tracks") <- c("F1[Hz]", "F2[Hz]", "B1[Hz]", "B2[Hz]")
```

**Why**:
- Titze 2015 compliant
- SSFF file format standard
- Unit extraction via regex: `\\[(.+)\\]`

### Layer 2: Data Frame Level (R-Friendly)

**Column names with underscores**:
```r
as.data.frame(obj, clean_names = TRUE)  # Default
# Columns: frame_time, fo_Hz, F1_Hz, F2_Hz

as.data.frame(obj, clean_names = FALSE)
# Columns: frame_time, fo[Hz], F1[Hz], F2[Hz]
```

**Implementation** (already proposed in NAMING_AESTHETICS_ANALYSIS.md:143):
```r
as.data.frame.AsspDataObj <- function(x, ...,
                                      convert_units = TRUE,
                                      clean_names = TRUE) {
  # ... existing code ...

  if (clean_names) {
    # fo[Hz] → fo_Hz
    names(df) <- gsub("\\[([^\\]]+)\\]", "_\\1", names(df))
    # H1-H2c[dB] → H1_H2c_dB
    names(df) <- gsub("-", "_", names(df))
  }

  # Assign units to columns
  if (convert_units) {
    for (col_name in names(df)) {
      unit <- .parse_unit_from_colname(col_name)  # Extract Hz, dB, etc.
      if (!is.na(unit)) {
        df[[col_name]] <- set_units(df[[col_name]], unit, mode = "standard")
      }
    }
  }

  df
}
```

### Layer 3: ggplot2 Aesthetics (Display)

**Problem**: ggplot2 doesn't auto-detect units for labels.

**Solution**: Custom labeling helper functions.

#### Option A: Manual Labels (Simple)

```r
library(superassp)
library(ggplot2)

df <- as.data.frame(trk_rapt("audio.wav", toFile = FALSE))

ggplot(df, aes(x = frame_time, y = fo_Hz)) +
  geom_line() +
  labs(y = "fo [Hz]")  # Manual label
```

**Pros**: Simple, explicit control
**Cons**: User must type labels manually

#### Option B: Smart Labels from Attributes (Recommended)

Store human-readable labels as data frame attributes:

```r
as.data.frame.AsspDataObj <- function(x, ..., clean_names = TRUE) {
  # ... conversion code ...

  # Store mapping: clean_name → display_label
  attr(df, "track_labels") <- list(
    fo_Hz = "fo [Hz]",
    F1_Hz = "F1 [Hz]",
    H1_H2c_dB = "H1-H2c [dB]"
  )

  # Store full descriptive labels (optional)
  attr(df, "track_descriptions") <- list(
    fo_Hz = "Frequency of oscillation [Hz]",
    F1_Hz = "First formant frequency [Hz]",
    H1_H2c_dB = "H1-H2 corrected [dB]"
  )

  df
}
```

Then provide helper function:

```r
#' Get plot label for a track
#' @param df Data frame from as.data.frame.AsspDataObj
#' @param col Column name
#' @param full Use full description (default: FALSE)
#' @export
get_track_label <- function(df, col, full = FALSE) {
  labels <- if (full) {
    attr(df, "track_descriptions")
  } else {
    attr(df, "track_labels")
  }

  if (!is.null(labels[[col]])) {
    return(labels[[col]])
  }

  # Fallback: convert underscore to bracket notation
  gsub("_([A-Za-z]+)$", " [\\1]", col)
}

#' Create ggplot with automatic track labels
#' @param df Data frame from as.data.frame.AsspDataObj
#' @param x X aesthetic (quoted or unquoted)
#' @param y Y aesthetic (quoted or unquoted)
#' @param ... Additional ggplot2 layers
#' @export
autoplot.AsspDataFrame <- function(df, x, y, ..., full_labels = FALSE) {
  x_var <- rlang::ensym(x)
  y_var <- rlang::ensym(y)

  x_label <- get_track_label(df, as.character(x_var), full = full_labels)
  y_label <- get_track_label(df, as.character(y_var), full = full_labels)

  ggplot2::ggplot(df, ggplot2::aes(x = !!x_var, y = !!y_var)) +
    ggplot2::labs(x = x_label, y = y_label) +
    ...
}
```

**Usage**:
```r
df <- as.data.frame(trk_rapt("audio.wav", toFile = FALSE))

# Automatic labels
autoplot(df, frame_time, fo_Hz) +
  geom_line()
# Y-axis shows: "fo [Hz]"

# Full descriptive labels
autoplot(df, frame_time, fo_Hz, full_labels = TRUE) +
  geom_line()
# Y-axis shows: "Frequency of oscillation [Hz]"

# Or manual with helper
ggplot(df, aes(x = frame_time, y = fo_Hz)) +
  geom_line() +
  labs(y = get_track_label(df, "fo_Hz"))
```

#### Option C: ggplot2 Scale Extensions

Use custom scale functions that read units from column metadata:

```r
#' Auto-label scales from units
#' @export
scale_y_track <- function(name = waiver(), ...) {
  # This would be called automatically if we register a ggplot2 method
  # for units class, but that's complex

  # Simpler: just provide a convenience function
  if (identical(name, waiver())) {
    # Try to get label from calling environment
    # This is tricky and might not be worth it
  }

  ggplot2::scale_y_continuous(name = name, ...)
}
```

**Verdict**: Too complex, Option B is better.

## Recommended Implementation Strategy

### Phase 1: Core Infrastructure (Immediate)

1. **Update `as.data.frame.AsspDataObj()`**:
   - Add `clean_names = TRUE` parameter (default)
   - Convert `fo[Hz]` → `fo_Hz` when `clean_names = TRUE`
   - Assign units using `units::set_units()` when `convert_units = TRUE`
   - Store `track_labels` and `track_descriptions` as attributes

2. **Update `as_tibble.AsspDataObj()`**:
   - Same behavior as `as.data.frame()`
   - Preserve attributes through tibble conversion

3. **Create label mapping system**:
   ```r
   .get_track_label_mapping <- function() {
     list(
       # Short labels (for plots)
       short = list(
         fo_Hz = "fo [Hz]",
         F1_Hz = "F1 [Hz]",
         F2_Hz = "F2 [Hz]",
         F3_Hz = "F3 [Hz]",
         F4_Hz = "F4 [Hz]",
         B1_Hz = "B1 [Hz]",
         B2_Hz = "B2 [Hz]",
         B3_Hz = "B3 [Hz]",
         B4_Hz = "B4 [Hz]",
         H1_dB = "H1 [dB]",
         H2_dB = "H2 [dB]",
         A1_dB = "A1 [dB]",
         A2_dB = "A2 [dB]",
         A3_dB = "A3 [dB]",
         H1_H2c_dB = "H1-H2c [dB]",
         H1_A1c_dB = "H1-A1c [dB]",
         H1_A2c_dB = "H1-A2c [dB]",
         H1_A3c_dB = "H1-A3c [dB]",
         CPP_dB = "CPP [dB]",
         HNR05_dB = "HNR05 [dB]",
         HNR15_dB = "HNR15 [dB]",
         HNR25_dB = "HNR25 [dB]",
         HNR35_dB = "HNR35 [dB]"
       ),

       # Full descriptive labels (for documentation/papers)
       full = list(
         fo_Hz = "Frequency of oscillation [Hz]",
         F1_Hz = "First formant frequency [Hz]",
         F2_Hz = "Second formant frequency [Hz]",
         F3_Hz = "Third formant frequency [Hz]",
         F4_Hz = "Fourth formant frequency [Hz]",
         B1_Hz = "First formant bandwidth [Hz]",
         B2_Hz = "Second formant bandwidth [Hz]",
         B3_Hz = "Third formant bandwidth [Hz]",
         B4_Hz = "Fourth formant bandwidth [Hz]",
         H1_dB = "First harmonic amplitude [dB]",
         H2_dB = "Second harmonic amplitude [dB]",
         A1_dB = "Amplitude of harmonic nearest F1 [dB]",
         A2_dB = "Amplitude of harmonic nearest F2 [dB]",
         A3_dB = "Amplitude of harmonic nearest F3 [dB]",
         H1_H2c_dB = "H1-H2 corrected for formants [dB]",
         H1_A1c_dB = "H1-A1 corrected for formants [dB]",
         H1_A2c_dB = "H1-A2 corrected for formants [dB]",
         H1_A3c_dB = "H1-A3 corrected for formants [dB]",
         CPP_dB = "Cepstral Peak Prominence [dB]",
         HNR05_dB = "Harmonics-to-Noise Ratio 0-500 Hz [dB]",
         HNR15_dB = "Harmonics-to-Noise Ratio 0-1500 Hz [dB]",
         HNR25_dB = "Harmonics-to-Noise Ratio 0-2500 Hz [dB]",
         HNR35_dB = "Harmonics-to-Noise Ratio 0-3500 Hz [dB]"
       )
     )
   }
   ```

### Phase 2: Helper Functions (Immediate)

```r
#' Get display label for track
#' @param df Data frame from AsspDataObj
#' @param col Column name
#' @param full Use full description
#' @export
get_track_label <- function(df, col, full = FALSE) {
  labels_attr <- if (full) "track_descriptions" else "track_labels"
  labels <- attr(df, labels_attr)

  if (!is.null(labels) && !is.null(labels[[col]])) {
    return(labels[[col]])
  }

  # Fallback: intelligent conversion
  .infer_track_label(col, full = full)
}

#' Infer track label from column name
#' @keywords internal
.infer_track_label <- function(col, full = FALSE) {
  # Try to get from mapping
  mapping <- .get_track_label_mapping()
  labels <- if (full) mapping$full else mapping$short

  if (!is.null(labels[[col]])) {
    return(labels[[col]])
  }

  # Fallback: convert fo_Hz → "fo [Hz]"
  if (grepl("_[A-Za-z]+$", col)) {
    return(gsub("_([A-Za-z]+)$", " [\\1]", col))
  }

  # No conversion
  col
}

#' Auto-plot with track labels
#' @param data Data frame from as.data.frame.AsspDataObj
#' @param mapping ggplot2 aesthetics
#' @param ... Additional layers
#' @param full_labels Use full descriptive labels
#' @export
ggtrack <- function(data, mapping = aes(), ..., full_labels = FALSE) {
  # Extract x and y variables from mapping
  x_var <- rlang::as_name(mapping$x)
  y_var <- rlang::as_name(mapping$y)

  # Get labels
  x_label <- get_track_label(data, x_var, full = full_labels)
  y_label <- get_track_label(data, y_var, full = full_labels)

  # Create plot
  ggplot2::ggplot(data, mapping) +
    ggplot2::labs(x = x_label, y = y_label) +
    ...
}
```

### Phase 3: User Workflows (Examples)

#### Workflow 1: Quick plotting with auto-labels

```r
library(superassp)
library(ggplot2)

# Load and process
df <- as.data.frame(trk_rapt("audio.wav", toFile = FALSE))

# Plot with automatic labels
ggtrack(df, aes(x = frame_time, y = fo_Hz)) +
  geom_line()
# Y-axis automatically shows: "fo [Hz]"
```

#### Workflow 2: Publication-quality with full labels

```r
df <- as.data.frame(trk_forest("audio.wav", toFile = FALSE))

ggtrack(df, aes(x = frame_time, y = F1_Hz), full_labels = TRUE) +
  geom_line() +
  theme_bw()
# Y-axis shows: "First formant frequency [Hz]"
```

#### Workflow 3: Multiple tracks with faceting

```r
library(tidyr)

df <- as.data.frame(trk_forest("audio.wav", toFile = FALSE))

# Reshape to long format
df_long <- df %>%
  select(frame_time, F1_Hz, F2_Hz, F3_Hz, F4_Hz) %>%
  pivot_longer(-frame_time, names_to = "formant", values_to = "frequency")

# Plot all formants
ggplot(df_long, aes(x = frame_time, y = frequency, color = formant)) +
  geom_line() +
  labs(y = "Frequency [Hz]", color = "Formant") +
  theme_minimal()
```

#### Workflow 4: Units-aware calculations

```r
df <- as.data.frame(trk_rapt("audio.wav", toFile = FALSE),
                    convert_units = TRUE)

# Units are preserved
class(df$fo_Hz)  # "units"

# Statistical analysis with units
mean(df$fo_Hz, na.rm = TRUE)
# 150 [Hz]

# Convert to other units
library(units)
df$fo_kHz <- set_units(df$fo_Hz, "kHz")
# Automatic conversion: 0.150 [kHz]

# Convert to psychoacoustic scales
df$fo_Bark <- hz_to_bark(as.numeric(df$fo_Hz), as_units = TRUE)
```

## Benefits of This Strategy

### 1. **Scientific Compliance**
- SSFF files use Titze 2015 notation: `fo[Hz]`
- Published standard, citable in papers

### 2. **R-Friendly Workflow**
- Data frames use clean names: `fo_Hz`
- No backticks required
- Tab-completion works

### 3. **Automatic Plotting**
- `ggtrack()` helper for auto-labels
- Fallback to intelligent inference
- Full descriptions available

### 4. **Units Integration**
- R `units` package for type safety
- Automatic unit checking
- Easy conversion between units

### 5. **Flexible Usage**
- Opt-out: `clean_names = FALSE`, `convert_units = FALSE`
- Manual labels still possible
- Progressive enhancement

## Implementation Checklist

- [ ] Update `as.data.frame.AsspDataObj()` with `clean_names` parameter
- [ ] Add unit assignment logic in conversion
- [ ] Create `.get_track_label_mapping()` with all track labels
- [ ] Implement `get_track_label()` helper
- [ ] Implement `ggtrack()` plotting function
- [ ] Add `track_labels` and `track_descriptions` attributes
- [ ] Update documentation with plotting examples
- [ ] Create vignette: "Plotting AsspDataObj with ggplot2"
- [ ] Add tests for label inference
- [ ] Add tests for units preservation

## Open Questions

1. **Should we create a custom S3 class for the data frame?**
   - E.g., `class(df) <- c("assp_df", "data.frame")`
   - Allows custom `plot()` and `print()` methods
   - Trade-off: More complexity vs better UX

2. **Should `ggtrack()` be in superassp or separate package?**
   - Could be in `superassp.plot` extension package
   - Or just include in main package (simpler)

3. **How to handle multi-track SSFF files in plotting?**
   - Forest has 8 tracks (F1-F4, B1-B4)
   - Automatic long-format conversion?
   - Helper: `reshape_tracks_long()`?

4. **Should we provide ggplot2 themes?**
   - E.g., `theme_phonetics()` with sensible defaults
   - Probably not in Phase 1, keep it simple

## Comparison to VoiceSauce/Praat

### VoiceSauce
- Exports CSV with headers like "F0(Hz)", "H1H2c(dB)"
- Users manually create ggplot labels
- **superassp improvement**: Automatic labels via `ggtrack()`

### Praat
- TextGrid annotation, formant tables
- Separate from plotting (use Praat's built-in plotting)
- **superassp improvement**: Direct R/ggplot2 integration

### Montreal Forced Aligner
- Outputs TextGrid, no direct signal processing results
- **superassp improvement**: End-to-end DSP + plotting

## Conclusion

The recommended strategy balances:
- **Scientific correctness**: Titze 2015 in SSFF files
- **R ergonomics**: Clean names in data frames
- **User experience**: Auto-labels in ggplot2
- **Type safety**: R units integration

**Next step**: Implement Phase 1 (core infrastructure) as part of track name migration Phase 2.
