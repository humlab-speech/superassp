# Dynamic Track Expansion Strategy: Fi[Hz] → F1[Hz], F2[Hz], F3[Hz], ...

**Date**: 2025-10-22
**Status**: Phase 1.5 - Pre-Implementation
**Issue**: Multi-column tracks (formants, bandwidths) have unknown column count at function definition time

## Problem Statement

### The Issue

For formant tracking functions like `trk_forest()`, the number of formants is determined at **runtime** by the `numFormants` parameter, not at function definition time.

```r
# User can request any number of formants
fms <- trk_forest("audio.wav", numFormants = 3, toFile = FALSE)  # 3 formants
fms <- trk_forest("audio.wav", numFormants = 6, toFile = FALSE)  # 6 formants
```

### Current State

```r
attr(trk_forest, "tracks") <- c("F[Hz]", "B[Hz]")  # Generic placeholders
```

**Problem**: When converting to data.frame, we need actual column names:
- `F1[Hz]`, `F2[Hz]`, `F3[Hz]`, `F4[Hz]` (not `F[Hz]_1`, `F[Hz]_2`)
- `B1[Hz]`, `B2[Hz]`, `B3[Hz]`, `B4[Hz]` (not `B[Hz]_1`, `B[Hz]_2`)

### AsspDataObj Structure

From wrassp `forest()` output:
```r
fms <- forest("audio.wav", numFormants = 4, toFile = FALSE)

# Structure
class(fms)  # "AsspDataObj"
names(fms)  # c("fm", "bw")

dim(fms$fm)  # 289 rows × 4 columns (matrix)
colnames(fms$fm)  # NULL (no column names!)

dim(fms$bw)  # 289 rows × 4 columns (matrix)
colnames(fms$bw)  # NULL
```

**Key finding**: The matrices have NO column names - they're just numeric matrices.

## Affected Functions

### Functions with Dynamic Multi-Column Tracks

| Function | Track | Parameter | Range | Default |
|----------|-------|-----------|-------|---------|
| `trk_forest` | `F[Hz]`, `B[Hz]` | `numFormants` | 1-8 | 4 |
| `formant_burg` | `F[Hz]`, `B[Hz]` | `numFormants` | 1-8 | 4 |
| `formant_path_burg` | `F[Hz]`, `B[Hz]` | `numFormants` | 1-8 | 4 |
| `trk_formantp` | `Fi[Hz]`, `Bi[Hz]` | `numFormants` | 1-8 | 4 |
| `trk_formantpathp` | `Fi[Hz]`, `Bi[Hz]` | `numFormants` | 1-8 | 4 |
| `trk_snackf` | `Fi[Hz]`, `Bi[Hz]` | `numFormants` | 1-8 | 4 |
| `trk_spectral_momentsp` | `cog`, `sd`, `skewness`, `kurtosis` | (fixed) | 4 | 4 |
| `arfana`, `larana`, `lpcana`, `rfcana` | Various LP coefficients | `order` | 1-50 | 20 |

### Template Notation Needed

**Formant frequency**: `Fi[Hz]` where `i` = placeholder for 1, 2, 3, 4, ...
**Formant bandwidth**: `Bi[Hz]` where `i` = placeholder for 1, 2, 3, 4, ...
**LP coefficients**: `LPC_i` where `i` = 1 to `order`
**ARF coefficients**: `ARF_i` where `i` = 1 to `order`

## Proposed Solution: Template-Based Expansion

### Step 1: Define Track Templates

Use subscript notation with placeholders in `attr(*, "tracks")`:

```r
# In function definition
attr(trk_forest, "tracks") <- c("Fi[Hz]", "Bi[Hz]")  # 'i' is placeholder
attr(arfana, "tracks") <- c("ARF_i", "gain[dB]", "RMS[dB]")  # 'i' expands
attr(lpcana, "tracks") <- c("LPC_i", "gain[dB]", "RMS[dB]")  # 'i' expands
```

### Step 2: Runtime Expansion in as.data.frame()

When converting AsspDataObj → data.frame, expand templates based on actual matrix dimensions:

```r
as.data.frame.AsspDataObj <- function(x, ...,
                                      convert_units = TRUE,
                                      clean_names = TRUE) {

  # Get track metadata
  track_names <- attr(x, "tracks")
  if (is.null(track_names)) {
    track_names <- names(x)  # Fallback to list names
  }

  # Expand each track
  expanded_cols <- list()

  for (i in seq_along(names(x))) {
    track_name <- names(x)[i]
    track_data <- x[[track_name]]

    # Get template from track_names
    template <- track_names[i]
    if (is.null(template) || is.na(template)) {
      template <- track_name  # Fallback
    }

    # Check if matrix (multi-column)
    if (is.matrix(track_data) && ncol(track_data) > 1) {
      # EXPAND TEMPLATE
      col_names <- .expand_track_template(template, ncol(track_data))

      # Split matrix into individual columns
      for (j in seq_len(ncol(track_data))) {
        expanded_cols[[col_names[j]]] <- track_data[, j]
      }
    } else {
      # Single column - use template as-is (or remove 'i')
      col_name <- gsub("i", "1", template)  # Fi → F1 for single column
      expanded_cols[[col_name]] <- as.vector(track_data)
    }
  }

  # Convert to data frame
  df <- as.data.frame(expanded_cols, stringsAsFactors = FALSE)

  # Add time column
  sample_rate <- attr(x, "sampleRate")
  start_time <- attr(x, "startTime")
  n_frames <- nrow(df)

  df <- data.frame(
    frame_time = start_time + (seq_len(n_frames) - 1) / sample_rate,
    df,
    stringsAsFactors = FALSE
  )

  # Clean names if requested
  if (clean_names) {
    names(df) <- .clean_track_names(names(df))
  }

  # Assign units if requested
  if (convert_units) {
    df <- .assign_track_units(df)
  }

  # Store label mappings as attributes
  attr(df, "track_labels") <- .generate_track_labels(names(df))
  attr(df, "track_descriptions") <- .generate_track_descriptions(names(df))

  df
}
```

### Step 3: Template Expansion Function

```r
#' Expand track template to column names
#' @param template Track name template (e.g., "Fi[Hz]", "Bi[Hz]", "LPC_i")
#' @param n_cols Number of columns to generate
#' @return Character vector of expanded column names
#' @keywords internal
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

#' Examples:
#' .expand_track_template("Fi[Hz]", 4)
#' # [1] "F1[Hz]" "F2[Hz]" "F3[Hz]" "F4[Hz]"
#'
#' .expand_track_template("Bi[Hz]", 4)
#' # [1] "B1[Hz]" "B2[Hz]" "B3[Hz]" "B4[Hz]"
#'
#' .expand_track_template("LPC_i", 20)
#' # [1] "LPC_1"  "LPC_2"  "LPC_3"  ... "LPC_20"
#'
#' .expand_track_template("ARF_i", 12)
#' # [1] "ARF_1"  "ARF_2"  ... "ARF_12"
```

### Step 4: Update Function Definitions

Update all affected functions to use template notation:

```r
# trk_forest (R/ssff_c_assp_forest.R:206)
attr(trk_forest, "tracks") <- c("Fi[Hz]", "Bi[Hz]")  # Changed from F[Hz], B[Hz]

# trk_formantp
attr(trk_formantp, "tracks") <- c("Fi[Hz]", "Bi[Hz]", "lv")

# arfana
attr(arfana, "tracks") <- c("ARF_i", "gain[dB]", "RMS[dB]")

# lpcana
attr(lpcana, "tracks") <- c("LPC_i", "gain[dB]", "RMS[dB]")

# larana
attr(larana, "tracks") <- c("LAR_i", "gain[dB]", "RMS[dB]")

# rfcana
attr(rfcana, "tracks") <- c("RFC_i", "gain[dB]", "RMS[dB]")
```

## Complete Workflow Example

### Example 1: Formant Tracking with 4 Formants

```r
library(superassp)
library(ggplot2)

# Get formants (default: 4 formants)
fms <- trk_forest("vowel.wav", numFormants = 4, toFile = FALSE)

# Check attr before conversion
attr(fms, "tracks")
# [1] "Fi[Hz]" "Bi[Hz]"

# Convert to data frame
df <- as.data.frame(fms)

# Check column names
names(df)
# [1] "frame_time" "F1_Hz" "F2_Hz" "F3_Hz" "F4_Hz"
# [6] "B1_Hz" "B2_Hz" "B3_Hz" "B4_Hz"

# Plot F1 vs F2 (vowel space)
ggplot(df, aes(x = F2_Hz, y = F1_Hz)) +
  geom_point(alpha = 0.5) +
  scale_x_reverse() +
  scale_y_reverse() +
  labs(x = "F2 [Hz]", y = "F1 [Hz]", title = "Vowel space") +
  theme_minimal()
```

### Example 2: Variable Number of Formants

```r
# Request 6 formants
fms6 <- trk_forest("speech.wav", numFormants = 6, toFile = FALSE)

df6 <- as.data.frame(fms6)

names(df6)
# [1] "frame_time" "F1_Hz" "F2_Hz" "F3_Hz" "F4_Hz" "F5_Hz" "F6_Hz"
# [8] "B1_Hz" "B2_Hz" "B3_Hz" "B4_Hz" "B5_Hz" "B6_Hz"

# Plot formant tracks over time
library(tidyr)

df6_long <- df6 %>%
  select(frame_time, starts_with("F")) %>%
  pivot_longer(-frame_time, names_to = "formant", values_to = "frequency")

ggplot(df6_long, aes(x = frame_time, y = frequency, color = formant)) +
  geom_line() +
  labs(y = "Frequency [Hz]", color = "Formant") +
  theme_minimal()
```

### Example 3: LP Coefficients (Variable Order)

```r
# LPC analysis with order = 12
lpc <- lpcana("audio.wav", order = 12, toFile = FALSE)

# Template in attr
attr(lpc, "tracks")
# [1] "LPC_i" "gain[dB]" "RMS[dB]"

# Convert to data frame
df <- as.data.frame(lpc)

names(df)
# [1] "frame_time" "LPC_1" "LPC_2" ... "LPC_12" "gain_dB" "RMS_dB"

# Plot LPC coefficients heatmap
library(tidyr)

lpc_long <- df %>%
  select(frame_time, starts_with("LPC_")) %>%
  pivot_longer(-frame_time, names_to = "coefficient", values_to = "value")

ggplot(lpc_long, aes(x = frame_time, y = coefficient, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "LPC Coefficients over time") +
  theme_minimal()
```

## Benefits of Template-Based Approach

### 1. **Flexibility**
- Works with any number of formants/coefficients
- Determined at runtime, not function definition

### 2. **Consistency**
- Same template notation across all functions
- Predictable expansion rules

### 3. **Standards Compliance**
- Expands to Titze notation: `F1[Hz]`, `F2[Hz]`, etc.
- Maintains subscript semantics from templates

### 4. **User Experience**
- Column names are intuitive: `F1_Hz`, `F2_Hz`, etc.
- No confusing suffixes like `F_1`, `F_2`

### 5. **Backwards Compatibility**
- Old SSFF files still readable
- Template expansion happens only at conversion time

## Implementation Details

### Clean Names Conversion

After template expansion, apply clean names transformation:

```r
# After expansion
names(df) <- c("frame_time", "F1[Hz]", "F2[Hz]", "F3[Hz]", "F4[Hz]",
                               "B1[Hz]", "B2[Hz]", "B3[Hz]", "B4[Hz]")

# With clean_names = TRUE
names(df) <- .clean_track_names(names(df))
# [1] "frame_time" "F1_Hz" "F2_Hz" "F3_Hz" "F4_Hz"
# [6] "B1_Hz" "B2_Hz" "B3_Hz" "B4_Hz"
```

### Units Assignment

After clean names, assign units:

```r
# Detect Hz columns
hz_cols <- grep("_Hz$", names(df), value = TRUE)
for (col in hz_cols) {
  df[[col]] <- set_units(df[[col]], "Hz")
}

# Detect dB columns
db_cols <- grep("_dB$", names(df), value = TRUE)
for (col in db_cols) {
  df[[col]] <- set_units(df[[col]], "dB")
}
```

### Label Generation

Generate automatic labels for plotting:

```r
.generate_track_labels <- function(col_names) {
  labels <- list()

  for (col in col_names) {
    if (col == "frame_time") {
      labels[[col]] <- "Time [s]"
    } else if (grepl("^F[0-9]+_Hz$", col)) {
      # F1_Hz → "F1 [Hz]"
      labels[[col]] <- gsub("_", " [", col) %>% paste0("]")
    } else if (grepl("^B[0-9]+_Hz$", col)) {
      # B1_Hz → "B1 [Hz]"
      labels[[col]] <- gsub("_", " [", col) %>% paste0("]")
    } else if (grepl("_([A-Za-z]+)$", col)) {
      # Generic: param_unit → "param [unit]"
      labels[[col]] <- gsub("_([A-Za-z]+)$", " [\\1]", col)
    } else {
      labels[[col]] <- col
    }
  }

  labels
}

.generate_track_descriptions <- function(col_names) {
  descriptions <- list()

  for (col in col_names) {
    if (col == "frame_time") {
      descriptions[[col]] <- "Frame time [s]"
    } else if (grepl("^F([0-9]+)_Hz$", col)) {
      # F1_Hz → "First formant frequency [Hz]"
      n <- as.integer(gsub("^F([0-9]+)_Hz$", "\\1", col))
      ord <- c("First", "Second", "Third", "Fourth", "Fifth", "Sixth", "Seventh", "Eighth")
      descriptions[[col]] <- paste0(ord[n], " formant frequency [Hz]")
    } else if (grepl("^B([0-9]+)_Hz$", col)) {
      # B1_Hz → "First formant bandwidth [Hz]"
      n <- as.integer(gsub("^B([0-9]+)_Hz$", "\\1", col))
      ord <- c("First", "Second", "Third", "Fourth", "Fifth", "Sixth", "Seventh", "Eighth")
      descriptions[[col]] <- paste0(ord[n], " formant bandwidth [Hz]")
    } else {
      descriptions[[col]] <- labels[[col]]  # Fallback to short label
    }
  }

  descriptions
}
```

## Edge Cases and Handling

### Case 1: Single Formant Requested

```r
fms <- trk_forest("audio.wav", numFormants = 1, toFile = FALSE)
df <- as.data.frame(fms)

names(df)
# [1] "frame_time" "F1_Hz" "B1_Hz"  # Correct - still uses subscript 1
```

### Case 2: No Template (Fixed Columns)

```r
# Function with fixed number of columns
fo <- trk_ksvfo("audio.wav", toFile = FALSE)
attr(fo, "tracks")
# [1] "fo[Hz]"  # No 'i' placeholder - not a template

df <- as.data.frame(fo)
names(df)
# [1] "frame_time" "fo_Hz"  # Single column, no expansion
```

### Case 3: Mixed Templates and Fixed

```r
# arfana has both template (ARF_i) and fixed (gain, RMS)
attr(arfana, "tracks") <- c("ARF_i", "gain[dB]", "RMS[dB]")

arf <- arfana("audio.wav", order = 10, toFile = FALSE)
df <- as.data.frame(arf)

names(df)
# [1] "frame_time" "ARF_1" "ARF_2" ... "ARF_10" "gain_dB" "RMS_dB"
```

### Case 4: Template Without 'i' Placeholder

```r
# Fallback: just append numbers
attr(some_func, "tracks") <- c("coef[dB]")  # No 'i'

# If matrix has 5 columns
df <- as.data.frame(some_func(...))
names(df)
# [1] "frame_time" "coef[dB]_1" "coef[dB]_2" ... "coef[dB]_5"
```

## Update Track Names Mapping

Update CSV to document template notation:

```csv
file,line,function_name,old_name,category,new_name,unit,notes
ssff_c_assp_forest.R,206,trk_forest,F[Hz],formant,Fi[Hz],Hz,"Template: expands to F1[Hz], F2[Hz], ..."
ssff_c_assp_forest.R,206,trk_forest,B[Hz],bandwidth,Bi[Hz],Hz,"Template: expands to B1[Hz], B2[Hz], ..."
ssff_python_pm_pformantb.R,237,trk_formantp,fm,formant,Fi[Hz],Hz,"Template: expands to F1[Hz], F2[Hz], ..."
ssff_python_pm_pformantb.R,237,trk_formantp,bw,bandwidth,Bi[Hz],Hz,"Template: expands to B1[Hz], B2[Hz], ..."
ssff_c_assp_lp_analysis.R,340,arfana,ARF,other,ARF_i,,"Template: expands to ARF_1, ARF_2, ..."
ssff_c_assp_lp_analysis.R,514,lpcana,LPC,other,LPC_i,,"Template: expands to LPC_1, LPC_2, ..."
ssff_c_assp_lp_analysis.R,427,larana,LAR,other,LAR_i,,"Template: expands to LAR_1, LAR_2, ..."
ssff_c_assp_lp_analysis.R,253,rfcana,RFC,other,RFC_i,,"Template: expands to RFC_1, RFC_2, ..."
```

## Validation Updates

Update validation script to handle templates:

```r
validate_formant <- function(name) {
  errors <- c()

  # Check for template notation
  if (grepl("^Fi\\[Hz\\]$", name)) {
    # Template is valid
    return(errors)
  }

  # Or specific formant numbers
  if (grepl("^F[0-9]+\\[Hz\\]$", name)) {
    return(errors)
  }

  # Check for issues
  if (grepl("^f[0-9]", name)) {
    errors <- c(errors, "Should use uppercase 'F' for formants (Titze 2015)")
  }

  if (!grepl("\\[Hz\\]", name) && !grepl("i", name)) {
    errors <- c(errors, "Missing unit [Hz] or template placeholder 'i'")
  }

  errors
}
```

## Documentation Updates

### Function documentation example:

```r
#' Forest formant tracker
#'
#' @param numFormants Number of formants to track (1-8). Default: 4
#'
#' @return AsspDataObj with tracks:
#'   - `Fi[Hz]`: Formant frequencies (expands to F1[Hz], F2[Hz], ...)
#'   - `Bi[Hz]`: Formant bandwidths (expands to B1[Hz], B2[Hz], ...)
#'
#'   When converted to data.frame with `clean_names = TRUE`:
#'   - `F1_Hz`, `F2_Hz`, `F3_Hz`, `F4_Hz`, ...
#'   - `B1_Hz`, `B2_Hz`, `B3_Hz`, `B4_Hz`, ...
#'
#' @examples
#' # Get 4 formants (default)
#' fms <- trk_forest("vowel.wav", toFile = FALSE)
#' df <- as.data.frame(fms)
#' names(df)  # "frame_time", "F1_Hz", "F2_Hz", "F3_Hz", "F4_Hz",
#'            # "B1_Hz", "B2_Hz", "B3_Hz", "B4_Hz"
#'
#' # Get 6 formants
#' fms6 <- trk_forest("speech.wav", numFormants = 6, toFile = FALSE)
#' df6 <- as.data.frame(fms6)
#' names(df6)  # Includes F5_Hz, F6_Hz, B5_Hz, B6_Hz
```

## Implementation Checklist

- [ ] Implement `.expand_track_template()` function
- [ ] Update `as.data.frame.AsspDataObj()` with template expansion
- [ ] Update all formant functions with `Fi[Hz]`, `Bi[Hz]` notation
- [ ] Update all LP functions with template notation (`LPC_i`, `ARF_i`, etc.)
- [ ] Update `.generate_track_labels()` to handle numbered formants
- [ ] Update `.generate_track_descriptions()` for full labels
- [ ] Update TRACK_NAMES_MAPPING.csv with template notation
- [ ] Update validation script to accept templates
- [ ] Add tests for template expansion
- [ ] Add tests for variable formant counts
- [ ] Update documentation with template notation
- [ ] Create vignette: "Working with Multi-Column Tracks"

## Conclusion

Template-based expansion solves the dynamic column problem elegantly:

1. **At function definition**: Use template `Fi[Hz]`, `Bi[Hz]`, `LPC_i`
2. **In AsspDataObj**: Templates stored in `attr(*, "tracks")`
3. **At conversion**: Templates expanded based on actual matrix dimensions
4. **In data.frame**: Proper column names `F1_Hz`, `F2_Hz`, etc.
5. **In plots**: Automatic labels "F1 [Hz]", "First formant frequency [Hz]"

**Result**: Flexible, standards-compliant, user-friendly multi-column track handling.

---

**Status**: Ready for implementation in Phase 2
**Dependencies**: Requires `as.data.frame.AsspDataObj()` implementation
