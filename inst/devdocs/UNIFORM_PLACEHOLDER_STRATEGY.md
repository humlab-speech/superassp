# Uniform Placeholder Strategy: Consistent 'i' Notation

**Date**: 2025-10-22
**Status**: Phase 1 Refinement
**Issue**: Inconsistency between `Fi[Hz]`, `Bi[Hz]` vs `LPC_i`, `ARF_i`

## Problem Identified

Current strategy has inconsistent placeholder notation:

```r
# Formants: no underscore
attr(trk_forest, "tracks") <- c("Fi[Hz]", "Bi[Hz]")

# LP coefficients: with underscore
attr(lpcana, "tracks") <- c("LPC_i", "gain[dB]", "RMS[dB]")
```

**Issue**: Why the difference? Should be uniform!

## Proposed Uniform Strategy

**Rule**: The last `i` before either `[<unit>]` or end-of-string is the placeholder for sequence expansion, **without underscore**.

### Pattern Detection

```r
# Regex pattern: uppercase letter + 'i' + (bracket or end-of-string)
pattern <- "[A-Z]i(\\[|$)"
```

This correctly identifies:
- ✅ `Fi[Hz]` → `F1[Hz], F2[Hz], F3[Hz], ...`
- ✅ `Bi[Hz]` → `B1[Hz], B2[Hz], B3[Hz], ...`
- ✅ `LPCi` → `LPC1, LPC2, LPC3, ...`
- ✅ `ARFi` → `ARF1, ARF2, ARF3, ...`
- ✅ `LARi` → `LAR1, LAR2, LAR3, ...`
- ✅ `RFCi` → `RFC1, RFC2, RFC3, ...`
- ✅ `Hi[dB]` → `H1[dB], H2[dB], H3[dB], ...`
- ✅ `Ai[dB]` → `A1[dB], A2[dB], A3[dB], ...`

And correctly **excludes**:
- ✅ `foi[Hz]` → no expansion (lowercase before 'i')
- ✅ `intensity` → no expansion (lowercase, no bracket/end after 'i')
- ✅ `pitch[Hz]` → no expansion (no 'i' before bracket)
- ✅ `gain[dB]` → no expansion (no 'i' before bracket)

## Updated Notation

### Before (Inconsistent)

```r
# Formants
attr(trk_forest, "tracks") <- c("Fi[Hz]", "Bi[Hz]")

# LP coefficients
attr(lpcana, "tracks") <- c("LPC_i", "gain[dB]", "RMS[dB]")
attr(arfana, "tracks") <- c("ARF_i", "gain[dB]", "RMS[dB]")
attr(larana, "tracks") <- c("LAR_i", "gain[dB]", "RMS[dB]")
attr(rfcana, "tracks") <- c("RFC_i", "gain[dB]", "RMS[dB]")
```

### After (Uniform)

```r
# Formants (no change)
attr(trk_forest, "tracks") <- c("Fi[Hz]", "Bi[Hz]")

# LP coefficients (remove underscore)
attr(lpcana, "tracks") <- c("LPCi", "gain[dB]", "RMS[dB]")
attr(arfana, "tracks") <- c("ARFi", "gain[dB]", "RMS[dB]")
attr(larana, "tracks") <- c("LARi", "gain[dB]", "RMS[dB]")
attr(rfcana, "tracks") <- c("RFCi", "gain[dB]", "RMS[dB]")
```

## Expansion Examples

### Formants (4 formants)

```r
template <- "Fi[Hz]"
expand_track_template(template, 4)
# [1] "F1[Hz]" "F2[Hz]" "F3[Hz]" "F4[Hz]"

# With clean_names = TRUE:
# [1] "F1_Hz" "F2_Hz" "F3_Hz" "F4_Hz"
```

### LP Coefficients (order = 12)

```r
template <- "LPCi"
expand_track_template(template, 12)
# [1] "LPC1" "LPC2" "LPC3" "LPC4" "LPC5" "LPC6"
# [7] "LPC7" "LPC8" "LPC9" "LPC10" "LPC11" "LPC12"

# With clean_names = TRUE (no change, no unit):
# [1] "LPC1" "LPC2" "LPC3" ... "LPC12"
```

### Harmonics

```r
# Note: Harmonics can also use template notation if needed
template <- "Hi[dB]"
expand_track_template(template, 5)
# [1] "H1[dB]" "H2[dB]" "H3[dB]" "H4[dB]" "H5[dB]"

# With clean_names = TRUE:
# [1] "H1_dB" "H2_dB" "H3_dB" "H4_dB" "H5_dB"
```

## Updated Implementation

### Detection Function

```r
#' Detect if track name has placeholder 'i'
#' @param name Track name template
#' @return Logical: TRUE if has placeholder
#' @keywords internal
.has_placeholder <- function(name) {
  # Pattern: uppercase letter + 'i' + (bracket or end-of-string)
  grepl("[A-Z]i(\\[|$)", name)
}
```

### Expansion Function

```r
#' Expand track template to column names
#' @param template Track name template (e.g., "Fi[Hz]", "LPCi", "ARFi")
#' @param n_cols Number of columns to generate
#' @return Character vector of expanded column names
#' @keywords internal
#' @examples
#' .expand_track_template("Fi[Hz]", 4)
#' # [1] "F1[Hz]" "F2[Hz]" "F3[Hz]" "F4[Hz]"
#'
#' .expand_track_template("LPCi", 12)
#' # [1] "LPC1" "LPC2" ... "LPC12"
.expand_track_template <- function(template, n_cols) {

  # Check if template has placeholder
  if (!.has_placeholder(template)) {
    # No placeholder - just number the columns with underscore
    # (This is fallback for non-standard templates)
    return(paste0(template, "_", seq_len(n_cols)))
  }

  # Expand by substituting last 'i' before [ or end with numbers
  col_names <- character(n_cols)

  for (j in seq_len(n_cols)) {
    # Replace: (uppercase)i(bracket or end) → (uppercase)(number)(bracket or end)
    col_names[j] <- sub("([A-Z])i(\\[|$)", paste0("\\1", j, "\\2"), template)
  }

  col_names
}
```

### Test Cases

```r
# Test all cases
.expand_track_template("Fi[Hz]", 3)
# [1] "F1[Hz]" "F2[Hz]" "F3[Hz]"

.expand_track_template("Bi[Hz]", 3)
# [1] "B1[Hz]" "B2[Hz]" "B3[Hz]"

.expand_track_template("LPCi", 5)
# [1] "LPC1" "LPC2" "LPC3" "LPC4" "LPC5"

.expand_track_template("ARFi", 10)
# [1] "ARF1" "ARF2" ... "ARF10"

.expand_track_template("Hi[dB]", 4)
# [1] "H1[dB]" "H2[dB]" "H3[dB]" "H4[dB]"

.expand_track_template("Ai[dB]", 3)
# [1] "A1[dB]" "A2[dB]" "A3[dB]"

# Non-templates (no expansion)
.has_placeholder("fo[Hz]")     # FALSE
.has_placeholder("gain[dB]")   # FALSE
.has_placeholder("intensity")  # FALSE
```

## Benefits of Uniform Strategy

### 1. **Consistency**
- Same rule for all templates: `[A-Z]i` pattern
- No need to remember when to use underscore

### 2. **Simplicity**
- Single regex pattern detects all placeholders
- Clear rule: "uppercase + i before bracket/end"

### 3. **Titze Compliance**
- `Fi` matches F_i subscript notation
- `Bi` matches B_i subscript notation
- `Hi` matches H_i subscript notation
- `Ai` matches A_i subscript notation

### 4. **Unambiguous**
- Lowercase 'i' in middle of word ignored: `intensity`, `pitch`
- Only uppercase + i pattern triggers expansion
- No false positives

### 5. **Aesthetics**
- `LPCi` cleaner than `LPC_i`
- `ARFi` cleaner than `ARF_i`
- Consistent with formant notation

## Edge Cases

### Case 1: Multiple 'i' in Name

```r
# Hypothetical: what if track name is "CPPi[dB]"?
.expand_track_template("CPPi[dB]", 3)
# [1] "CPP1[dB]" "CPP2[dB]" "CPP3[dB]"

# Only the LAST uppercase + i is replaced
# This is correct behavior for templates
```

### Case 2: Lowercase Before 'i'

```r
# "foi[Hz]" - should NOT expand (lowercase 'o' before 'i')
.has_placeholder("foi[Hz]")  # FALSE

# This protects against accidental expansion
```

### Case 3: 'i' in Middle of Word

```r
# "pitch[Hz]" - no 'i' before bracket
.has_placeholder("pitch[Hz]")  # FALSE

# "intensity" - 'i' not at end
.has_placeholder("intensity")  # FALSE

# Both correctly excluded
```

### Case 4: No Unit Bracket

```r
# "LPCi" - ends with i, no unit
.has_placeholder("LPCi")  # TRUE
.expand_track_template("LPCi", 3)
# [1] "LPC1" "LPC2" "LPC3"

# Pattern works with or without unit bracket
```

## Updated Track Mapping

Update `TRACK_NAMES_MAPPING.csv`:

```csv
file,line,function_name,old_name,category,new_name,unit,notes
ssff_c_assp_forest.R,206,trk_forest,F[Hz],formant,Fi[Hz],Hz,"Template: expands to F1[Hz], F2[Hz], ..."
ssff_c_assp_forest.R,206,trk_forest,B[Hz],bandwidth,Bi[Hz],Hz,"Template: expands to B1[Hz], B2[Hz], ..."
ssff_c_assp_lp_analysis.R,340,arfana,ARF,other,ARFi,,"Template: expands to ARF1, ARF2, ..."
ssff_c_assp_lp_analysis.R,514,lpcana,LPC,other,LPCi,,"Template: expands to LPC1, LPC2, ..."
ssff_c_assp_lp_analysis.R,427,larana,LAR,other,LARi,,"Template: expands to LAR1, LAR2, ..."
ssff_c_assp_lp_analysis.R,253,rfcana,RFC,other,RFCi,,"Template: expands to RFC1, RFC2, ..."
```

**Note**: Changed from `LPC_i`, `ARF_i` to `LPCi`, `ARFi` (removed underscore)

## Validation Update

Update validation script to use uniform pattern:

```r
validate_template <- function(name) {
  errors <- c()

  # Check for old-style underscore notation
  if (grepl("_i(\\[|$)", name)) {
    errors <- c(errors, "Use uniform notation: 'LPCi' not 'LPC_i'")
  }

  # Check if follows pattern
  if (.has_placeholder(name)) {
    # Valid template - ensure uppercase before i
    if (!grepl("[A-Z]i", name)) {
      errors <- c(errors, "Placeholder 'i' must be preceded by uppercase letter")
    }
  }

  errors
}
```

## Documentation Examples

### Function Documentation Template

```r
#' Linear Predictive Coding Analysis
#'
#' @param order LP order (number of coefficients). Default: 20
#'
#' @return AsspDataObj with template tracks:
#'   - `LPCi`: LP coefficients (expands to LPC1, LPC2, ..., LPC<order>)
#'   - `gain[dB]`: Prediction gain
#'   - `RMS[dB]`: RMS amplitude
#'
#'   When converted to data.frame:
#'   - Column names: `LPC1`, `LPC2`, ..., `LPC20` (for order=20)
#'   - Fixed columns: `gain_dB`, `RMS_dB`
#'
#' @examples
#' # Get LPC coefficients (order = 12)
#' lpc <- lpcana("audio.wav", order = 12, toFile = FALSE)
#' df <- as.data.frame(lpc)
#' names(df)
#' # [1] "frame_time" "LPC1" "LPC2" ... "LPC12" "gain_dB" "RMS_dB"
lpcana <- function(listOfFiles, order = 20, ...) {
  # ...
}

attr(lpcana, "tracks") <- c("LPCi", "gain[dB]", "RMS[dB]")
```

## Complete Affected Functions

### Formant Tracking (Fi, Bi)

| Function | Old Tracks | New Tracks |
|----------|-----------|------------|
| `trk_forest` | `F[Hz]`, `B[Hz]` | `Fi[Hz]`, `Bi[Hz]` |
| `trk_formant` | `fm`, `bw` | `Fi[Hz]`, `Bi[Hz]` |
| `trk_snackf` | (varies) | `Fi[Hz]`, `Bi[Hz]` |

### LP Analysis (Coefficient Templates)

| Function | Old Tracks | New Tracks |
|----------|-----------|------------|
| `lpcana` | `LPC` | `LPCi`, `gain[dB]`, `RMS[dB]` |
| `arfana` | `ARF` | `ARFi`, `gain[dB]`, `RMS[dB]` |
| `larana` | `LAR` | `LARi`, `gain[dB]`, `RMS[dB]` |
| `rfcana` | `RFC` | `RFCi`, `gain[dB]`, `RMS[dB]` |

### Potential Future Uses

If we ever need dynamic harmonic expansion:

```r
# Dynamic number of harmonics
attr(some_func, "tracks") <- c("Hi[dB]", "Ai[dB]")

# Expands to:
# H1[dB], H2[dB], H3[dB], ..., Hn[dB]
# A1[dB], A2[dB], A3[dB], ..., An[dB]
```

## Migration Impact

### Files to Update

1. **DYNAMIC_TRACK_EXPANSION_STRATEGY.md**
   - Update all `_i` to just `i`
   - Update expansion examples

2. **TRACK_NAMES_MAPPING.csv**
   - Change `LPC_i` → `LPCi`
   - Change `ARF_i` → `ARFi`
   - Change `LAR_i` → `LARi`
   - Change `RFC_i` → `RFCi`

3. **Implementation code**
   - Update `.expand_track_template()` regex
   - Add `.has_placeholder()` helper

### Breaking Changes

**None** - this is still in Phase 1 (documentation only)

## Conclusion

The uniform placeholder strategy is:

✅ **Consistent**: Same rule for all templates
✅ **Simple**: Single regex pattern
✅ **Titze-compliant**: Matches subscript notation
✅ **Unambiguous**: No false positives
✅ **Aesthetic**: Cleaner without underscores

**Recommendation**: Adopt uniform `[A-Z]i` pattern for all placeholders.

**Rule**: The last `i` after an uppercase letter and before `[` or end-of-string is the sequence placeholder.

---

**Status**: ✅ Refinement Complete
**Next**: Update documentation to use uniform notation
