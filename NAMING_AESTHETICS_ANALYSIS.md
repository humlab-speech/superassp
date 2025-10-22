# Track Name Aesthetics for tibble/data.frame

## Problem Statement

The proposed bracket notation (e.g., `fo[Hz]`, `F1[Hz]`) creates syntactically awkward column names in R data frames and tibbles:

```r
# Awkward: requires backticks
df$`fo[Hz]`
df$`F1[Hz]`
df$`H1-H2c[dB]`

# Compare to:
df$fo_Hz
df$F1_Hz
df$H1_H2c_dB
```

## Alternative Notation Systems

### Option 1: Underscore Separator (Recommended)
**Format**: `measure_unit` (e.g., `fo_Hz`, `F1_Hz`, `H1_H2c_dB`)

**Pros**:
- ✅ No backticks required
- ✅ R-friendly syntax
- ✅ Works seamlessly with tidyverse
- ✅ Tab-completion friendly
- ✅ Clean in tibble output

**Cons**:
- ⚠️ Less visually distinct from measure name
- ⚠️ Deviates slightly from standard notation

**Examples**:
```r
# F0
fo_Hz

# Formants
F1_Hz, F2_Hz, F3_Hz, F4_Hz

# Bandwidths
B1_Hz, B2_Hz, B3_Hz

# Harmonics
H1_dB, H2_dB, A1_dB, A2_dB

# Harmonic differences
H1_H2_dB, H1_H2c_dB, H1_A1u_dB

# Voice quality
CPP_dB, HNR05_dB, SHR_dB
Jitter_local_pct, Shimmer_local_dB
```

**Tibble appearance**:
```r
# A tibble: 100 × 8
   frame_time fo_Hz F1_Hz F2_Hz F3_Hz H1_H2c_dB CPP_dB HNR05_dB
        <dbl> <dbl> <dbl> <dbl> <dbl>     <dbl>  <dbl>    <dbl>
 1      0.010   120   800  1200  2400       5.2   12.3     18.5
 2      0.020   125   820  1220  2380       5.5   12.1     18.2
```

### Option 2: Dot Separator
**Format**: `measure.unit` (e.g., `f0.Hz`, `F1.Hz`)

**Pros**:
- ✅ No backticks required
- ✅ Compact notation
- ✅ Common in other languages

**Cons**:
- ⚠️ Conflicts with S3 methods (e.g., `summary.f0` vs `summary` method for `fo` class)
- ⚠️ Can be confused with nested data structures
- ⚠️ Less common in R conventions

**Examples**:
```r
f0.Hz, F1.Hz, H1.H2c.dB
```

**Not recommended** due to S3 method conflicts.

### Option 3: Bracket + Attribute (Hybrid)
**Format**: Track name with brackets in track definition, but converted to underscore in data.frame

**Implementation**:
```r
# In AsspDataObj:
attr(obj, "tracks") <- c("fo[Hz]", "F1[Hz]")

# In as.data.frame.AsspDataObj:
# Convert brackets to underscores
names(df) <- gsub("\\[(.+)\\]", "_\\1", names(df))
# Result: fo_Hz, F1_Hz
```

**Pros**:
- ✅ Maintains bracket notation in SSFF files (standard)
- ✅ User-friendly column names in data frames
- ✅ Best of both worlds
- ✅ Automatic conversion

**Cons**:
- ⚠️ Naming differs between AsspDataObj and data.frame
- ⚠️ Slight complexity in conversion

**This is the RECOMMENDED approach**

### Option 4: Keep Brackets (Current Proposal)
**Format**: `measure[unit]` (e.g., `fo[Hz]`, `F1[Hz]`)

**Pros**:
- ✅ Matches scientific notation exactly
- ✅ Clear unit distinction
- ✅ Titze 2015 compliant
- ✅ Works with units package

**Cons**:
- ⚠️ Requires backticks in R
- ⚠️ Not tidyverse-friendly
- ⚠️ Awkward for users
- ⚠️ Tab-completion issues

**Only use if scientific exactness outweighs usability**

## Recommendation: Hybrid Approach

### Implementation Strategy

**1. Track Definition (SSFF/AsspDataObj)**
Use bracket notation for scientific compliance:
```r
attr(trk_rapt, "tracks") <- c("fo[Hz]")
attr(trk_forest, "tracks") <- c("F1[Hz]", "F2[Hz]", "B1[Hz]", "B2[Hz]")
```

**2. Automatic Conversion in as.data.frame/as_tibble**
Convert brackets to underscores for R-friendly names:
```r
as.data.frame.AsspDataObj <- function(x, ..., convert_units = TRUE,
                                      clean_names = TRUE) {
  # ... existing code ...

  if (clean_names) {
    # Convert brackets to underscores: fo[Hz] → fo_Hz
    names(df) <- gsub("\\[([^\\]]+)\\]", "_\\1", names(df))

    # Also clean up: H1-H2c[dB] → H1_H2c_dB
    names(df) <- gsub("-", "_", names(df))
  }

  df
}
```

**3. Alias System Supports Both**
```r
.track_name_aliases <- function() {
  list(
    # Old names
    "f0" = "fo[Hz]",
    "F0" = "fo[Hz]",

    # User-friendly names
    "fo_Hz" = "fo[Hz]",
    "F1_Hz" = "F1[Hz]",

    # Both work!
  )
}
```

## Updated Track Name Mapping

### F0 / Fundamental Frequency

| Old | SSFF/Attr | DataFrame | Notes |
|-----|-----------|-----------|-------|
| `fo` | `fo[Hz]` | `fo_Hz` | Lowercase, auto-convert |
| `F0` | `fo[Hz]` | `fo_Hz` | Normalize case |
| `fo` | `fo[Hz]` | `fo_Hz` | Add zero |
| `fo[Hz]` | `fo[Hz]` | `fo_Hz` | Add zero |

### Formants

| Old | SSFF/Attr | DataFrame | Notes |
|-----|-----------|-----------|-------|
| `F1` | `F1[Hz]` | `F1_Hz` | Add unit |
| `F2` | `F2[Hz]` | `F2_Hz` | Add unit |
| `F3` | `F3[Hz]` | `F3_Hz` | Add unit |
| `F4` | `F4[Hz]` | `F4_Hz` | Add unit |
| `fm` | `Fi[Hz]` | `Fi_Hz` | Generic formant |

### Bandwidths

| Old | SSFF/Attr | DataFrame | Notes |
|-----|-----------|-----------|-------|
| `B1` | `B1[Hz]` | `B1_Hz` | Add unit |
| `B2` | `B2[Hz]` | `B2_Hz` | Add unit |
| `B3` | `B3[Hz]` | `B3_Hz` | Add unit |
| `bw` | `Bi[Hz]` | `Bi_Hz` | Generic bandwidth |

### Harmonics

| Old | SSFF/Attr | DataFrame | Notes |
|-----|-----------|-----------|-------|
| `H1` | `H1[dB]` | `H1_dB` | Add unit |
| `H2` | `H2[dB]` | `H2_dB` | Add unit |
| `A1` | `A1[dB]` | `A1_dB` | Add unit |
| `H2K` | `H2k[dB]` | `H2k_dB` | Lowercase k |

### Harmonic Differences

| Old | SSFF/Attr | DataFrame | Notes |
|-----|-----------|-----------|-------|
| `H1H2` | `H1-H2[dB]` | `H1_H2_dB` | Hyphen + unit |
| `H1H2c` | `H1-H2c[dB]` | `H1_H2c_dB` | Corrected |
| `H1H2u` | `H1-H2u[dB]` | `H1_H2u_dB` | Uncorrected |
| `H1A1c` | `H1-A1c[dB]` | `H1_A1c_dB` | Corrected |

## Code Examples

### User Experience (Hybrid Approach)

```r
library(superassp)

# Generate F0 track
result <- trk_rapt("audio.wav", toFile = FALSE)

# AsspDataObj level: bracket notation
names(result)
# [1] "fo[Hz]"

attr(result, "tracks")
# [1] "fo[Hz]"

# Convert to data frame: automatic clean names
df <- as.data.frame(result)
names(df)
# [1] "frame_time" "fo_Hz"

# User-friendly access (no backticks!)
mean_f0 <- mean(df$fo_Hz, na.rm = TRUE)

# Tibble: clean output
library(tibble)
tbl <- as_tibble(result)
tbl
# # A tibble: 100 × 2
#    frame_time fo_Hz
#         <dbl> <dbl>
#  1      0.010   120
#  2      0.020   125
#  ...
```

### VoiceSauce Example

```r
# Generate VoiceSauce parameters
result <- trk_praat_sauce("audio.wav", toFile = FALSE)

# AsspDataObj: bracket notation (scientific)
attr(result, "tracks")
# [1] "fo[Hz]" "F1[Hz]" "F2[Hz]" "H1-H2c[dB]" ...

# Data frame: underscore notation (user-friendly)
df <- as.data.frame(result)
names(df)
# [1] "frame_time" "fo_Hz" "F1_Hz" "F2_Hz" "H1_H2c_dB" ...

# Easy access
df %>%
  mutate(
    spectral_tilt = H1_A1c_dB,
    breathiness = H1_H2c_dB
  ) %>%
  ggplot(aes(x = frame_time, y = fo_Hz)) +
  geom_line()
```

### Clean Names Parameter

```r
# Option 1: Clean names (default)
df <- as.data.frame(result, clean_names = TRUE)
names(df)
# [1] "fo_Hz" "F1_Hz"

# Option 2: Keep brackets (for scientific writing)
df <- as.data.frame(result, clean_names = FALSE)
names(df)
# [1] "fo[Hz]" "F1[Hz]"

# Access with backticks (less convenient)
df$`fo[Hz]`
```

## Implementation Changes

### Update as.data.frame.AsspDataObj

```r
#' @param clean_names Logical; if TRUE (default), convert bracket notation
#'   to underscore notation for R-friendly column names. If FALSE, keep
#'   bracket notation for scientific exactness.
as.data.frame.AsspDataObj <- function(x, ...,
                                      convert_units = TRUE,
                                      clean_names = TRUE) {
  # ... existing conversion code ...

  # Clean column names for R usability
  if (clean_names) {
    names(df) <- .clean_track_names(names(df))
  }

  df
}

#' Clean track names for R data frames
#' @keywords internal
.clean_track_names <- function(names) {
  # Convert brackets to underscores: [Hz] → _Hz
  names <- gsub("\\[([^\\]]+)\\]", "_\\1", names)

  # Convert hyphens to underscores: H1-H2 → H1_H2
  names <- gsub("-", "_", names)

  # Convert parentheses to underscores: (local) → _local
  names <- gsub("\\s*\\(([^)]+)\\)", "_\\1", names)

  # Remove spaces
  names <- gsub("\\s+", "_", names)

  names
}
```

### Update as_tibble.AsspDataObj

```r
#' @param clean_names Logical; if TRUE (default), convert bracket notation
#'   to underscore notation.
as_tibble.AsspDataObj <- function(x, ...,
                                  convert_units = TRUE,
                                  clean_names = TRUE) {
  df <- as.data.frame(x, convert_units = convert_units,
                      clean_names = clean_names, ...)
  tibble::as_tibble(df)
}
```

## Benefits of Hybrid Approach

1. **Scientific Compliance**: SSFF files use bracket notation (Titze 2015)
2. **User-Friendly**: Data frames use underscore notation (R conventions)
3. **Automatic**: Conversion happens transparently
4. **Flexible**: Users can choose via `clean_names` parameter
5. **Backwards Compatible**: Alias system supports all variants
6. **Publication-Ready**: Can export with scientific notation
7. **Analysis-Friendly**: Clean names for tidyverse workflows

## Updated Validation

All proposed names remain scientifically valid:
- ✅ SSFF level: `fo[Hz]`, `F1[Hz]`, `H1-H2c[dB]` (Titze/Nylén compliant)
- ✅ DataFrame level: `fo_Hz`, `F1_Hz`, `H1_H2c_dB` (R-friendly)
- ✅ Best of both worlds

## Conclusion

**Recommendation**: Adopt the **Hybrid Approach**

- Use bracket notation in track definitions (scientific)
- Auto-convert to underscores in data frames (practical)
- Provide `clean_names` parameter for user control
- Update all documentation to show both forms
- Implement `.clean_track_names()` helper function

This balances scientific rigor with R usability, making superassp both standards-compliant and user-friendly.
