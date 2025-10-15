# Track Name Differences: SuperASP vs wrassp

## Overview

SuperASP and wrassp packages produce functionally equivalent outputs but use different capitalization conventions for track names.

## Differences by Function

| Function | SuperASP Track Names | wrassp Track Names |
|----------|---------------------|-------------------|
| acfana   | ACF                 | acf               |
| rfcana   | RMS[dB], gain[dB], RFC | rms, gain, rfc    |
| rmsana   | RMS[dB]             | rms               |
| zcrana   | ZCR[Hz]             | zcr               |

## Patterns

### Case Difference
- **SuperASP**: Uses UPPERCASE track names (e.g., "ACF", "RMS", "RFC", "ZCR")
- **wrassp**: Uses lowercase track names (e.g., "acf", "rms", "rfc", "zcr")

### Bracket Notation with Units
- **SuperASP**: Includes units in square brackets (e.g., "RMS[dB]", "ZCR[Hz]", "gain[dB]")
- **wrassp**: No bracket notation, just the base name (e.g., "rms", "zcr", "gain")

## Impact

### No Impact On:
- ✓ Data values - identical or nearly identical
- ✓ Dimensions - same number of frames and columns
- ✓ Data structure - both return AsspDataObj
- ✓ Functionality - both work correctly

### Minor Impact On:
- Track name matching - requires normalized comparison (case + bracket notation)
- Code that hardcodes track names - needs to handle both cases

## Solution

When comparing or accessing tracks from both packages, use a normalization function:

```r
# Helper function to normalize track names
# Removes bracket notation like [dB], [Hz], etc. and converts to lowercase
normalize_track_name <- function(name) {
  tolower(gsub("\\[.*?\\]", "", name))
}

# Normalized track name comparison
expect_equal(
  normalize_track_name(names(result_superassp)),
  normalize_track_name(names(result_wrassp))
)

# Normalized track access
track_name_normalized <- normalize_track_name("RMS[dB]")  # becomes "rms"

# Find matching track
superassp_tracks <- names(result_superassp)
track <- superassp_tracks[normalize_track_name(superassp_tracks) == track_name_normalized][1]

# Access data
data <- result_superassp[[track]]
```

## Recommendation

When writing code that should work with both packages:

```r
# Helper function to normalize track names
normalize_track_name <- function(name) {
  tolower(gsub("\\[.*?\\]", "", name))
}

# Good: Normalized access
get_track <- function(result, track_name) {
  tracks <- names(result)
  target_normalized <- normalize_track_name(track_name)
  matching_track <- tracks[normalize_track_name(tracks) == target_normalized]
  if (length(matching_track) > 0) {
    return(result[[matching_track[1]]])
  }
  return(NULL)
}

# Usage works with both packages and all naming variations
rms_data <- get_track(result, "rms")      # Works for both "RMS[dB]" and "rms"
rms_data <- get_track(result, "RMS[dB]")  # Also works
acf_data <- get_track(result, "acf")      # Works for both "ACF" and "acf"
```

## Test Suite Implementation

The equivalence test suites (`test-equivalence-comprehensive.R` and `test-equivalence-simple.R`) handle this by:

1. Normalizing track names before comparison (removing bracket notation and converting to lowercase)
2. Matching tracks using normalized names when comparing values
3. Testing that data is equivalent regardless of naming conventions

## Why the Difference?

- **Historical**: Different development teams, different conventions
- **Not a bug**: Both are valid, just different styles
- **Maintained**: Each package maintains its convention for backward compatibility

## Conclusion

These are cosmetic differences that don't affect functionality:
1. **Case differences**: SuperASP uses UPPERCASE, wrassp uses lowercase
2. **Unit notation**: SuperASP includes units in brackets (e.g., [dB], [Hz]), wrassp omits them

The test suite handles these differences transparently using a normalization function that removes bracket notation and converts to lowercase. User code can easily accommodate both conventions with the same approach.
