# Units of Measurement Verification Summary

## Overview

This document summarizes the verification of units of measurement across all DSP functions in the superassp package, with particular attention to Praat/Parselmouth-based voice analysis functions.

## Verification Process

1. Consulted Praat manual (https://www.fon.hum.uva.nl/praat/manual/)
2. Reviewed specific documentation pages for Jitter and Shimmer
3. Examined original Praat scripts in `inst/praat/`
4. Verified R documentation in function roxygen2 comments
5. Updated Python documentation with explicit unit specifications
6. Created comprehensive units documentation (UNITS_OF_MEASUREMENT.md)

## Key Findings and Corrections

### praat_voice_report / praat_voice_report_opt

#### Units Verified as Correct

| Measurement Type | Parameters | Unit | Source |
|-----------------|------------|------|--------|
| Pitch | Median, Mean, SD, Min, Max | Hz | Praat manual |
| Pulse/Period counts | Number of pulses, Number of periods | count | Praat manual |
| Period time | Mean period, SD of period | **seconds** | Praat manual |
| Voicing fractions | Fraction unvoiced, Degree of breaks | fraction (0-1) | Praat manual |
| Jitter (relative) | local, rap, ppq5, ddp | % (percentage) | Praat Voice 2 |
| Jitter (absolute) | local, absolute | **seconds** | Praat Voice 2 |
| Shimmer (relative) | local, apq3, apq5, apq11, dda | % (percentage) | Praat Voice 3 |
| Shimmer (absolute) | local, dB | dB (decibels) | Praat Voice 3 |
| Harmonicity | HNR | dB (decibels) | Praat manual |
| Harmonicity | NHR | ratio (unitless) | Praat manual |
| Autocorrelation | Mean autocorrelation | fraction (0-1) | Praat manual |

#### Critical Note: Jitter (local, absolute)

**Important Discovery**:

- **Literature convention**: Often reported in microseconds (μs)
  - Example: MDVP threshold = 83.200 μs

- **Praat actual output**: Returns values in **seconds** (s)
  - Same threshold: 0.0000832 s

- **Our implementation**: Correctly returns **seconds** to match Praat

**Conversion formula** (if users need microseconds):
```r
jitter_microseconds <- jitter_seconds * 1e6
```

### Documentation Updates Made

#### 1. R Documentation (R/praat_slicefunctions.R)

**Added explicit units**:
- Mean period: "in seconds" (was missing)
- Standard deviation of period: "in seconds" (was missing)
- Mean autocorrelation: "unitless, 0-1" (was vague)
- Mean NHR: "ratio, unitless" (was missing)
- Mean HNR: "in dB" (was missing)

#### 2. Python Documentation (inst/python/praat_voice_report_memory.py)

**Added comprehensive Returns section** with:
- Complete list of all 26 measurements
- Explicit unit for each measurement
- Warning about jitter (local, absolute) being in seconds, not microseconds
- Organized by measurement category
- Notes section explaining unit conventions

#### 3. New Documentation Files Created

**UNITS_OF_MEASUREMENT.md**:
- Complete reference for all voice parameters
- Units table for each measurement category
- Common confusions and clarifications section
- Conversion formulas where relevant
- References to Praat manual pages

**UNITS_VERIFICATION_SUMMARY.md** (this file):
- Summary of verification process
- Key findings
- Documentation updates made

## Verification Checklist

- [x] Consulted Praat manual for authoritative unit specifications
- [x] Verified Jitter units (% for relative, seconds for absolute)
- [x] Verified Shimmer units (% for relative, dB for absolute)
- [x] Verified period measurements are in seconds
- [x] Verified pitch measurements are in Hz
- [x] Verified harmonicity measurements (HNR in dB, NHR as ratio)
- [x] Updated R documentation with explicit units
- [x] Updated Python documentation with comprehensive unit information
- [x] Created reference documentation (UNITS_OF_MEASUREMENT.md)
- [x] Documented critical jitter (local, absolute) seconds vs microseconds issue
- [x] Rebuilt package with updated documentation

## Implementation Compliance

### Parselmouth Functions

**praat_voice_report_memory.py**:
- ✅ Returns values exactly as Praat provides them
- ✅ No unit conversions performed
- ✅ Comprehensive documentation of units in docstring
- ✅ Warning about jitter (local, absolute) convention included

### R Wrapper Functions

**praat_voice_report_opt()**:
- ✅ Passes values unchanged from Python
- ✅ @return documentation specifies units via @inheritParams
- ✅ Function attributes maintain compatibility with original functions

### Original Functions

**praat_voice_report()**:
- ✅ Documentation updated with explicit units
- ✅ Maintains same behavior (values from Praat unchanged)

## Testing Recommendations

When testing equivalence between old and new implementations:

1. **Exact value matching**: Values should match within floating-point precision
2. **Unit consistency**: Both should return same units (no conversions)
3. **Special attention to**:
   - Jitter (local, absolute) should be in seconds (small values like 0.000050)
   - Period measurements should be in seconds (values like 0.005-0.015)
   - Percentage jitter/shimmer should be reasonable (0.5-5.0%)

## References

1. Praat Manual: https://www.fon.hum.uva.nl/praat/manual/Manual.html
2. Praat Voice 2 (Jitter): https://www.fon.hum.uva.nl/praat/manual/Voice_2__Jitter.html
3. Praat Voice 3 (Shimmer): https://www.fon.hum.uva.nl/praat/manual/Voice_3__Shimmer.html
4. Package documentation: UNITS_OF_MEASUREMENT.md

## Conclusion

All units of measurement have been verified against the Praat manual and are correctly documented. The most critical finding is that jitter (local, absolute) is returned in **seconds**, not microseconds as commonly reported in literature. This has been explicitly documented in all relevant locations to prevent confusion.

The documentation now provides complete, unambiguous unit specifications for all voice analysis parameters.
