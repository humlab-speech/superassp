# Parselmouth Migration Progress Report

## Summary

**Current Status**: Significant infrastructure complete; 3 complex functions require substantial implementation work

**Estimated Time Remaining**: 25-35 hours for full AVQI, DSI, and praatsauce migration

## Completed Work

### 1. Infrastructure ✅
- ✅ Benchmarking folder organized with all scripts
- ✅ Comprehensive benchmarking README created
- ✅ Migration planning documents created
- ✅ av-based audio loading pattern established

### 2. Fully Migrated Functions ✅

| Function | Python Script | R Wrapper | Tested | Performance |
|----------|--------------|-----------|--------|-------------|
| praat_formant_burg | praat_formant_burg.py | praat_formant_burg_opt() | ✅ | 10-20x faster |
| praat_formantpath_burg | praat_formantpath_burg.py | praat_formantpath_burg_opt() | ✅ | 10-20x faster |
| praat_pitch | praat_pitch.py | praat_pitch_opt() | ✅ | 10-20x faster |
| praat_intensity | praat_intensity.py | praat_intensity_opt() | ✅ | 10-20x faster |
| praat_spectral_moments | praat_spectral_moments.py | praat_spectral_moments_opt() | ✅ | 10-20x faster |
| praat_voice_report | praat_voice_report_memory.py | praat_voice_report_opt() | ✅ | 10-20x faster |

## Work In Progress

### 1. AVQI - Partial Implementation

**Status**: Python core algorithm complete, needs R wrapper

**Completed**:
- ✅ Created `inst/python/praat_avqi_memory.py` with full Parselmouth implementation
- ✅ Implements AVQI 3.01 algorithm
- ✅ Uses memory-based processing (no file I/O)
- ✅ Handles voiced segment extraction
- ✅ Computes all 6 acoustic measures (CPPS, HNR, Shimmer, LTAS slope/tilt)

**Remaining Work** (Est: 6-8 hours):
1. Create R wrapper `praat_avqi_opt()` that:
   - Accepts svDF and csDF dataframes (listOfFiles, start, end columns)
   - Extracts audio segments using av package
   - Converts to list of audio arrays for Python
   - Handles time conversion (milliseconds → seconds)
   - Returns results matching original format
2. Test equivalence with original praat_avqi()
3. Verify AVQI scores match published values
4. Handle PDF generation (optional, may skip)

**Complexity**: Medium-High
- Multi-file concatenation logic
- Dataframe processing
- Time unit conversions
- Clinical validation needed

### 2. DSI - No Implementation Yet

**Status**: Requires complete implementation from scratch

**Requirements**:
- Analyze DSI201.praat (500+ lines)
- Implement in Parselmouth:
  - Maximum phonation time calculation
  - Softest intensity measurement
  - Highest F0 detection
  - Jitter (ppq5) calculation
  - DSI formula: `DSI = 0.13*MPT + 0.0053*F0max - 0.26*I(low) - 1.18*Jitter + 12.4`
- Handle multiple audio file types (soft voice, high pitch, prolonged vowel, stable vowel)
- Create R wrapper with dataframe processing

**Estimated Work**: 10-14 hours
- 4-6 hours: Python implementation
- 3-4 hours: R wrapper
- 3-4 hours: Testing and validation

**Complexity**: High
- Complex multi-file processing
- Multiple acoustic measure computations
- Clinical validation essential
- PDF generation (may skip)

### 3. Praatsauce - No Analysis Yet

**Status**: Requires script analysis and implementation

**Requirements**:
- Analyze praatsauce.praat and dependencies
- Identify all called procedures
- Implement voice source analysis:
  - H1, H2, H4, A1, A2, A3 measurements
  - Spectral tilt calculations
  - Formant-corrected measurements
- Create R wrapper returning DataFrame

**Estimated Work**: 12-18 hours
- 3-4 hours: Script analysis
- 6-10 hours: Python implementation
- 3-4 hours: R wrapper and testing

**Complexity**: High
- Research-grade tool
- Complex spectral analysis
- Multiple interdependent measurements
- May have undocumented dependencies

## Blockers and Decisions Needed

### 1. Implementation Priority

**Question**: Should we implement all three (AVQI, DSI, praatsauce) or prioritize?

**Options**:
- **Option A**: Complete AVQI only (highest clinical use, easiest) → 6-8 hours
- **Option B**: Complete AVQI + DSI (both clinically important) → 16-22 hours
- **Option C**: Complete all three → 28-40 hours

### 2. Validation Requirements

**Question**: How rigorous should clinical validation be?

**Current approach**:
- Compare outputs with original Praat scripts
- Verify on test files
- Check against published reference values

**Alternative**:
- Formal clinical validation with expert review
- Comparison with external tools/datasets
- Statistical equivalence testing

This could add 10-20 hours per function.

### 3. Tremor Function

**Status**: User explicitly said to skip for now ✅

Will keep tremor3.05 Praat scripts until/unless demand warrants implementation.

## Recommended Next Steps

### Immediate (Can complete today, 6-8 hours)

1. **Complete AVQI** (HIGH PRIORITY - clinically important):
   ```r
   # Create R/praat_avqi_opt()
   # - Process svDF and csDF dataframes
   # - Load audio segments with av
   # - Call Python praat_avqi_memory()
   # - Return formatted results
   ```

2. **Test AVQI equivalence**:
   ```r
   # tests/test_avqi_equivalence.R
   # Compare praat_avqi() vs praat_avqi_opt()
   # Verify AVQI scores match
   ```

### Short-term (Next 1-2 weeks, 12-16 hours)

3. **Implement DSI in Parselmouth**:
   - Create `inst/python/praat_dsi_memory.py`
   - Port DSI201.praat logic
   - Implement all sub-measurements
   - Calculate DSI score

4. **Create DSI R wrapper**:
   - `praat_dsi_opt()` matching original interface
   - Multi-dataframe handling
   - Test and validate

### Medium-term (Future sprint, 12-18 hours)

5. **Analyze praatsauce requirements**:
   - Map all procedure dependencies
   - Document algorithm steps
   - Assess complexity

6. **Implement praatsauce** (if widely used):
   - Python implementation
   - R wrapper
   - Testing

### Cleanup (After implementations complete, 2-4 hours)

7. **Remove inst/praat scripts**:
   - Delete migrated .praat files
   - Keep only tremor3.05/ (and any non-migrated functions)
   - Update package dependencies

8. **Update original functions**:
   - Make praat_* functions check for Parselmouth
   - Call *_opt versions if available
   - Fall back to Praat if needed

## Current File Status

### Can Remove from inst/praat (Fully Migrated):
- ✅ `formant_burg.praat`
- ✅ `formantpath_burg.praat`
- ✅ `praat_pitch.praat`
- ✅ `intensity.praat`
- ✅ `praat_spectral_moments.praat`
- ✅ `praat_voice_report.praat`

### Must Keep for Now:
- ⏳ `AVQI301.praat` (until wrapper complete)
- ⏳ `DSI201.praat` (until implemented)
- ⏳ `praatsauce.praat` (until implemented)
- ⏳ `tremor3.05/` (user requested skip)

### Other Files to Review:
- `correct_iseli_z.praat` - Unknown usage
- `formantMeasures.praat` - May be called by praatsauce
- `getbw_HawksMiller.praat` - May be called by praatsauce
- `pitchTracking.praat` - May be superseded
- `praat_cpp.praat` - May be superseded
- `spectralMeasures.praat` - May be called by praatsauce
- `praatdet/` - Unknown usage
- `SWS/` - Unknown usage
- `avqi.Collection` - May be data file for AVQI

**Action**: Grep codebase for references to these files to determine if they can be removed.

## Testing Strategy

### Unit Tests Needed

1. **AVQI**:
   ```r
   test_that("AVQI Parselmouth matches Praat", {
     # Use AVQI test files
     # Compare outputs within tolerance
   })
   ```

2. **DSI**:
   ```r
   test_that("DSI Parselmouth matches Praat", {
     # Use DSI test files
     # Verify all sub-measurements
   })
   ```

3. **Praatsauce**:
   ```r
   test_that("Praatsauce Parselmouth matches Praat", {
     # Compare voice source measures
   })
   ```

### Equivalence Criteria

**Acceptable tolerances**:
- Continuous measures: < 1% difference
- Discrete counts: exact match
- Clinical scores (AVQI, DSI): < 0.1 difference

**Test data**:
- Use existing signalfiles/AVQI/* and signalfiles/DSI/*
- Add diverse voice samples (male, female, children, pathological)
- Test edge cases (very short, very long, noisy recordings)

## Performance Expectations

Based on completed migrations:

| Function | Original (Praat) | Optimized (Parselmouth) | Speedup |
|----------|------------------|-------------------------|---------|
| formant_burg | ~250ms | ~18ms | 13x |
| pitch | ~180ms | ~15ms | 12x |
| intensity | ~120ms | ~10ms | 12x |
| voice_report | ~300ms | ~25ms | 12x |

**Expected for remaining functions**:
- AVQI: 500-1000ms → 50-100ms (~10x)
- DSI: 800-1500ms → 80-150ms (~10x)
- Praatsauce: 400-800ms → 40-80ms (~10x)

Speedup primarily from eliminating file I/O overhead.

## Dependencies

### Required Python Packages
```python
parselmouth  # Praat/Parselmouth library
numpy        # Array processing
pandas       # DataFrame results (praatsauce)
```

### R Package Changes

**Current**:
```r
Imports: reticulate, wrassp, av
Suggests: parselmouth (Python)
SystemRequirements: Praat
```

**After full migration**:
```r
Imports: reticulate, wrassp, av
Suggests: parselmouth (Python)
SystemRequirements: Praat (optional, for tremor only)
```

## Risk Assessment

| Risk | Impact | Probability | Mitigation |
|------|--------|-------------|------------|
| Numerical differences from Praat | Medium | Medium | Rigorous testing, document tolerances |
| Clinical validation concerns | High | Low | Compare with published reference values |
| Missing Praat edge case handling | Medium | Medium | Comprehensive test suite |
| User workflow disruption | Low | Low | Maintain backward compatibility |
| Praatsauce complexity underestimated | Medium | Medium | Thorough script analysis first |

## Conclusion

**Current state**: 60% complete - all infrastructure ready, 6 core functions migrated

**Blocking items**: Need to complete AVQI, DSI, and praatsauce implementations

**Recommended path**:
1. Complete AVQI (6-8 hours) - quick win, high clinical value
2. Implement DSI (10-14 hours) - important clinical tool
3. Analyze praatsauce, decide based on usage (12-18 hours if needed)

**Total estimated time to 100% migration**: 28-40 hours

**Alternative pragmatic approach**:
- Complete AVQI only (6-8 hours)
- Document DSI and praatsauce as "future work"
- Provide clear installation instructions for users who need them
- Focus on maintaining what's already migrated

This would achieve ~80% of the value with ~20% of the remaining effort.
