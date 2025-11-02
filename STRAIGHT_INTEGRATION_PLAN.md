# STRAIGHT Integration Plan - Superassp Package

**Date:** November 1, 2025  
**Goal:** Integrate improved STRAIGHT implementation into superassp package  
**Current Accuracy:** 91.9% frame, 99.0% mean F0  
**Target Accuracy:** ≥91% frame accuracy (maintain current level)

---

## Current State Analysis

### Legacy_STRAIGHT Project (`/Users/frkkan96/Documents/src/legacy_STRAIGHT`)
- **Location:** `straight_python/straight/`
- **Current Accuracy:** 91.9% frame, 99.0% mean F0
- **Files:**
  - `__init__.py` - Module exports (synthesis, spectral, aperiodicity)
  - `f0_extraction.py` - F0 extraction (MulticueF0v14 function) - 92KB, heavily optimized
  - `f0_extraction_optimized.py` - Numba-optimized version
  - `spectral.py` - Spectral analysis (99.9962% accurate)
  - `aperiodicity.py` - Aperiodicity analysis (99.83% accurate)
  - `synthesis.py` - Speech synthesis (99.99% accurate)

### Superassp Package (`/Users/frkkan96/Documents/src/superassp`)
- **Location:** `inst/python/legacy_STRAIGHT/`
- **Current Integration:**
  - Uses older version with lower accuracy (~91% frame, ~96.5% mean F0)
  - R wrappers: `R/ssff_python_straight_f0.R`, `R/ssff_python_straight_spec.R`, `R/ssff_python_straight_synth.R`
  - Installation helper: `R/install_legacy_straight.R`

---

## Integration Strategy

### Phase 1: Backup and Preparation (5 min)
1. ✅ Backup current superassp implementation
2. ✅ Document current function signatures
3. ✅ Identify dependencies

### Phase 2: Code Transfer (15 min)
1. Copy improved Python files from legacy_STRAIGHT to superassp
2. Update `__init__.py` to expose F0 extraction functions
3. Preserve existing API compatibility

### Phase 3: R Wrapper Updates (20 min)
1. Update `R/ssff_python_straight_f0.R` to:
   - Use new function signature
   - Update documentation with new accuracy metrics
   - Maintain backward compatibility
2. Update `R/install_legacy_straight.R` to include Numba option
3. Test R-Python interface

### Phase 4: Testing (30 min)
1. **Unit Tests:**
   - Test F0 extraction with sample files
   - Verify output format matches AsspDataObj spec
   - Test batch processing
   - Test time windowing

2. **Accuracy Tests:**
   - Run against MATLAB baseline
   - Verify ≥91% frame accuracy maintained
   - Check mean F0 accuracy ≥99%

3. **Integration Tests:**
   - Test with devtools::load_all()
   - Test with AVAudio objects
   - Test file output (toFile = TRUE)
   - Test batch processing

### Phase 5: Package Build Validation (20 min)
1. Run `devtools::document()`
2. Run `devtools::test()`
3. Run `devtools::check()`
4. Build package with `devtools::build()`
5. Test installation

### Phase 6: Documentation (15 min)
1. Update function documentation
2. Update CLAUDE.md
3. Update NEWS.md
4. Create migration notes if API changed

---

## Files to Modify

### Python Files (in superassp)
- `inst/python/legacy_STRAIGHT/__init__.py` - Add F0 exports
- `inst/python/legacy_STRAIGHT/f0_extraction.py` - Replace with new version
- `inst/python/legacy_STRAIGHT/f0_extraction_optimized.py` - Add Numba version (optional)

### R Files (in superassp)
- `R/ssff_python_straight_f0.R` - Update wrapper
- `R/install_legacy_straight.R` - Update installation helper
- `tests/testthat/test-straight.R` - Update/add tests

### Documentation Files
- `CLAUDE.md` - Update with new accuracy metrics
- `NEWS.md` - Document improvements
- `man/trk_straight_f0.Rd` - Auto-generated from roxygen2

---

## API Compatibility Check

### Current API (superassp)
```python
# Python side
result = MulticueF0v14(x, fs, f0floor, f0ceil)
# Returns: (f0_values, vuv, aux_data)
```

### New API (legacy_STRAIGHT)
```python
# Python side  
result = MulticueF0v14(x, fs, f0floor=71, f0ceil=800, frame_period=5.0)
# Returns: dict with keys: f0, vuv, if_score, ac_score, times, etc.
```

**Action Required:** Update R wrapper to handle dict return format

---

## Risk Assessment

### Low Risk
- ✅ Spectral and aperiodicity modules unchanged (99%+ accurate)
- ✅ Synthesis module unchanged (99.99% accurate)
- ✅ Python package structure similar

### Medium Risk
- ⚠️ F0 API return format changed (dict vs tuple)
- ⚠️ New dependencies (check if any)
- ⚠️ Numba optional (graceful degradation needed)

### Mitigation
- Thorough API compatibility testing
- Graceful fallback if Numba not available
- Comprehensive test suite before integration

---

## Success Criteria

### Must Have ✅
1. Package builds without errors
2. All existing tests pass
3. F0 accuracy ≥91% frame, ≥99% mean
4. Backward compatibility maintained
5. Documentation updated

### Nice to Have 🎯
1. Numba optimization working
2. Performance benchmarks documented
3. Migration guide for users
4. Additional tests for edge cases

---

## Execution Checklist

- [ ] Phase 1: Backup complete
- [ ] Phase 2: Code transferred
- [ ] Phase 3: R wrappers updated
- [ ] Phase 4: All tests passing
- [ ] Phase 5: Package builds successfully
- [ ] Phase 6: Documentation complete

---

## Notes

- Current implementation already has good accuracy (91.9% / 99.0%)
- Integration should maintain this level
- Focus on clean integration, not further optimization
- Document any breaking changes clearly
- Keep superassp package constraints in mind (from CLAUDE.md)

---

*Plan Created: November 1, 2025*  
*Estimated Time: 105 minutes*  
*Priority: High - Package Integration*
