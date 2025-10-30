# Migration Status Update - High Priority Functions

**Date**: 2025-10-30
**Context**: Post-audit migration request

## Status Check Results

After reviewing the 5 high-priority functions identified in the audit, here's the actual status:

### ✅ Already Migrated (1 function)

#### trk_yaapt
- **File**: `R/ssff_python_yaapt.R`
- **Status**: ✅ **ALREADY USES AV PACKAGE**
- **Implementation**: Lines 134-151 use `av::read_audio_bin()`
- **In-memory**: ✅ Yes
- **Action**: None needed - already compliant

### ⚠️ Needs Migration (4 functions)

#### 1. trk_snackp (Snack Pitch Tracking)
- **File**: `R/ssff_python_snack_pitch.R`
- **Current**: Calls external Python script via `system()` command
- **Issue**: Python script loads audio from file path
- **Python Script**: `inst/python/snack_pitch.py`
- **Migration Pattern**: Update Python script to accept numpy array instead of file path
- **Complexity**: Medium (requires Python script modification)

#### 2. trk_snackf (Snack Formant Tracking)
- **File**: `R/ssff_python_snack_formant.R`
- **Current**: Similar to trk_snackp, uses external Python script
- **Python Script**: `inst/python/snack_formant.py` (assumed)
- **Migration Pattern**: Update Python script to accept numpy array
- **Complexity**: Medium

#### 3. trk_formantp (Praat/Parselmouth Formant Tracking)
- **File**: `R/ssff_python_pm_pformantb.R`
- **Current**: Uses temporary files with Parselmouth
- **Migration Pattern**: Implement `av_load_for_parselmouth()` (Pattern 2)
- **Complexity**: Low-Medium (helper function exists)

#### 4. trk_formantpathp (Praat/Parselmouth Formant Path)
- **File**: `R/ssff_python_pm_pformantpathb.R`
- **Current**: Uses temporary files with Parselmouth
- **Migration Pattern**: Implement `av_load_for_parselmouth()` (Pattern 2)
- **Complexity**: Low-Medium

---

## Revised Migration Plan

### Complexity Assessment

| Function | Lines of Code | Python Script | Complexity | Est. Time |
|----------|---------------|---------------|------------|-----------|
| trk_yaapt | N/A | N/A | ✅ Done | 0 hours |
| trk_snackp | ~250 | Yes (external) | Medium | 2-3 hours |
| trk_snackf | ~250 | Yes (external) | Medium | 2-3 hours |
| trk_formantp | ~200 | No (inline) | Low-Med | 1-2 hours |
| trk_formantpathp | ~200 | No (inline) | Low-Med | 1-2 hours |

**Total Estimated Time**: 6-10 hours

### Migration Challenges

#### Snack Functions (trk_snackp, trk_snackf)

**Challenge**: These use external Python scripts called via `system()`:
```r
cmd <- sprintf("python3 '%s' '%s'", python_script, params_json)
result_json <- system(cmd, intern = TRUE)
```

**Options**:
1. **Rewrite to use reticulate** (Recommended)
   - Load audio with `av::read_audio_bin()`
   - Pass numpy array to Python via reticulate
   - Modify Python functions to accept arrays instead of file paths

2. **Keep external script but pre-load audio**
   - Load audio with av
   - Save to temp WAV
   - Pass temp path to script
   - ❌ Not ideal - still uses temp files

3. **Reimplement in R or C++**
   - Port Snack algorithms to native code
   - ❌ Too complex for this sprint

**Recommendation**: Option 1 (Rewrite to use reticulate)

#### Parselmouth Functions (trk_formantp, trk_formantpathp)

**Challenge**: Currently use temp files for audio

**Solution**: Use `av_load_for_parselmouth()` helper
```r
# OLD (temp files)
temp_file <- tempfile(fileext = ".wav")
av::av_audio_convert(file_path, temp_file)
sound <- pm$Sound(temp_file)
unlink(temp_file)

# NEW (in-memory)
sound <- av_load_for_parselmouth(
  file_path = file_path,
  start_time = if (beginTime > 0) beginTime else NULL,
  end_time = if (endTime > 0) endTime else NULL,
  channels = 1
)
```

**Recommendation**: Straightforward - use existing helper

---

## Immediate Actions

Given the scope, I recommend:

### Option A: Complete All 4 Migrations (6-10 hours)
**Pros**: Achieves 100% compliance for high-priority list
**Cons**: Requires significant time and testing

**Steps**:
1. Migrate trk_formantp (1-2 hours)
2. Migrate trk_formantpathp (1-2 hours)
3. Migrate trk_snackp (2-3 hours)
4. Migrate trk_snackf (2-3 hours)
5. Test all functions (1-2 hours)
6. Document and commit

### Option B: Migrate Parselmouth Functions Only (2-4 hours)
**Pros**: Quick wins, straightforward implementation
**Cons**: Leaves Snack functions for later

**Steps**:
1. Migrate trk_formantp
2. Migrate trk_formantpathp
3. Test both functions
4. Document and commit
5. Defer Snack functions to v0.9.0

### Option C: Document Current State, Plan Next Steps
**Pros**: Accurate assessment, clear roadmap
**Cons**: No immediate migration progress

**Steps**:
1. Update audit with corrected status
2. Create detailed migration specs
3. Schedule migrations for v0.9.0

---

## Recommendation

I recommend **Option B: Migrate Parselmouth Functions Only** for this session:

**Rationale**:
1. **Low complexity**: Helper function already exists
2. **High impact**: Parselmouth functions are actively used
3. **Quick completion**: 2-4 hours vs 6-10 hours
4. **Testable**: Existing test files available
5. **Defer complexity**: Snack migrations need more planning

**Next Session**:
- Tackle Snack functions with full Python script refactoring
- Estimated 4-6 hours with testing

---

## Files to Modify (Option B)

### R Files
1. `R/ssff_python_pm_pformantb.R` → Update to use `av_load_for_parselmouth()`
2. `R/ssff_python_pm_pformantpathb.R` → Update to use `av_load_for_parselmouth()`

### Test Files
1. Create `tests/testthat/test-formantp-memory.R`
2. Create `tests/testthat/test-formantpathp-memory.R`

### Documentation
1. Update function documentation (remove temp file mentions)
2. Update NEWS.md
3. Create migration completion report

---

## Decision Required

**Question**: Should we:
- **A**: Complete all 4 migrations (6-10 hours)
- **B**: Complete Parselmouth migrations only (2-4 hours)
- **C**: Document and plan for next session

**My recommendation**: **Option B** - Complete the Parselmouth migrations now, defer Snack functions to a dedicated session with proper Python refactoring.

---

## Files Created

- `MIGRATION_STATUS_UPDATE.md` (this document)

**Awaiting decision before proceeding...**
