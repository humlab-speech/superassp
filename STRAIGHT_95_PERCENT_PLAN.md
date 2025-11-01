# STRAIGHT Python Implementation - Path to 95% Accuracy

**Date**: November 1, 2025  
**Current State**: 91.3% frame accuracy (superassp baseline)  
**Target**: >95% frame accuracy  
**Gap**: ~4% improvement needed  

---

## Executive Summary

The legacy STRAIGHT implementation in superassp currently achieves **91.3% frame-level F0 accuracy**. Recent attempts (Jan 30) to reach 99% accuracy through incremental fixes **degraded performance to 83.6%**. Based on lessons learned, achieving 95% requires a systematic, test-driven approach rather than ad-hoc modifications.

---

## Current Status Assessment

### What We Have
| Component | Location | Accuracy | Status |
|-----------|----------|----------|--------|
| superassp baseline | `inst/python/legacy_STRAIGHT/` | 91.3% | ✅ Stable (Oct 29) |
| Modified version | `~/Documents/src/legacy_STRAIGHT/straight_python/` | 83.6% | ❌ Degraded (Jan 30) |
| Backup | `straight_python_modified_backup/` | Unknown | ❓ Untested |

### Known Issues
1. **Segfault in R wrapper** - Current superassp version crashes during testing
   - Occurs after "Step 2/8: AC-based F0 candidate extraction"
   - Needs immediate investigation
   - Possible causes: numpy array handling, memory corruption, parameter mismatch

2. **Octave errors in low F0 regions** (< 100 Hz)
   - Affects ~8-10% of frames
   - Primarily at utterance onset for male speakers
   - Root cause: backward tracking stops too early (frame 121 instead of frame 0)

3. **Modified version degradation**
   - Multiple experimental fixes applied without systematic validation
   - Broke high-F0 tracking while fixing low-F0
   - Demonstrates need for comprehensive test framework

---

## Recommended Approach: Test-Driven Development

### Phase 1: Fix Immediate Issues (1-2 days)

#### Step 1.1: Resolve Segfault
**Priority**: CRITICAL  
**Estimated time**: 2-4 hours

Actions:
1. Run tests under debugger to identify crash location
2. Check numpy array memory management in R wrapper
3. Verify Python function signatures match R calls
4. Test with smaller audio files to isolate issue
5. Compare working version (if backup is good) with current

**Success criteria**: Tests pass without segfault

#### Step 1.2: Establish Baseline
**Priority**: HIGH  
**Estimated time**: 2-4 hours

Actions:
1. Verify 91.3% accuracy claim with validation script
2. Document exact test conditions (audio file, parameters, MATLAB version)
3. Create reference output for regression testing
4. Ensure Numba optimization is working (20% speedup)

**Success criteria**: Reproducible 91.3% baseline established

### Phase 2: Build Test Infrastructure (2-3 days)

#### Step 2.1: Create Comprehensive Test Suite
**Estimated time**: 1 day

Components needed:
1. **Unit tests** for each F0 extraction stage:
   - IF candidate extraction
   - AC candidate extraction
   - Candidate fusion
   - Tracking algorithm
   - Gap filling
   - V/UV decision

2. **Accuracy tests** per frame range:
   - Frames 0-30 (onset)
   - Frames 31-120 (early tracking)
   - Frames 121-700 (stable region)
   - Frames 701+ (ending)

3. **Parameter sweep tests**:
   - F0 floor: 40, 60, 80 Hz
   - F0 ceiling: 400, 600, 800 Hz
   - Frame shift: 1, 2, 5 ms

4. **Multi-speaker validation**:
   - Male speakers (low F0)
   - Female speakers (high F0)
   - Children (very high F0)

**Success criteria**: Test framework reports accuracy per stage and frame range

#### Step 2.2: MATLAB Comparison Framework
**Estimated time**: 1-2 days

Tools needed:
1. **MATLAB bridge** for live comparison:
   ```matlab
   % Export MATLAB intermediate values
   save('matlab_if_candidates.mat', 'if_cands', 'if_scores');
   save('matlab_ac_candidates.mat', 'ac_cands', 'ac_scores');
   save('matlab_fusion.mat', 'fused_f0', 'fused_scores');
   save('matlab_tracking.mat', 'tracked_f0', 'segments');
   ```

2. **Python comparison script**:
   ```python
   def compare_with_matlab(stage_name, python_output, matlab_file):
       matlab_data = scipy.io.loadmat(matlab_file)
       correlation = np.corrcoef(python_output, matlab_data['variable'])[0,1]
       rmse = np.sqrt(np.mean((python_output - matlab_data['variable'])**2))
       return {'correlation': correlation, 'rmse': rmse}
   ```

3. **Automated comparison reports**:
   - Generate HTML reports with plots
   - Highlight discrepancies > threshold
   - Track accuracy metrics over time

**Success criteria**: Can identify exact stage where Python diverges from MATLAB

### Phase 3: Systematic Improvement (1-2 weeks)

#### Step 3.1: Fix Backward Tracking (High Priority)
**Target**: +3-4% accuracy  
**Estimated time**: 3-5 days

**Problem**: Tracking stops at frame 121, should reach frame 0

Current stopping criteria (line ~1920 in f0_extraction.py):
```python
if (best_distance > 0.1 or 
    pwsdb[ii] < noiselevel + 6 or 
    maskr[ii, idx] == 0 or 
    relv[ii, idx] < 0.17):
    break
```

**Proposed fix** (with testing):
```python
# Frame-dependent relaxation for early voiced regions
if ii < 150:
    distance_threshold = 0.15  # More lenient
    power_threshold = noiselevel + 3
    rel_threshold = 0.10
else:
    distance_threshold = 0.1
    power_threshold = noiselevel + 6
    rel_threshold = 0.17

if (best_distance > distance_threshold or 
    pwsdb[ii] < power_threshold or 
    maskr[ii, idx] == 0 or 
    relv[ii, idx] < rel_threshold):
    break
```

**Testing approach**:
1. Test on 10 different speakers
2. Measure accuracy improvement per frame range
3. Ensure no degradation in frames 121-700
4. Validate that tracking reaches frame 0
5. Check for over-relaxation (spurious tracking)

**Success criteria**: 
- Tracking reaches frame 0 for appropriate speech
- Accuracy in frames 0-120 improves by 2-3%
- No degradation in other regions

#### Step 3.2: Improve AC Candidate Weighting
**Target**: +1-2% accuracy  
**Estimated time**: 2-3 days

**Problem**: IF candidates dominate over AC in fusion, causing octave errors

**Approach**:
1. Analyze IF vs AC reliability in different F0 ranges
2. Test frequency-dependent weighting:
   ```python
   if f0_candidate < 100:
       ac_weight = 0.7  # Favor AC for low F0
       if_weight = 0.3
   elif f0_candidate < 200:
       ac_weight = 0.5  # Balanced
       if_weight = 0.5
   else:
       ac_weight = 0.4  # Favor IF for high F0
       if_weight = 0.6
   ```

3. Validate across speakers and F0 ranges
4. A/B test different weighting schemes

**Testing metrics**:
- Octave error rate (f0_python / f0_matlab ratio)
- Accuracy in 50-100 Hz range
- Accuracy in 100-300 Hz range

**Success criteria**:
- Octave errors reduced by 50%
- Overall accuracy improves by 1-2%

#### Step 3.3: Refine Gap Filling (Low Priority)
**Target**: +0.5-1% accuracy  
**Estimated time**: 1-2 days

**Current approach**: Linear interpolation or candidate-based filling

**Improvements to test**:
1. Use AC candidates for initial gap (before first anchor)
2. Preserve octave consistency
3. Smooth transitions at segment boundaries

**Success criteria**:
- Gaps filled with appropriate F0
- No artifacts at segment boundaries

### Phase 4: Validation and Documentation (2-3 days)

#### Step 4.1: Multi-Speaker Validation
**Datasets to test**:
- VaiUEO (Japanese vowels) - current test
- TIMIT (English speech) - diverse speakers
- Arctic (voice conversion corpus)
- Custom male/female/child samples

**Metrics to report**:
- Mean accuracy per speaker
- Accuracy by F0 range
- Octave error rate
- V/UV decision accuracy

#### Step 4.2: Performance Optimization
- Profile with Numba JIT
- Identify bottlenecks
- Optimize critical loops
- Target: < 1.0x real-time

#### Step 4.3: Documentation Updates
1. Update ACCURACY_REPORT.md with new results
2. Document parameter tuning recommendations
3. Create usage examples for different speech types
4. Update CLAUDE.md with implementation details

---

## Timeline Estimate

| Phase | Duration | Deliverable |
|-------|----------|-------------|
| Phase 1: Fix Issues | 1-2 days | Working baseline, no segfaults |
| Phase 2: Test Infrastructure | 2-3 days | Comprehensive test suite |
| Phase 3: Improvements | 1-2 weeks | 95% accuracy achieved |
| Phase 4: Validation | 2-3 days | Multi-speaker validation, docs |
| **Total** | **2-3 weeks** | Production-ready 95% implementation |

---

## Risk Assessment

### Low Risk (High Confidence)
✅ Fixing backward tracking to reach frame 0  
✅ Establishing test infrastructure  
✅ Multi-speaker validation  

### Medium Risk (Moderate Confidence)
⚠️ AC/IF weighting improvements  
⚠️ Achieving exactly 95% (might be 93-97%)  
⚠️ Maintaining performance with changes  

### High Risk (Low Confidence)
❌ Reaching 99% accuracy (not recommended)  
❌ Perfect MATLAB replication  
❌ Zero octave errors  

---

## Alternative: Hybrid Approach

If systematic improvement proves too time-consuming, consider:

**Use different algorithms for different tasks**:
```r
# High-accuracy F0 for analysis
f0_accurate <- trk_rapt(audio, toFile = FALSE)  # >98% accuracy

# STRAIGHT spectral analysis (99.996% accuracy)
spec_straight <- trk_straight_spec(audio, toFile = FALSE)

# STRAIGHT synthesis with external F0
audio_synth <- straight_synth(
  f0 = f0_accurate$f0[,1],
  spec = spec_straight$spec,
  sample_rate = 22050
)
```

**Benefits**:
- Immediate >98% F0 accuracy
- Best-in-class spectral analysis
- High-quality synthesis
- No development time required

**Trade-offs**:
- Not a pure STRAIGHT implementation
- Two separate F0 algorithms to maintain

---

## Immediate Next Steps

### TODAY:
1. ✅ Fix segfault in superassp R wrapper
2. ✅ Verify 91.3% baseline
3. ✅ Run existing test suite

### THIS WEEK:
1. Build comprehensive test framework
2. Set up MATLAB comparison infrastructure
3. Begin systematic backward tracking fixes

### DECISION POINT:
After fixing segfault and establishing baseline, decide:
- **Option A**: Proceed with systematic improvement (2-3 weeks to 95%)
- **Option B**: Use hybrid approach (RAPT + STRAIGHT spectral)
- **Option C**: Accept 91.3% as sufficient for use case

---

## Success Criteria

### Minimum Viable (Must Have):
- ✅ No segfaults or crashes
- ✅ 91.3% baseline reproducible
- ✅ All tests pass
- ✅ Package builds successfully

### Target (Should Have):
- ✅ 95% frame accuracy
- ✅ < 5% octave error rate
- ✅ Multi-speaker validation
- ✅ < 1.0x real-time performance

### Stretch (Nice to Have):
- ✅ 97%+ accuracy
- ✅ < 2% octave error rate
- ✅ Numba optimization functional
- ✅ Comprehensive documentation

---

## Conclusion

The path to 95% accuracy is clear and achievable through systematic, test-driven development. The key lessons from the Jan 30 session are:
1. **Test incrementally** - validate each change before moving on
2. **Use MATLAB as ground truth** - compare at each stage
3. **Don't break what works** - regression tests are critical
4. **Frame-specific analysis** - understand errors by region
5. **Document everything** - for reproducibility

**Recommendation**: Fix the segfault first, then decide whether to invest 2-3 weeks in systematic improvement or use the hybrid approach for immediate results.

The decision should be based on:
- **Time available**: Do we have 2-3 weeks?
- **Accuracy requirements**: Is 91.3% sufficient? Is 98% needed?
- **Use case**: Research (needs perfect F0) vs Production (needs good synthesis)?

---

*This plan is based on the Jan 30 session findings and represents a more realistic, systematic approach to achieving the 95% accuracy target.*
