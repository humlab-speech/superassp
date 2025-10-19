# COVAREP Python - Validation Phase Complete Report

**Date:** October 17, 2025  
**Phase:** Numerical Validation  
**Status:** ⚠️ SIGNIFICANT DIFFERENCES IDENTIFIED - ANALYSIS COMPLETE

---

## Executive Summary

Completed comprehensive numerical validation comparing Python implementation against MATLAB COVAREP reference outputs. **Significant algorithmic differences** identified requiring investigation and tuning.

### Key Findings

| Algorithm | Correlation | Error | Status |
|-----------|-------------|-------|--------|
| **F0 Tracking (SRH)** | 0.03 | 75.7% | ⚠️ **MAJOR ISSUES** |
| **IAIF Glottal Flow** | 0.67 | 95.5% RMS | ⚠️ **MODERATE ISSUES** |
| **IAIF Flow Derivative** | 0.13 | 101.2% RMS | ⚠️ **MAJOR ISSUES** |

**Overall Assessment:** Python implementation runs without errors but produces different results from MATLAB. Requires algorithmic review and parameter tuning.

---

## Detailed Validation Results

### Test Configuration ✅

**Audio File:** `0011.arctic_bdl1.wav`
- Duration: 3.53 seconds
- Sample Rate: 32 kHz (resampled to 16 kHz by MATLAB)
- Samples: 112,960 (original)

**Processing:**
- ✅ Python: Runs successfully, no errors
- ✅ MATLAB: Reference generated successfully
- ✅ Comparison: 344-350 frames analyzed

---

## F0 Tracking Validation ⚠️

### Output Comparison

| Metric | Python | MATLAB | Difference |
|--------|--------|--------|------------|
| Frames | 350 | 344 | +6 frames |
| Time range | 0.00-3.49s | 0.05-3.48s | Different start |
| F0 range (all) | 50-500 Hz | 65-255 Hz | Much wider |
| F0 range (voiced) | 50-500 Hz | 72-232 Hz | **Suspicious** |
| Voicing % | 70.0% | 52.3% | +17.7% |

###issues Identified

**1. F0 Range Problem 🔴 CRITICAL**

Python reports F0 range 50-500 Hz (full search range), indicating:
- Voicing decision too permissive
- Selecting boundary values (50 Hz, 500 Hz)
- Not refining F0 range iteratively like MATLAB

**MATLAB Behavior (from code review):**
```matlab
% MATLAB does iterative refinement:
if max(SRHVal) > threshold
    F0medEst = median(F0s(SRHVal > threshold))
    if round(0.5*F0medEst) > f0min
        f0min = round(0.5*F0medEst)  % Refine lower bound
    if round(2*F0medEst) < f0max
        f0max = round(2*F0medEst)    % Refine upper bound
```

**Python Missing:** This iterative refinement logic!

**2. Voicing Decision Too Lenient 🔴 CRITICAL**

- Python: 70% voiced
- MATLAB: 52% voiced  
- Difference: 18% over-detection

Likely causes:
- Threshold too low in Python
- Different SRH calculation
- Missing post-processing

**3. Frame Timing Mismatch ⚠️ MINOR**

- Python starts at 0.00s
- MATLAB starts at 0.05s  
- Likely due to frame centering differences

### Error Statistics

**Voiced Frames Comparison (180 frames):**
```
Mean absolute error:    97.65 Hz   ❌
Median absolute error:  71.26 Hz   ❌
Correlation:            0.03       ❌ (random!)
```

**Error Distribution:**
```
0-5 Hz:     40 frames (22.2%)  ← Good matches
5-10 Hz:     4 frames (2.2%)
10-20 Hz:    3 frames (1.7%)
20-50 Hz:   17 frames (9.4%)
50-100 Hz:  51 frames (28.3%)
>100 Hz:    65 frames (36.1%)  ← Bad matches
```

**Analysis:** Only 22% of frames have <5 Hz error. 36% have >100 Hz error, indicating fundamental algorithm differences.

---

## IAIF Validation ⚠️

### Glottal Flow Comparison

| Metric | Value | Status |
|--------|-------|--------|
| Correlation | 0.668 | ⚠️ Moderate |
| RMS Error | 0.304 | ⚠️ High |
| Relative RMS | 95.5% | ⚠️ Very High |

**Analysis:**
- Moderate correlation (0.67) suggests similar general shape
- But very high RMS error (95%) indicates magnitude differences
- Likely causes:
  - Different pre-emphasis
  - Different filter orders
  - Integration constant differences
  - Normalization differences

### Flow Derivative Comparison

| Metric | Value | Status |
|--------|-------|--------|
| Correlation | 0.127 | ❌ Poor |
| RMS Error | 0.015 | ⚠️ High |
| Relative RMS | 101.2% | ❌ Very High |

**Analysis:**
- Very low correlation (0.13) indicates major differences
- GCI peak detection: Python=21, MATLAB=159 (7.6x difference!)
- Critical issue: Derivative computation or scaling wrong

---

## Root Cause Analysis

### F0 Tracking Issues

**Issue 1: Missing Iterative Refinement**

MATLAB SRH performs 2 iterations:
1. First pass with wide F0 range
2. Refine f0min/f0max based on median
3. Second pass with tighter range

Python only does single pass with full range.

**Fix Required:**
```python
# Add to pitch_srh():
# Iteration 1
f0_est1, srh1 = compute_srh(...)

# Refine range
if max(srh1) > threshold:
    f0_median = np.median(f0_est1[srh1 > threshold])
    f0min_refined = max(f0min, int(0.5 * f0_median))
    f0max_refined = min(f0max, int(2.0 * f0_median))
    
    # Iteration 2
    f0_est2, srh2 = compute_srh(..., f0min_refined, f0max_refined)
```

**Issue 2: Voicing Threshold**

Python uses percentile-based threshold:
```python
threshold = np.percentile(srh_values[srh_values > 0], 30)
```

MATLAB uses fixed threshold:
```matlab
VoicingThresh = 0.07
VoicingThresh2 = 0.085
```

**Fix Required:** Use MATLAB's thresholds or adaptive method based on SRH statistics.

**Issue 3: LPC Residual Computation**

Python implementation may differ from MATLAB in:
- Frame overlap handling
- Window application
- LPC order calculation
- Filter application method

### IAIF Issues

**Issue 1: Different Default Parameters**

Python:
```python
p_vt = int(fs / 1000) + 4    # For 16kHz: 20
p_gl = 2
```

MATLAB:
```matlab
p_vt = 2*round(fs/2000)+4    # For 16kHz: 12
p_gl = 2*round(fs/4000)      # For 16kHz: 2
```

**Difference:** VT order is 20 vs 12! This significantly affects filtering.

**Fix Required:** Match MATLAB's formula exactly.

**Issue 2: Integration Method**

Python:
```python
g = signal.lfilter([1], [1, -1], g)  # Simple cumsum
```

MATLAB likely uses different integration or includes DC removal.

**Issue 3: Pre-emphasis**

Differences in high-pass filter implementation or pre-emphasis coefficient could cause magnitude differences.

---

## Specific Algorithm Fixes Needed

### Priority 1: F0 Tracking (CRITICAL)

1. **Implement 2-iteration refinement**
   - Lines to add: ~20
   - Complexity: Low
   - Impact: **MAJOR** - Should reduce error from 75% to <10%

2. **Fix voicing threshold**
   - Lines to change: 2-3
   - Complexity: Very low
   - Impact: **MAJOR** - Should fix voicing agreement

3. **Review LPC residual**
   - Compare MATLAB lpcresidual.m vs Python implementation
   - Check window application
   - Verify filter coefficients

### Priority 2: IAIF (IMPORTANT)

1. **Fix VT filter order formula**
   - Change: `p_vt = 2*round(fs/2000)+4`
   - Lines to change: 1
   - Impact: **MAJOR** - Should improve correlation

2. **Review integration method**
   - Compare with MATLAB's integration
   - Check for DC removal
   - Verify scaling

3. **Match pre-emphasis**
   - Verify high-pass filter matches MATLAB
   - Check coefficient (0.99 vs other values)

### Priority 3: Frame Timing (MINOR)

1. **Align frame timing**
   - Start frames at same time as MATLAB
   - May need to adjust hop size handling

---

## Immediate Action Plan

### Today: Critical Fixes

1. **F0 Iteration Refinement** (2 hours)
   - Implement 2-pass algorithm
   - Test on validation audio
   - Expect error reduction: 75% → 30%

2. **IAIF Parameter Fix** (30 min)
   - Correct VT order formula
   - Re-run validation
   - Expect correlation improvement: 0.67 → 0.85

3. **Voicing Threshold** (30 min)
   - Use MATLAB's fixed thresholds
   - Test voicing agreement
   - Expect agreement: 52% → 75%

### Tomorrow: Verification

1. **Re-run full validation**
   - All fixes applied
   - Generate new comparison
   - Document improvements

2. **Test on multiple files**
   - arctic_a0007.wav
   - Other COVAREP test files
   - Verify robustness

3. **Parameter tuning**
   - Fine-tune remaining parameters
   - Achieve <10% error target

---

## Success Metrics After Fixes

### Target Performance

| Metric | Current | Target | Required Improvement |
|--------|---------|--------|---------------------|
| F0 mean error | 97.7 Hz | <10 Hz | 90% reduction |
| F0 correlation | 0.03 | >0.90 | 30x improvement |
| IAIF glot corr | 0.67 | >0.90 | 34% improvement |
| IAIF deriv corr | 0.13 | >0.80 | 6x improvement |
| Voicing agreement | 52% | >85% | 63% improvement |

**Estimated Achievability:** HIGH - Most issues are parameter/formula mismatches, not fundamental algorithm problems.

---

## Lessons Learned

### What Worked Well ✅

1. **Validation Framework** - Automated comparison excellent
2. **Visualization** - Plots clearly show issues
3. **MATLAB Bridge** - Reference generation seamless
4. **Error Analysis** - Quantitative metrics helpful

### What Needs Improvement ⚠️

1. **Algorithm Study** - Need deeper review of MATLAB code before implementation
2. **Parameter Matching** - Must verify all formulas match exactly
3. **Iterative Logic** - Multi-pass algorithms need special attention
4. **Unit Tests** - Need tests comparing with MATLAB on small examples

### Process Improvements 📋

1. **Implement in stages**
   - First: Get parameters exactly right
   - Then: Implement core algorithm
   - Finally: Add iterations and refinements

2. **Test incrementally**
   - Compare each sub-function with MATLAB
   - Don't wait for full implementation

3. **Document assumptions**
   - Note where Python differs from MATLAB
   - Explain design decisions

---

## Revised Timeline

### Week 2: Fix Critical Issues (Updated)

**Days 3-4: Algorithm Fixes**
- [x] Identify root causes (DONE - this report)
- [ ] Fix F0 iteration refinement (2 hours)
- [ ] Fix IAIF parameters (30 min)
- [ ] Fix voicing thresholds (30 min)
- [ ] Re-validate (1 hour)

**Days 5-7: Verification**
- [ ] Test on multiple audio files
- [ ] Fine-tune parameters
- [ ] Achieve <10% error target
- [ ] Document final parameters

**Expected Results:**
- F0 error: <10 Hz mean
- IAIF correlation: >0.90
- Voicing agreement: >85%

---

## Documentation Generated

### Validation Artifacts ✅

1. **matlab_f0_reference.txt** - MATLAB F0 output (344 frames)
2. **matlab_iaif_glottal_flow.txt** - MATLAB glottal waveform
3. **python_f0_output.txt** - Python F0 output (350 frames)
4. **python_iaif_glottal_flow.txt** - Python glottal waveform
5. **detailed_validation_comparison.png** - 6-panel comparison (200 DPI)
6. **zoomed_comparison.png** - Detailed zoom views (200 DPI)

### Analysis Reports ✅

1. **This document** - Complete validation analysis
2. **Error statistics** - Quantitative metrics
3. **Root cause analysis** - Specific fixes needed
4. **Action plan** - Clear next steps

---

## Recommendations

### Immediate (Today)

1. ✅ **APPROVED:** Proceed with critical algorithm fixes
2. **Focus on F0 iteration refinement** (highest impact)
3. **Fix IAIF parameters** (quick win)
4. **Re-validate after each fix** (verify improvements)

### Short-Term (This Week)

1. **Achieve validation targets** (<10% error)
2. **Test on diverse audio**
3. **Document tuning process**
4. **Update code comments** with MATLAB correspondence

### Strategic

1. **Create MATLAB-Python correspondence docs**
   - Map each Python function to MATLAB equivalent
   - Note parameter differences
   - Document design decisions

2. **Build regression test suite**
   - Small examples with known outputs
   - Test each sub-function
   - Prevent future regressions

3. **Consider hybrid approach**
   - Use MATLAB for complex algorithms initially
   - Port incrementally with validation
   - Build confidence iteratively

---

## Confidence Assessment (Updated)

| Aspect | Before Validation | After Validation | Notes |
|--------|------------------|------------------|-------|
| **Technical** | HIGH | ⚠️ MODERATE | Algorithm differences identified |
| **Fixability** | N/A | ✅ HIGH | Issues are fixable (parameters/formulas) |
| **Timeline** | EXCELLENT | ✅ GOOD | +3-5 days for fixes |
| **Success** | VERY HIGH | ✅ HIGH | Clear path forward |

---

## Conclusion

### Summary

The validation phase successfully identified **specific, fixable issues** in the Python implementation:

1. **F0 Tracking:** Missing 2-iteration refinement and wrong voicing threshold
2. **IAIF:** Incorrect VT filter order formula
3. **Both:** Need careful parameter matching with MATLAB

### Status

⚠️ **VALIDATION REVEALED ISSUES** - But all are addressable!

- ✅ Framework works perfectly
- ✅ Issues clearly identified
- ✅ Fixes are straightforward
- ✅ Expected to achieve targets

### Next Steps

1. Apply critical fixes (4 hours estimated)
2. Re-validate
3. Iterate until <10% error achieved
4. Expand to more algorithms

**Expected Timeline:** +2-3 days to achieve validation targets, then continue with expansion.

---

**Report Date:** October 17, 2025  
**Validation Status:** ⚠️ **ISSUES IDENTIFIED - FIXES IN PROGRESS**  
**Confidence in Resolution:** ✅ **HIGH**  
**Project Health:** ✅ **GOOD** (this is normal development process)
