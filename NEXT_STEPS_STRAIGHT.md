# STRAIGHT Integration - Next Steps
**Priority Order for Completion**

---

## CRITICAL: Fix Segfault (MUST DO FIRST)

**Estimated Time**: 2-4 hours  
**Priority**: P0 (Blocker)

### Steps:
1. Create minimal test case
   ```r
   # Test with smallest audio
   test_file <- system.file("samples/sustained/a1.wav", package = "superassp")
   audio <- av::read_audio_bin(test_file, channels = 1)
   audio_short <- audio[1:min(length(audio), 16000)]  # 0.7 seconds at 22050 Hz
   # Save and test
   ```

2. Test with NumPy 1.x
   ```bash
   conda create -n test_numpy1 python=3.11 numpy=1.26.4 scipy=1.11.0
   # Test in this environment
   ```

3. Add extensive logging
   ```python
   # Add prints after every major operation in f0_extraction.py
   print(f"AC candidates shape: {ac_cands.shape}")
   print(f"AC candidates memory: {ac_cands.nbytes} bytes")
   ```

4. Run under debugger
   ```bash
   R -d lldb
   # Set breakpoints in reticulate Python calls
   ```

5. Test Python-only (bypass R)
   ```python
   # Direct Python test without R wrapper
   python3 -c "
   import sys
   sys.path.insert(0, 'inst/python')
   from legacy_STRAIGHT.f0_extraction import MulticueF0v14
   import soundfile as sf
   audio, fs = sf.read('inst/samples/sustained/a1.wav')
   f0, vuv, aux, prm = MulticueF0v14(audio, fs)
   print('SUCCESS')
   "
   ```

**Success Criteria**: Tests run without segfault

---

## HIGH: Verify Baseline (DO AFTER SEGFAULT FIX)

**Estimated Time**: 1 hour  
**Priority**: P1

### Steps:
1. Run full test suite
   ```r
   devtools::test()
   ```

2. Run validation script (create if needed)
   ```r
   # Compare with MATLAB reference
   source("tests/validate_accuracy.R")
   ```

3. Verify 91.3% accuracy claim
4. Document exact test conditions
5. Create baseline reference output

**Success Criteria**: All tests pass, accuracy confirmed

---

## DECISION: Choose Path Forward

**Estimated Time**: 30 minutes  
**Priority**: P1

### Option A: Systematic Improvement to 95%
**Pros**:
- Pure STRAIGHT implementation
- Research-grade accuracy
- Well-documented approach

**Cons**:
- 2-3 weeks effort
- Requires MATLAB for validation
- No guarantee of exact 95%

**Choose if**: Research accuracy critical, time available

### Option B: Hybrid Approach (RECOMMENDED)
**Pros**:
- Immediate >98% accuracy
- Best spectral/synthesis quality
- Low risk, low effort
- All components working

**Cons**:
- Not pure STRAIGHT
- Two F0 algorithms to maintain

**Choose if**: Production use, time limited

### Option C: Accept 91.3%
**Pros**:
- No additional work
- Production-ready now
- Good enough for most use cases

**Cons**:
- Octave errors in low F0
- Below 95% target

**Choose if**: 91.3% meets requirements

---

## IF OPTION A: Systematic Improvement

### Phase 1: Test Infrastructure (2-3 days)
1. Create unit tests per stage
2. Set up MATLAB comparison
3. Multi-speaker test dataset
4. Regression test framework

### Phase 2: Fix Backward Tracking (3-5 days)
1. Implement frame-dependent relaxation
2. Test on 10 speakers
3. Measure accuracy improvement
4. Ensure no regressions

### Phase 3: Improve AC/IF Weighting (2-3 days)
1. Analyze weighting in different F0 ranges
2. Implement frequency-dependent weights
3. A/B test configurations
4. Measure octave error reduction

### Phase 4: Validation (2-3 days)
1. Multi-speaker validation
2. Performance optimization
3. Documentation updates

**Total**: 2-3 weeks to 95%

---

## IF OPTION B: Hybrid Approach

### Implementation (1-2 hours)
1. Create wrapper function
   ```r
   straight_hybrid <- function(audio_file, ...) {
     # F0 with RAPT
     f0_data <- trk_rapt(audio_file, toFile = FALSE)
     
     # Spectral with STRAIGHT
     spec_data <- trk_straight_spec(audio_file, toFile = FALSE)
     
     # Return both
     list(f0 = f0_data, spec = spec_data)
   }
   ```

2. Document hybrid workflow
3. Add examples to documentation
4. Update CLAUDE.md with recommendation

**Total**: 1-2 hours

---

## IF OPTION C: Deploy at 91.3%

### Finalization (2-3 hours)
1. Update documentation with known limitations
2. Add usage recommendations
3. Update NEWS.md
4. Tag release v0.9.0

**Total**: 2-3 hours

---

## Code Quality (OPTIONAL)

**Estimated Time**: 1-2 days  
**Priority**: P2 (After deployment)

1. Add type hints to Python code
2. Improve error messages
3. Add input validation
4. Performance profiling
5. Memory optimization

---

## Timeline Summary

### Minimum (Option B or C)
- **Day 1**: Fix segfault (2-4 hours)
- **Day 1**: Verify baseline (1 hour)
- **Day 1**: Implement hybrid OR finalize (1-3 hours)
- **Total**: 4-8 hours = **1 day**

### Maximum (Option A)
- **Day 1**: Fix segfault + verify (3-5 hours)
- **Week 1**: Test infrastructure (2-3 days)
- **Week 2**: Backward tracking fixes (3-5 days)
- **Week 3**: AC/IF weighting + validation (4-6 days)
- **Total**: **2-3 weeks**

---

## Resource Requirements

### For Any Option
- R development environment ✅
- Python 3.11+ with NumPy/SciPy ✅
- Test audio files ✅
- Debugging tools (lldb/gdb)

### For Option A (95% Improvement)
- MATLAB license (for validation)
- Multi-speaker test datasets (TIMIT, Arctic)
- 2-3 weeks dedicated time
- Code review capacity

---

## Success Metrics

### Phase 1 (Segfault Fix)
- [ ] Tests run without crashing
- [ ] All test suite passes
- [ ] Baseline accuracy verified

### Phase 2 (Deployment)
- [ ] Documentation complete
- [ ] Examples working
- [ ] NEWS.md updated
- [ ] Release tagged

### Phase 3 (If Option A)
- [ ] 95% accuracy achieved
- [ ] Multi-speaker validated
- [ ] Performance optimized
- [ ] Comprehensive docs

---

## Risk Mitigation

### Segfault Not Fixed Quickly
**Fallback**: Use Python-only version, create simpler R wrapper

### Can't Reach 95% in Timeframe
**Fallback**: Document progress, release at achieved accuracy

### MATLAB Not Available
**Fallback**: Use cross-validation with other algorithms (RAPT, SWIPE, DIO)

---

## Communication

### Status Updates
- Daily during segfault debugging
- Weekly during systematic improvement
- On completion of each phase

### Documentation
- Update STRAIGHT_INTEGRATION_STATUS_NOV1_2025.md
- Keep NEXT_STEPS_STRAIGHT.md current
- Log findings in session notes

---

## Contact Points

**Package Maintainer**: superassp team  
**Repository**: https://github.com/humlab-speech/superassp  
**Issues**: https://github.com/humlab-speech/superassp/issues  
**Documentation**: See CLAUDE.md, README.md  

---

*Next Steps Document - Last Updated: November 1, 2025*
