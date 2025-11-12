# SIMD Phase 2 Optimization Roadmap

**Date**: 2025-11-12
**Status**: Planning Phase
**Phase 1 Status**: ✅ Complete (YIN + ESTK PDA)

---

## Executive Summary

Phase 1 of SIMD optimization successfully delivered production-ready implementations for YIN and ESTK PDA pitch trackers, achieving 13.3x and 80x real-time performance respectively. Phase 2 focuses on additional ESTK PDA optimizations and TANDEM pitch tracker improvements, with expected cumulative speedups of 2-8x.

### Phase 1 Achievements ✅

1. **YIN difference function** - 13.3x real-time (304ms for 4s audio)
2. **ESTK PDA super-resolution** - 80x real-time (50.5ms for 4s audio)
3. Both algorithms production-ready and fully tested

### Phase 2 Candidates

Three high-value optimization targets identified:

1. **ESTK PDA Peak Scoring** - Expected 4-5x speedup (HIGH priority)
2. **ESTK PDA Refinement** - Expected 4x speedup (MEDIUM-HIGH priority)
3. **TANDEM RMS/Scaling** - Expected 3-8x speedup (MEDIUM priority)

---

## Phase 2 Optimization Targets

### 1. ESTK PDA Peak Scoring Loop

**Priority**: HIGH
**Expected Speedup**: 4-5x
**Difficulty**: EASY
**File**: `src/estk_pda.cpp`

#### Current Implementation

Located in peak scoring function (search for "yy += y_val * y_val"):

```cpp
// Compute cross-correlation for peak scoring
for (size_t i = 0; i < window_size; i++) {
    float y_val = y[i];
    float z_val = z[i];

    yy += y_val * y_val;
    zz += z_val * z_val;
    yz += y_val * z_val;
}
```

#### Why This Is High Value

- **Same pattern as super-resolution loop** (already implemented)
- 3 independent accumulations (yy, zz, yz)
- Contiguous memory access
- No dependencies between iterations
- Textbook SIMD case

#### Proposed SIMD Implementation

```cpp
#ifdef RCPPXSIMD_AVAILABLE
using batch_type = xsimd::simd_type<float>;
constexpr size_t simd_size = batch_type::size;

batch_type yy_vec(0.0f), zz_vec(0.0f), yz_vec(0.0f);
size_t i = 0;

// SIMD loop: process simd_size elements at once
for (; i + simd_size <= window_size; i += simd_size) {
    alignas(32) float buf_y[simd_size];
    alignas(32) float buf_z[simd_size];

    for (size_t j = 0; j < simd_size; j++) {
        buf_y[j] = y[i + j];
        buf_z[j] = z[i + j];
    }

    batch_type y_batch, z_batch;
    y_batch.load_aligned(buf_y);
    z_batch.load_aligned(buf_z);

    yy_vec += y_batch * y_batch;
    zz_vec += z_batch * z_batch;
    yz_vec += y_batch * z_batch;
}

// Horizontal reductions
yy = xsimd::hadd(yy_vec);
zz = xsimd::hadd(zz_vec);
yz = xsimd::hadd(yz_vec);

// Scalar tail loop
for (; i < window_size; i++) {
    float y_val = y[i];
    float z_val = z[i];
    yy += y_val * y_val;
    zz += z_val * z_val;
    yz += y_val * z_val;
}
#else
// Existing scalar implementation
#endif
```

#### Implementation Checklist

- [ ] Locate peak scoring function in `src/estk_pda.cpp`
- [ ] Add `#ifdef RCPPXSIMD_AVAILABLE` guard
- [ ] Implement SIMD loop with aligned buffers
- [ ] Add horizontal reductions with `xsimd::hadd()`
- [ ] Add scalar tail loop for remaining elements
- [ ] Test for correctness (determinism check)
- [ ] Benchmark performance (100 iterations)
- [ ] Update documentation

#### Expected Impact

- **Speedup**: 4-5x (ARM NEON), 6-8x (x86 AVX2)
- **Effort**: Low (copy pattern from super-resolution loop)
- **Risk**: Very low (well-tested pattern)

---

### 2. ESTK PDA Refinement Loop (Dot Product)

**Priority**: MEDIUM-HIGH
**Expected Speedup**: 4x
**Difficulty**: EASY
**File**: `src/estk_pda.cpp`

#### Current Implementation

Located in refinement function (search for "dot += x[i] * y[i]"):

```cpp
// Pure dot product
float dot = 0.0f;
for (size_t i = 0; i < length; i++) {
    dot += x[i] * y[i];
}
```

#### Why This Is High Value

- **Classic SIMD textbook case** (vectorized dot product)
- Single accumulation (simpler than peak scoring)
- Contiguous memory access
- No dependencies
- Well-studied optimization pattern

#### Proposed SIMD Implementation

```cpp
#ifdef RCPPXSIMD_AVAILABLE
using batch_type = xsimd::simd_type<float>;
constexpr size_t simd_size = batch_type::size;

batch_type dot_vec(0.0f);
size_t i = 0;

// SIMD loop
for (; i + simd_size <= length; i += simd_size) {
    alignas(32) float buf_x[simd_size];
    alignas(32) float buf_y[simd_size];

    for (size_t j = 0; j < simd_size; j++) {
        buf_x[j] = x[i + j];
        buf_y[j] = y[i + j];
    }

    batch_type x_batch, y_batch;
    x_batch.load_aligned(buf_x);
    y_batch.load_aligned(buf_y);

    dot_vec += x_batch * y_batch;
}

// Horizontal reduction
float dot = xsimd::hadd(dot_vec);

// Scalar tail loop
for (; i < length; i++) {
    dot += x[i] * y[i];
}
#else
// Existing scalar implementation
#endif
```

#### Implementation Checklist

- [ ] Locate refinement function in `src/estk_pda.cpp`
- [ ] Add `#ifdef RCPPXSIMD_AVAILABLE` guard
- [ ] Implement SIMD dot product
- [ ] Add horizontal reduction with `xsimd::hadd()`
- [ ] Add scalar tail loop
- [ ] Test for correctness (determinism check)
- [ ] Benchmark performance
- [ ] Update documentation

#### Expected Impact

- **Speedup**: ~4x (ARM NEON), ~6-8x (x86 AVX2)
- **Effort**: Very low (simplest SIMD pattern)
- **Risk**: Very low (standard dot product)

---

### 3. TANDEM RMS Calculation

**Priority**: MEDIUM
**Expected Speedup**: 3-5x
**Difficulty**: EASY
**File**: `src/tandem_memory.cpp`

#### Current Implementation

Located in RMS calculation (search for "sum += sample * sample"):

```cpp
// RMS calculation
float sum = 0.0f;
for (int i = 0; i < frame_size; i++) {
    float sample = audio_frame[i];
    sum += sample * sample;
}
float rms = sqrt(sum / frame_size);
```

#### Why This Is Valuable

- **Simple squared sum** with single accumulation
- Appears in every frame (called frequently)
- Contiguous memory access
- No dependencies
- Easy to vectorize

#### Proposed SIMD Implementation

```cpp
#ifdef RCPPXSIMD_AVAILABLE
using batch_type = xsimd::simd_type<float>;
constexpr size_t simd_size = batch_type::size;

batch_type sum_vec(0.0f);
int i = 0;

// SIMD loop
for (; i + static_cast<int>(simd_size) <= frame_size; i += simd_size) {
    alignas(32) float buf[simd_size];

    for (size_t j = 0; j < simd_size; j++) {
        buf[j] = audio_frame[i + j];
    }

    batch_type samples;
    samples.load_aligned(buf);

    sum_vec += samples * samples;
}

// Horizontal reduction
float sum = xsimd::hadd(sum_vec);

// Scalar tail loop
for (; i < frame_size; i++) {
    float sample = audio_frame[i];
    sum += sample * sample;
}

float rms = sqrt(sum / frame_size);
#else
// Existing scalar implementation
#endif
```

#### Implementation Checklist

- [ ] Locate RMS calculation in `src/tandem_memory.cpp`
- [ ] Add `#ifdef RCPPXSIMD_AVAILABLE` guard
- [ ] Implement SIMD squared sum
- [ ] Add horizontal reduction
- [ ] Add scalar tail loop
- [ ] Test for correctness (compare RMS values)
- [ ] Benchmark performance
- [ ] Update documentation

#### Expected Impact

- **Speedup**: 3-5x for RMS calculation
- **Overall impact**: Moderate (RMS is small part of TANDEM)
- **Effort**: Very low (simplest pattern)
- **Risk**: Very low (single accumulation)

---

## Secondary Candidates (Optional)

### 4. TANDEM Scaling Operations

**Priority**: MEDIUM-LOW
**Expected Speedup**: 3-4x
**Difficulty**: EASY

Vectorize array scaling operations:
```cpp
for (int i = 0; i < size; i++) {
    output[i] = input[i] * scale_factor;
}
```

### 5. TANDEM Max Reduction

**Priority**: MEDIUM-LOW
**Expected Speedup**: 4-6x
**Difficulty**: EASY

Vectorize max-finding:
```cpp
float max_val = -INFINITY;
for (int i = 0; i < size; i++) {
    if (input[i] > max_val) max_val = input[i];
}
```

---

## Implementation Strategy

### Recommended Order

1. **ESTK PDA Peak Scoring** (HIGH priority, easy)
   - Reuse super-resolution pattern
   - Expected 4-5x speedup
   - Low risk, high reward

2. **ESTK PDA Refinement** (MEDIUM-HIGH priority, easiest)
   - Classic dot product pattern
   - Expected 4x speedup
   - Lowest risk

3. **TANDEM RMS** (MEDIUM priority, easy)
   - Simple squared sum
   - Expected 3-5x speedup
   - Moderate overall impact

### Testing Protocol

For each optimization:

1. **Correctness**:
   - Run 3 iterations, verify 100% identical outputs
   - Compare with scalar baseline
   - Check edge cases (small arrays, odd sizes)

2. **Performance**:
   - Benchmark with 100 iterations
   - Measure median, mean, min, max
   - Calculate real-time factor
   - Compare to expected speedup

3. **Portability**:
   - Test on ARM NEON (if available)
   - Test on x86 AVX2 (if available)
   - Verify scalar fallback works

---

## Expected Cumulative Impact

### Phase 1 Baseline (Already Achieved)

- **YIN**: 13.3x real-time
- **ESTK PDA super-resolution**: 80x real-time

### Phase 2 Additions (Projected)

If all three optimizations implemented:

- **ESTK PDA peak scoring**: +4-5x contribution
- **ESTK PDA refinement**: +4x contribution
- **TANDEM RMS**: +3-5x contribution

**Cumulative speedup**: 2-4x overall DSP performance improvement

### Phase 3+ Potential

- Additional TANDEM optimizations: +2-3x
- OpenSMILE SIMD: +2-4x (if not already optimized)
- Other DSP functions: +1-3x

**Total potential**: 5-10x overall DSP performance vs original baseline

---

## Technical Considerations

### Memory Alignment

All Phase 2 optimizations use the same pattern:
- Fixed 32-byte alignment (`alignas(32)`)
- Temporary aligned buffers for loading
- Horizontal reduction with `xsimd::hadd()`
- Scalar tail loops for remaining elements

### xsimd v7 API Consistency

All implementations follow Phase 1 patterns:
```cpp
using batch_type = xsimd::simd_type<float>;
batch_type vec(0.0f);
vec.load_aligned(buffer);  // Member function
float result = xsimd::hadd(vec);  // Free function
```

### Compilation Flags

Already configured in `src/Makevars`:
```makefile
PKG_CPPFLAGS = ... -DRCPPXSIMD_AVAILABLE
```

No additional changes required.

---

## Risk Assessment

### Low Risk Optimizations

✅ **ESTK PDA Peak Scoring** - Reuses proven pattern
✅ **ESTK PDA Refinement** - Classic dot product
✅ **TANDEM RMS** - Simple squared sum

All three optimizations are:
- Well-understood SIMD patterns
- No complex data dependencies
- Easy to test and verify
- Backward compatible (scalar fallback)

### Medium Risk (Not Recommended for Phase 2)

⚠️ **ESTK PDA Main Correlation** - Strided access patterns
⚠️ **Complex TANDEM loops** - Multiple dependencies

These require more analysis and testing.

---

## Resource Estimates

### Time Estimates (per optimization)

- **Implementation**: 1-2 hours
- **Testing**: 1 hour
- **Documentation**: 30 minutes
- **Total per optimization**: 2.5-3.5 hours

### Total Phase 2 Estimate

- **3 optimizations**: 7.5-10.5 hours
- **Documentation**: 2 hours
- **Total**: ~10-12 hours

### Recommended Phasing

**Quick wins (4-6 hours)**:
1. ESTK PDA Peak Scoring (2 hours)
2. ESTK PDA Refinement (2 hours)

**Extended Phase 2 (10-12 hours)**:
3. TANDEM RMS (2 hours)
4. Additional testing and documentation (4 hours)

---

## Deliverables

### Code

1. Modified `src/estk_pda.cpp` with peak scoring and refinement SIMD
2. Modified `src/tandem_memory.cpp` with RMS SIMD (optional)
3. Updated test scripts for each optimization

### Documentation

1. Individual completion documents (pattern: `*_SIMD_COMPLETE.md`)
2. Updated `SIMD_INTEGRATION_SUMMARY.md`
3. Updated `SIMD_IMPLEMENTATION_STATUS.md`
4. Final session summary

### Testing

1. Correctness tests (3 iterations each)
2. Performance benchmarks (100 iterations each)
3. Cross-platform verification (if applicable)

---

## Success Criteria

### Phase 2 Complete When:

✅ All optimizations implemented and tested
✅ 100% deterministic outputs verified
✅ Performance benchmarks show expected speedups
✅ Documentation complete
✅ No regressions in existing tests
✅ Code compiles cleanly on all platforms

### Minimum Success:

✅ ESTK PDA peak scoring implemented and tested
✅ At least 4x speedup measured
✅ Documentation updated

---

## Recommendations

### For Immediate Next Steps

1. **Start with ESTK PDA Peak Scoring**
   - Reuse super-resolution pattern
   - Lowest risk, highest impact
   - Quick win to build momentum

2. **Then ESTK PDA Refinement**
   - Simplest pattern (dot product)
   - Completes all ESTK PDA optimizations
   - High confidence of success

3. **TANDEM RMS (Optional)**
   - If time permits
   - Smaller overall impact
   - Still valuable

### For Long-Term Strategy

- **Focus on high-impact, low-risk optimizations**
- **Prioritize user-facing DSP functions**
- **Build library of SIMD patterns for reuse**
- **Document thoroughly for future work**

---

## Conclusion

Phase 2 SIMD optimization offers 2-4x additional speedups with low risk and moderate effort. The three primary targets (ESTK PDA peak scoring, ESTK PDA refinement, TANDEM RMS) all use proven patterns from Phase 1 and can be implemented incrementally with clear testing and success criteria.

**Recommendation**: Proceed with Phase 2, starting with ESTK PDA peak scoring as the highest-value, lowest-risk optimization.

---

**Created**: 2025-11-12
**Status**: Planning Complete
**Ready**: Yes - all patterns proven in Phase 1
