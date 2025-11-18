# SIMD Optimization Plan for superassp

**Date**: 2025-11-10
**Analysis Tool**: Gemini CLI (gemini-2.0-flash-exp)
**Target Library**: RcppXsimd (https://github.com/OHDSI/RcppXsimd)
**Status**: Assessment Complete - Ready for Implementation

---

## Executive Summary

Analysis of the superassp C++ codebase has identified **8 high-value SIMD vectorization opportunities** across 3 core DSP implementation files. The top 5 candidates offer **2-8x performance improvements** through RcppXsimd vectorization.

**Impact**: Optimizing the identified hotspots could reduce pitch tracking computation time by 50-75%, benefiting all downstream DSP functions.

---

## Top Priority Candidates

### 1. YIN Difference Function ⭐ HIGHEST IMPACT
**File**: `src/yin_wrapper.cpp`
**Function**: YIN autocorrelation difference function
**Estimated Speedup**: **4-8x**
**Priority**: **CRITICAL**

#### Current Implementation (Non-Vectorized)
```cpp
// Compute difference function (autocorrelation-based)
for (int tau = 1; tau < halfBufferSize; tau++) {
    yinBuffer[tau] = 0.0f;
    for (int i = 0; i < halfBufferSize; i++) {
        float delta = buffer[i] - buffer[i + tau];
        yinBuffer[tau] += delta * delta;
    }
}
```

#### Optimization Strategy
**Algorithm**: Vectorized sum of squared differences
- Inner loop: `(buffer[i] - buffer[i + tau])^2` is a pure vector operation
- Independent across `i` - perfect for SIMD
- Reduction (sum) can be vectorized with horizontal add

#### Proposed RcppXsimd Implementation
```cpp
#include <RcppXsimd.h>

// Compute difference function with SIMD
for (int tau = 1; tau < halfBufferSize; tau++) {
    using batch = xsimd::batch<float>;
    constexpr size_t simd_size = batch::size;

    float sum = 0.0f;
    int i = 0;

    // SIMD loop (process simd_size elements at once)
    for (; i + simd_size <= halfBufferSize; i += simd_size) {
        batch b1 = batch::load_unaligned(&buffer[i]);
        batch b2 = batch::load_unaligned(&buffer[i + tau]);
        batch delta = b1 - b2;
        batch sq = delta * delta;
        sum += xsimd::reduce_add(sq);
    }

    // Scalar tail loop (remaining elements)
    for (; i < halfBufferSize; i++) {
        float delta = buffer[i] - buffer[i + tau];
        sum += delta * delta;
    }

    yinBuffer[tau] = sum;
}
```

#### Performance Estimate
- **Baseline**: 200-300ms for 3s audio (1.5M operations)
- **Optimized**: 25-75ms (4-8x speedup depending on CPU)
- **Benefit**: All functions using YIN pitch tracking accelerate

---

### 2. ESTK PDA Super-Resolution Refinement
**File**: `src/estk_pda.cpp`
**Function**: Super-resolution pitch period refinement
**Estimated Speedup**: **4-6x**
**Priority**: **HIGH**

#### Current Implementation
```cpp
// Super-resolution refinement (6 independent accumulations)
for (size_t N = 0; N < fft_size; N++) {
    float y_N = y[N];
    float z_N = z[N];

    xx_N += y_N * y_N;
    yy_N += z_N * z_N;
    xy_N += y_N * z_N;

    // ... more accumulations
}
```

#### Optimization Strategy
**Algorithm**: Vectorized multiply-accumulate with horizontal reduction
- 6 independent accumulations can be computed in parallel
- Pure vector operations (multiply, add)
- Final reduction with horizontal sum

#### Proposed RcppXsimd Implementation
```cpp
using batch = xsimd::batch<float>;
constexpr size_t simd_size = batch::size;

batch xx_vec = batch(0.0f), yy_vec = batch(0.0f), xy_vec = batch(0.0f);
batch xz_vec = batch(0.0f), yz_vec = batch(0.0f), zz_vec = batch(0.0f);

size_t N = 0;
for (; N + simd_size <= fft_size; N += simd_size) {
    batch y_N = batch::load_unaligned(&y[N]);
    batch z_N = batch::load_unaligned(&z[N]);

    xx_vec += y_N * y_N;
    yy_vec += z_N * z_N;
    xy_vec += y_N * z_N;
    // ... other accumulations
}

// Horizontal reduction
float xx_N = xsimd::reduce_add(xx_vec);
float yy_N = xsimd::reduce_add(yy_vec);
float xy_N = xsimd::reduce_add(xy_vec);
// ... other reductions

// Scalar tail
for (; N < fft_size; N++) {
    float y_N = y[N];
    float z_N = z[N];
    xx_N += y_N * y_N;
    yy_N += z_N * z_N;
    xy_N += y_N * z_N;
}
```

#### Performance Estimate
- **Baseline**: ~50ms per refinement iteration
- **Optimized**: ~8-12ms (4-6x speedup)

---

### 3. ESTK PDA Peak Scoring Loop
**File**: `src/estk_pda.cpp`
**Function**: Peak scoring via correlation
**Estimated Speedup**: **4-5x**
**Priority**: **HIGH**

#### Current Implementation
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

#### Optimization Strategy
Same pattern as super-resolution - 3 independent accumulations, vectorizable with horizontal reduction.

#### Proposed Implementation
```cpp
using batch = xsimd::batch<float>;
batch yy_vec = batch(0.0f), zz_vec = batch(0.0f), yz_vec = batch(0.0f);

size_t i = 0;
for (; i + batch::size <= window_size; i += batch::size) {
    batch y_val = batch::load_unaligned(&y[i]);
    batch z_val = batch::load_unaligned(&z[i]);

    yy_vec += y_val * y_val;
    zz_vec += z_val * z_val;
    yz_vec += y_val * z_val;
}

float yy = xsimd::reduce_add(yy_vec);
float zz = xsimd::reduce_add(zz_vec);
float yz = xsimd::reduce_add(yz_vec);

// Scalar tail
for (; i < window_size; i++) {
    float y_val = y[i];
    float z_val = z[i];
    yy += y_val * y_val;
    zz += z_val * z_val;
    yz += y_val * z_val;
}
```

---

### 4. ESTK PDA Refinement Loop (Dot Product)
**File**: `src/estk_pda.cpp`
**Function**: Dot product for refinement
**Estimated Speedup**: **~4x**
**Priority**: **MEDIUM-HIGH**

#### Current Implementation
```cpp
// Pure dot product
float dot = 0.0f;
for (size_t i = 0; i < length; i++) {
    dot += x[i] * y[i];
}
```

#### Optimization Strategy
**Algorithm**: Vectorized dot product with horizontal reduction
- Textbook SIMD case
- Contiguous memory, no dependencies

#### Proposed Implementation
```cpp
using batch = xsimd::batch<float>;
batch dot_vec = batch(0.0f);

size_t i = 0;
for (; i + batch::size <= length; i += batch::size) {
    batch x_val = batch::load_unaligned(&x[i]);
    batch y_val = batch::load_unaligned(&y[i]);
    dot_vec += x_val * y_val;
}

float dot = xsimd::reduce_add(dot_vec);

// Scalar tail
for (; i < length; i++) {
    dot += x[i] * y[i];
}
```

---

### 5. ESTK PDA Main Correlation Loop
**File**: `src/estk_pda.cpp`
**Function**: Main autocorrelation computation
**Estimated Speedup**: **2-3x**
**Priority**: **MEDIUM**

#### Current Implementation
```cpp
// Strided autocorrelation
for (int lag = 0; lag < max_lag; lag++) {
    float sum = 0.0f;
    for (int i = 0; i < window_size - lag; i++) {
        sum += signal[i] * signal[i + lag];
    }
    correlation[lag] = sum;
}
```

#### Optimization Strategy
**Algorithm**: Vectorized strided dot product
- More challenging due to stride (`i + lag`)
- Still benefits from SIMD with gather/scatter or pointer arithmetic
- Lower speedup (2-3x) due to memory access pattern

#### Proposed Implementation
```cpp
for (int lag = 0; lag < max_lag; lag++) {
    using batch = xsimd::batch<float>;
    batch sum_vec = batch(0.0f);

    int i = 0;
    int length = window_size - lag;

    for (; i + batch::size <= length; i += batch::size) {
        batch s1 = batch::load_unaligned(&signal[i]);
        batch s2 = batch::load_unaligned(&signal[i + lag]);
        sum_vec += s1 * s2;
    }

    float sum = xsimd::reduce_add(sum_vec);

    for (; i < length; i++) {
        sum += signal[i] * signal[i + lag];
    }

    correlation[lag] = sum;
}
```

---

## Secondary Candidates

### 6. TANDEM RMS Calculation
**File**: `src/tandem_memory.cpp`
**Function**: Root Mean Square energy computation
**Estimated Speedup**: **3-4x**

```cpp
// Current
double sumE = 0.0;
for (int i = 0; i < numSamples; i++) {
    sumE += samples[i] * samples[i];
}

// Vectorized
using batch = xsimd::batch<double>;
batch sum_vec = batch(0.0);
size_t i = 0;
for (; i + batch::size <= numSamples; i += batch::size) {
    batch s = batch::load_unaligned(&samples[i]);
    sum_vec += s * s;
}
double sumE = xsimd::reduce_add(sum_vec);
for (; i < numSamples; i++) {
    sumE += samples[i] * samples[i];
}
```

---

### 7. TANDEM Scaling Loop
**File**: `src/tandem_memory.cpp`
**Function**: Vector scaling (y = scale * x)
**Estimated Speedup**: **3-5x**

```cpp
// Current
for (int i = 0; i < numSamples; i++) {
    signal[i] = scale * samples[i];
}

// Vectorized
using batch = xsimd::batch<double>;
batch scale_vec = batch(scale);
size_t i = 0;
for (; i + batch::size <= numSamples; i += batch::size) {
    batch s = batch::load_unaligned(&samples[i]);
    batch result = scale_vec * s;
    result.store_unaligned(&signal[i]);
}
for (; i < numSamples; i++) {
    signal[i] = scale * samples[i];
}
```

---

### 8. TANDEM Max Reduction
**File**: `src/tandem_memory.cpp`
**Function**: Find maximum value in array
**Estimated Speedup**: **2-3x**

```cpp
// Current
float max_val = -INFINITY;
for (int i = 0; i < length; i++) {
    if (array[i] > max_val) max_val = array[i];
}

// Vectorized
using batch = xsimd::batch<float>;
batch max_vec = batch(-INFINITY);
size_t i = 0;
for (; i + batch::size <= length; i += batch::size) {
    batch vals = batch::load_unaligned(&array[i]);
    max_vec = xsimd::max(max_vec, vals);
}
float max_val = xsimd::reduce_max(max_vec);
for (; i < length; i++) {
    if (array[i] > max_val) max_val = array[i];
}
```

---

## Implementation Roadmap

### Phase 1: Critical Path Optimization (Week 1-2)
**Goal**: 4-8x speedup on most-used functions

1. **YIN difference function** (`yin_wrapper.cpp`)
   - Effort: 2-3 hours
   - Impact: All YIN-based pitch tracking functions
   - Testing: Compare against reference implementation

2. **ESTK PDA super-resolution** (`estk_pda.cpp`)
   - Effort: 3-4 hours
   - Impact: ESTK pitch detection accuracy and speed
   - Testing: Verify pitch detection accuracy unchanged

### Phase 2: Secondary Optimizations (Week 3)
**Goal**: Complete ESTK PDA vectorization

3. **Peak scoring loop** (`estk_pda.cpp`)
   - Effort: 2 hours

4. **Refinement loop** (`estk_pda.cpp`)
   - Effort: 1-2 hours

5. **Main correlation loop** (`estk_pda.cpp`)
   - Effort: 2-3 hours

### Phase 3: TANDEM Optimizations (Week 4)
**Goal**: Improve auxiliary DSP operations

6. **RMS calculation** (`tandem_memory.cpp`)
   - Effort: 1 hour

7. **Scaling loop** (`tandem_memory.cpp`)
   - Effort: 30 minutes

8. **Max reduction** (`tandem_memory.cpp`)
   - Effort: 30 minutes

---

## Testing Strategy

### 1. Correctness Verification
```r
# Test that SIMD and scalar versions produce identical results
test_that("SIMD YIN matches reference implementation", {
  audio <- read_avaudio("test.wav")

  # Force scalar version
  result_scalar <- trk_yin(audio, use_simd = FALSE, toFile = FALSE)

  # Use SIMD version
  result_simd <- trk_yin(audio, use_simd = TRUE, toFile = FALSE)

  # Results should be identical (within floating-point tolerance)
  expect_equal(result_simd$`F0[Hz]`, result_scalar$`F0[Hz]`, tolerance = 1e-5)
})
```

### 2. Performance Benchmarking
```r
# Benchmark SIMD vs scalar
library(microbenchmark)

audio <- read_avaudio("speech_3s.wav")

results <- microbenchmark(
  scalar = trk_yin(audio, use_simd = FALSE, toFile = FALSE),
  simd = trk_yin(audio, use_simd = TRUE, toFile = FALSE),
  times = 100
)

print(results)
# Expected: SIMD 4-8x faster
```

### 3. Cross-Platform Testing
- **x86_64**: SSE4.2, AVX2, AVX-512 (automatic via xsimd)
- **ARM**: NEON (automatic via xsimd)
- **Fallback**: Scalar path if SIMD unavailable

---

## RcppXsimd Integration

### Package Setup
**Already installed**: RcppXsimd package available

### Build Configuration
Add to `src/Makevars`:
```makefile
PKG_CPPFLAGS += $(shell "${R_HOME}/bin/Rscript" -e "RcppXsimd:::CxxFlags()")
PKG_LIBS += $(shell "${R_HOME}/bin/Rscript" -e "RcppXsimd:::LdFlags()")
```

### Header Inclusion
In each optimized `.cpp` file:
```cpp
#include <RcppXsimd.h>

// Optional: Enable specific SIMD instruction sets
// #define XSIMD_WITH_AVX2
// #define XSIMD_WITH_AVX512
```

### Runtime Detection
```cpp
// Check available SIMD instruction sets
bool has_avx2 = xsimd::avx2::available();
bool has_neon = xsimd::neon::available();

// Use best available implementation
if (has_avx2) {
    // Use AVX2 batch size
} else if (has_neon) {
    // Use NEON batch size
} else {
    // Fallback to scalar
}
```

---

## Expected Performance Impact

### Overall Package Performance
| Scenario | Before | After | Speedup |
|----------|--------|-------|---------|
| **YIN pitch tracking** | 300ms | 40-75ms | **4-8x** |
| **ESTK PDA pitch detection** | 150ms | 30-50ms | **3-5x** |
| **TANDEM processing** | 100ms | 30-40ms | **2.5-3x** |
| **Batch processing (100 files)** | 45s | 8-15s | **3-6x** |

### Function-Level Impact
Functions that will see speedups:
- `trk_yin()`: 4-8x faster
- `trk_estk_pda()`: 3-5x faster
- `trk_tandem()`: 2.5-3x faster
- All functions using these as dependencies

### User-Facing Impact
- **Real-time processing**: More files processed per second
- **Interactive workflows**: Faster feedback during exploration
- **Large-scale analysis**: Batch processing 3-6x faster

---

## Risks and Mitigations

### Risk 1: Numerical Precision
**Issue**: SIMD operations may have slightly different rounding behavior
**Mitigation**:
- Test with `tolerance = 1e-5` (floating-point epsilon)
- Verify DSP output quality unchanged (pitch detection accuracy, etc.)
- Document any differences

### Risk 2: Platform Compatibility
**Issue**: Different SIMD instruction sets across platforms
**Mitigation**:
- xsimd handles platform abstraction automatically
- Always provide scalar fallback path
- Test on x86_64 (SSE/AVX), ARM (NEON), and fallback

### Risk 3: Maintenance Complexity
**Issue**: SIMD code is harder to read and maintain
**Mitigation**:
- Add extensive comments explaining vectorization
- Keep scalar reference implementation
- Provide clear `use_simd` parameter for debugging

### Risk 4: Build Complexity
**Issue**: RcppXsimd adds build dependency
**Mitigation**:
- RcppXsimd is header-only (no linking required)
- Gracefully degrade if not available
- Document build requirements in README

---

## Implementation Guidelines

### Code Style
```cpp
// Always provide scalar fallback
if (use_simd && xsimd::default_arch::available()) {
    // SIMD implementation
    // ...
} else {
    // Original scalar implementation
    // ...
}
```

### Documentation
```cpp
/**
 * Compute YIN difference function with optional SIMD acceleration.
 *
 * @param buffer Input audio buffer
 * @param yinBuffer Output difference function
 * @param halfBufferSize Buffer size / 2
 * @param use_simd Enable SIMD vectorization (default: TRUE)
 *
 * Performance:
 * - Scalar: ~300ms for 3s audio
 * - SIMD (AVX2): ~40-75ms (4-8x speedup)
 */
void compute_yin_difference(float* buffer, float* yinBuffer,
                           int halfBufferSize, bool use_simd = true);
```

### Testing Pattern
```cpp
// Unit test for each SIMD function
TEST(YinTest, SIMDMatchesScalar) {
    std::vector<float> buffer(2048);
    std::vector<float> yin_scalar(1024);
    std::vector<float> yin_simd(1024);

    // Fill buffer with test data
    for (int i = 0; i < buffer.size(); i++) {
        buffer[i] = std::sin(2.0 * M_PI * 440.0 * i / 48000.0);
    }

    // Compute both versions
    compute_yin_difference(buffer.data(), yin_scalar.data(), 1024, false);
    compute_yin_difference(buffer.data(), yin_simd.data(), 1024, true);

    // Compare results (allow 1e-5 tolerance for floating-point)
    for (int i = 0; i < 1024; i++) {
        EXPECT_NEAR(yin_scalar[i], yin_simd[i], 1e-5);
    }
}
```

---

## Success Metrics

### Quantitative
- ✅ 4-8x speedup on YIN pitch tracking
- ✅ 3-5x speedup on ESTK PDA
- ✅ All tests pass with `tolerance = 1e-5`
- ✅ Builds successfully on macOS, Linux, Windows
- ✅ Builds successfully on x86_64 and ARM

### Qualitative
- ✅ Code readability maintained (comments, clear structure)
- ✅ Backward compatibility preserved (all existing tests pass)
- ✅ No regression in DSP quality (pitch detection accuracy unchanged)
- ✅ Documentation updated with performance notes

---

## Next Steps

1. **Immediate**: Start with YIN difference function (highest impact)
2. **Week 1**: Implement and test YIN + ESTK super-resolution
3. **Week 2**: Complete ESTK PDA optimizations
4. **Week 3**: TANDEM optimizations
5. **Week 4**: Performance benchmarking, documentation, release

---

## References

- **RcppXsimd**: https://github.com/OHDSI/RcppXsimd
- **xsimd library**: https://github.com/xtensor-stack/xsimd
- **YIN algorithm**: de Cheveigné & Kawahara (2002)
- **ESTK PDA**: Edinburgh Speech Tools pitch detection
- **Gemini CLI analysis**: Full codebase assessment (2025-11-10)

---

**Document Status**: Complete - Ready for Implementation
**Created**: 2025-11-10
**Last Updated**: 2025-11-10
