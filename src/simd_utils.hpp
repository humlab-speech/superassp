// simd_utils.hpp — Reusable double-precision SIMD primitives for DSP kernels.
//
// Design goal ordering: faithfulness of DSP output first, efficiency second.
// Every helper is written in DOUBLE precision (xsimd::batch<double>) so the
// vectorized path matches the scalar path to within double-precision rounding
// (differences arise only from summation order in the horizontal reduction and
// are < ~1e-12 relative; well inside testthat's default expect_equal tolerance).
// Contrast src/yin_wrapper.cpp, which vectorizes in float because its scalar
// reference was already float.
//
// Each helper follows the established pattern from src/yin_wrapper.cpp:
//   #ifdef RCPPXSIMD_AVAILABLE  -> SIMD body + scalar tail
//   #else                       -> plain scalar fallback
// No architecture flags (-march=native etc.) are required or wanted: xsimd uses
// the baseline ISA the compiler already targets, keeping builds reproducible and
// portable. Header-only; include from any translation unit.
#pragma once

#include <cstddef>
#include <vector>

#ifdef RCPPXSIMD_AVAILABLE
#include <xsimd/xsimd.hpp>
#endif

namespace sasp {

// Dot product: sum_{i=0}^{n-1} a[i] * b[i].
inline double simd_dot(const double* a, const double* b, int n) {
#ifdef RCPPXSIMD_AVAILABLE
  using batch_type = xsimd::simd_type<double>;
  constexpr std::size_t simd_size = batch_type::size;

  batch_type acc(0.0);
  int i = 0;
  for (; i + static_cast<int>(simd_size) <= n; i += static_cast<int>(simd_size)) {
    batch_type va = xsimd::load_unaligned(a + i);
    batch_type vb = xsimd::load_unaligned(b + i);
    acc += va * vb;
  }
  double sum = xsimd::hadd(acc);
  for (; i < n; i++) sum += a[i] * b[i];  // scalar tail
  return sum;
#else
  double sum = 0.0;
  for (int i = 0; i < n; i++) sum += a[i] * b[i];
  return sum;
#endif
}

// Signal energy: sum_{i=0}^{n-1} x[i]^2.
inline double simd_energy(const double* x, int n) {
#ifdef RCPPXSIMD_AVAILABLE
  using batch_type = xsimd::simd_type<double>;
  constexpr std::size_t simd_size = batch_type::size;

  batch_type acc(0.0);
  int i = 0;
  for (; i + static_cast<int>(simd_size) <= n; i += static_cast<int>(simd_size)) {
    batch_type v = xsimd::load_unaligned(x + i);
    acc += v * v;
  }
  double sum = xsimd::hadd(acc);
  for (; i < n; i++) sum += x[i] * x[i];  // scalar tail
  return sum;
#else
  double sum = 0.0;
  for (int i = 0; i < n; i++) sum += x[i] * x[i];
  return sum;
#endif
}

// FIR filter: y[n] = sum_{k=0}^{min(M-1,n)} b[k] * x[n-k], for n in [0, N).
// Equivalent to R's filter(b, 1, x) / stats::filter(..., sides = 1).
// The full-support region (n >= M-1) is a reversed-b dot product handled by
// simd_dot; the ramp-up boundary (n < M-1) stays scalar for exactness.
inline void simd_fir(const double* x, const double* b, double* y, int N, int M) {
  if (M <= 0) return;

  // Ramp-up: partial support, taps k = 0..n.
  int boundary = M - 1 < N ? M - 1 : N;
  for (int n = 0; n < boundary; n++) {
    double sum = 0.0;
    for (int k = 0; k <= n; k++) sum += b[k] * x[n - k];
    y[n] = sum;
  }

  // Full support: y[n] = sum_k b[k]*x[n-k] = sum_j b_rev[j]*x[n-M+1+j].
  // Reverse b once so both operands are contiguous & ascending for simd_dot.
  static thread_local std::vector<double> b_rev;
  b_rev.resize(M);
  for (int k = 0; k < M; k++) b_rev[k] = b[M - 1 - k];

  for (int n = M - 1; n < N; n++) {
    y[n] = simd_dot(b_rev.data(), x + (n - M + 1), M);
  }
}

// Autocorrelation: r[k] = sum_{i=0}^{n-k-1} x[i]*x[i+k], for k = 0..order.
// Writes order+1 lag values into r (caller-allocated buffer, size >= order+1).
// Implemented as one simd_dot() call per lag, so it inherits simd_dot's
// faithfulness contract (matches the scalar double loop to within
// double-precision summation-order rounding).
inline void simd_autocorr(const double* x, int n, int order, double* r) {
  for (int k = 0; k <= order; k++) {
    int len = n - k;
    r[k] = (len > 0) ? simd_dot(x, x + k, len) : 0.0;
  }
}

}  // namespace sasp
