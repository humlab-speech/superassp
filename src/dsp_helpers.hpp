// dsp_helpers.hpp — Shared DSP helpers for iaif.cpp and gfmiaif.cpp
#pragma once
#include <vector>
#include <algorithm>

// Levinson-Durbin recursion: returns LPC coefficients [1, a1, ..., a_order]
inline std::vector<double> levinson_durbin(const std::vector<double>& r, int order) {
  std::vector<double> a(order + 1, 0.0);
  a[0] = 1.0;
  if (r[0] <= 0.0) return a;

  double e = r[0];
  std::vector<double> a_prev(order + 1, 0.0);

  for (int m = 1; m <= order; m++) {
    double k = -r[m];
    for (int j = 1; j < m; j++) k -= a[j] * r[m - j];
    k /= e;
    std::copy(a.begin(), a.begin() + m, a_prev.begin());
    a[m] = k;
    for (int j = 1; j < m; j++) a[j] = a_prev[j] + k * a_prev[m - j];
    e *= (1.0 - k * k);
    if (e <= 0.0) break;
  }
  return a;
}

// FIR filter: y = filter(b, 1, x)  (FIR numerator b, denominator = 1)
inline std::vector<double> fir_filter(const std::vector<double>& b,
                                       const double* x, int N) {
  int M = (int)b.size();
  std::vector<double> y(N, 0.0);
  for (int n = 0; n < N; n++) {
    double sum = 0.0;
    for (int k = 0; k < M && k <= n; k++) sum += b[k] * x[n - k];
    y[n] = sum;
  }
  return y;
}

// IIR first-order leaky integrator: y[n] = x[n] + d * y[n-1]
inline std::vector<double> leaky_integrate(const double* x, int N, double d) {
  std::vector<double> y(N, 0.0);
  y[0] = x[0];
  for (int n = 1; n < N; n++) y[n] = x[n] + d * y[n - 1];
  return y;
}
