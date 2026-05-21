// vat_lpc.h — LPC primitives: autocorrelation method (lpcauto), Burg,
// and AR<->RF<->RR conversions. VOICEBOX-compatible sign conventions.

#ifndef VAT_LPC_H
#define VAT_LPC_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "vat_dsp.h"

namespace vat {

// Levinson-Durbin recursion from autocorrelation r[0..p].
// Returns ar (length p+1, ar(0)=1) and residual energy e.
inline void levinson(const arma::vec& r, int p, arma::vec& ar, double& e) {
  ar.zeros(p + 1);
  ar(0) = 1.0;
  if (r(0) == 0.0) { e = 0.0; return; }
  arma::vec a(p, arma::fill::zeros);
  e = r(0);
  for (int m = 0; m < p; ++m) {
    double num = r(m + 1);
    for (int j = 0; j < m; ++j) num += a(j) * r(m - j);
    double k = -num / e;
    arma::vec a_new = a;
    for (int j = 0; j < m; ++j) a_new(j) += k * a(m - 1 - j);
    a_new(m) = k;
    a = a_new;
    e *= (1.0 - k * k);
    if (e <= 0.0) { e = 1e-12; break; }
  }
  for (int j = 0; j < p; ++j) ar(j + 1) = a(j);
}

// VOICEBOX lpcauto: Hamming window per frame, biased autocorr, Levinson.
// For single-frame use (t == length(s)): returns ar (p+1) with ar(0)=1, residual energy e.
inline void lpcauto(const arma::vec& s, int p, arma::vec& ar, double& e) {
  int n = s.n_elem;
  arma::vec w = hamming(n);
  arma::vec wd = s % w;
  arma::vec r(p + 1, arma::fill::zeros);
  for (int k = 0; k <= p; ++k)
    for (int i = 0; i < n - k; ++i) r(k) += wd(i) * wd(i + k);
  levinson(r, p, ar, e);
  // VOICEBOX e = rr * ar' (with windowed signal):
  // recompute to match MATLAB exactly
  double ee = 0.0;
  for (int k = 0; k <= p; ++k) ee += r(k) * ar(k);
  e = ee;
}

// Burg's method on raw signal (no windowing). Length-p AR with ar(0)=1.
inline void burg(const arma::vec& s, int p, arma::vec& ar, double& e) {
  int N = s.n_elem;
  arma::vec ef = s, eb = s;
  arma::vec a(p + 1, arma::fill::zeros); a(0) = 1.0;
  e = arma::dot(s, s) / double(N);
  for (int m = 0; m < p; ++m) {
    double num = 0.0, den = 0.0;
    for (int n = m + 1; n < N; ++n) {
      num += ef(n) * eb(n - 1);
      den += ef(n) * ef(n) + eb(n - 1) * eb(n - 1);
    }
    double k = (den != 0.0) ? -2.0 * num / den : 0.0;
    arma::vec a_new = a;
    for (int j = 1; j <= m + 1; ++j) a_new(j) = a(j) + k * a(m + 1 - j);
    a = a_new;
    arma::vec ef_new = ef, eb_new = eb;
    for (int n = m + 1; n < N; ++n) {
      ef_new(n) = ef(n) + k * eb(n - 1);
      eb_new(n) = eb(n - 1) + k * ef(n);
    }
    ef = ef_new; eb = eb_new;
    e *= (1.0 - k * k);
  }
  ar = a;
}

// AR -> reflection coefficients (Schur recursion, descending).
inline arma::vec lpcar2rf(const arma::vec& ar) {
  int p = ar.n_elem - 1;
  arma::vec a = ar.subvec(1, p) / ar(0);
  arma::vec k(p, arma::fill::zeros);
  for (int m = p - 1; m >= 0; --m) {
    k(m) = a(m);
    if (std::fabs(k(m)) >= 1.0) k(m) = std::copysign(1.0 - 1e-12, k(m));
    arma::vec a_new(m, arma::fill::zeros);
    for (int j = 0; j < m; ++j)
      a_new(j) = (a(j) - k(m) * a(m - 1 - j)) / (1.0 - k(m) * k(m));
    if (m > 0) a = a_new;
  }
  return k;
}

// AR -> autocorrelation (p+1 lags, normalized so r(0)=1/e or arbitrary scale).
// Uses Levinson backward step.
inline arma::vec lpcar2ra(const arma::vec& ar) {
  int p = ar.n_elem - 1;
  arma::vec k = lpcar2rf(ar);
  arma::vec r(p + 1, arma::fill::zeros);
  r(0) = 1.0;
  arma::vec a(1, arma::fill::ones);
  double e = 1.0;
  for (int m = 0; m < p; ++m) {
    double s = 0.0;
    for (int j = 0; j < (int)a.n_elem; ++j) s += a(j) * r(m - j);
    // Should equal -k(m) * e
    r(m + 1) = -k(m) * e - s;
    arma::vec a_new(a.n_elem + 1, arma::fill::zeros);
    for (arma::uword j = 0; j < a.n_elem; ++j) a_new(j) = a(j);
    for (arma::uword j = 0; j < a.n_elem; ++j) a_new(j + 1) += k(m) * a(a.n_elem - 1 - j);
    a = a_new;
    e *= (1.0 - k(m) * k(m));
  }
  return r;
}

// RF -> RR (reflection coefficients to AR-form autocorrelation).
inline arma::vec lpcrf2rr(const arma::vec& k) {
  int p = k.n_elem;
  arma::vec a(p + 1, arma::fill::zeros); a(0) = 1.0;
  double e = 1.0;
  arma::vec r(p + 1, arma::fill::zeros); r(0) = 1.0;
  for (int m = 0; m < p; ++m) {
    double s = 0.0;
    for (int j = 0; j <= m; ++j) s += a(j) * r(m - j);
    r(m + 1) = -k(m) * e - s;
    arma::vec a_new = a;
    for (int j = 1; j <= m + 1; ++j) a_new(j) = a(j) + k(m) * a(m + 1 - j);
    a = a_new;
    e *= (1.0 - k(m) * k(m));
  }
  return r;
}

// Itakura distance (symmetric variant from VOICEBOX distitar).
inline double distitar(const arma::vec& ar1, const arma::vec& ar2, bool symmetric = false) {
  arma::vec r1 = lpcar2ra(ar1);
  arma::vec r2 = lpcar2ra(ar2);
  // d = log( ar1 * R(ar2) * ar1' / (ar2 * R(ar2) * ar2') )
  int p = ar1.n_elem;
  auto quad = [&](const arma::vec& a, const arma::vec& r) {
    double s = 0.0;
    for (int i = 0; i < p; ++i)
      for (int j = 0; j < p; ++j)
        s += a(i) * a(j) * r(std::abs(i - j));
    return s;
  };
  double d1 = std::log(quad(ar1, r2) / quad(ar2, r2));
  if (!symmetric) return d1;
  double d2 = std::log(quad(ar2, r1) / quad(ar1, r1));
  return 0.5 * (d1 + d2);
}

}  // namespace vat
#endif
