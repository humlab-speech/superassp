// audio_input.hpp — read one channel of an AsspDataObj audio matrix without
// materialising intermediate R copies.
//
// Replaces the pattern
//   NumericMatrix m = audio_list["audio"];          // INTSXP -> full REALSXP copy
//   NumericVector c = m.column(0);                  // column copy
//   std::vector<double> w(c.begin(), c.end());      // third copy
// with a single pass. Values are identical to Rcpp's coercion: int -> double is
// exact and NA_INTEGER maps to NA_REAL, as Rf_coerceVector does.
#pragma once
#include <Rcpp.h>
#include <vector>

namespace sasp {

// Number of rows of the audio matrix (same error as `NumericMatrix m = x` for a
// non-matrix input).
inline int audio_nrow(SEXP audio) {
  if (!Rf_isMatrix(audio)) throw Rcpp::not_a_matrix();
  return INTEGER(Rf_getAttrib(audio, R_DimSymbol))[0];
}

// Call f(i, value) for every sample i of channel `ch` (0-based), where value
// is the double Rcpp's NumericMatrix coercion would yield.
template <typename F>
inline void audio_channel_for_each(SEXP audio, int ch, F f) {
  switch (TYPEOF(audio)) {
  case INTSXP: {
    const int nrow = audio_nrow(audio);
    const int* src = INTEGER(audio) + (R_xlen_t)ch * nrow;
    for (int i = 0; i < nrow; i++)
      f(i, (src[i] == NA_INTEGER) ? NA_REAL : (double)src[i]);
    break;
  }
  case REALSXP: {
    const int nrow = audio_nrow(audio);
    const double* src = REAL(audio) + (R_xlen_t)ch * nrow;
    for (int i = 0; i < nrow; i++) f(i, src[i]);
    break;
  }
  default: {
    // Other storage modes: keep Rcpp's coercion semantics exactly.
    Rcpp::NumericMatrix m(audio);
    const int nrow = m.nrow();
    for (int i = 0; i < nrow; i++) f(i, m(i, ch));
  }
  }
}

// Channel `ch` (0-based) of the audio matrix as doubles.
inline std::vector<double> audio_channel_double(SEXP audio, int ch = 0) {
  if (TYPEOF(audio) == REALSXP) {
    const int nrow = audio_nrow(audio);
    const double* src = REAL(audio) + (R_xlen_t)ch * nrow;
    return std::vector<double>(src, src + nrow);
  }
  std::vector<double> out;
  if (TYPEOF(audio) == INTSXP) out.resize(audio_nrow(audio));
  else out.resize(Rcpp::NumericMatrix(audio).nrow());
  audio_channel_for_each(audio, ch, [&](int i, double v) { out[i] = v; });
  return out;
}

} // namespace sasp
