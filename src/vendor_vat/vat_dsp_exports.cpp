// vat_dsp_exports.cpp — Rcpp-exposed wrappers around vat_dsp.h + vat_lpc.h
// primitives, used by R-level unit tests against MATLAB reference outputs.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "vat_dsp.h"
#include "vat_lpc.h"

using namespace Rcpp;

// Convert arma::vec to flat NumericVector (no dim attribute).
static inline NumericVector vec2nv(const arma::vec& v) {
  NumericVector out(v.n_elem);
  for (arma::uword i = 0; i < v.n_elem; ++i) out[i] = v(i);
  return out;
}

// [[Rcpp::export]]
NumericVector vat_hamming_cpp(int N, bool symmetric = true) {
  return vec2nv(vat::hamming(N, symmetric));
}

// [[Rcpp::export]]
NumericVector vat_hanning_cpp(int N, bool symmetric = true) {
  return vec2nv(vat::hanning(N, symmetric));
}

// [[Rcpp::export]]
NumericVector vat_kaiser_cpp(int N, double beta) {
  return vec2nv(vat::kaiser(N, beta));
}

// [[Rcpp::export]]
NumericVector vat_filter_cpp(NumericVector b, NumericVector a, NumericVector x,
                              Nullable<NumericVector> zi = R_NilValue) {
  arma::vec zi_v;
  if (zi.isNotNull()) zi_v = as<arma::vec>(NumericVector(zi));
  return vec2nv(vat::filter(as<arma::vec>(b), as<arma::vec>(a),
                          as<arma::vec>(x), zi_v));
}

// [[Rcpp::export]]
NumericVector vat_filtfilt_cpp(NumericVector b, NumericVector a, NumericVector x) {
  return vec2nv(vat::filtfilt(as<arma::vec>(b), as<arma::vec>(a), as<arma::vec>(x)));
}

// [[Rcpp::export]]
NumericVector vat_fir1_cpp(int n, NumericVector Wn, std::string type = "low") {
  return vec2nv(vat::fir1(n, as<arma::vec>(Wn), type));
}

// [[Rcpp::export]]
List vat_butter_cpp(int n, NumericVector Wn, std::string type = "low") {
  arma::vec b, a;
  vat::butter(n, as<arma::vec>(Wn), type, b, a);
  return List::create(_["b"] = wrap(b), _["a"] = wrap(a));
}

// [[Rcpp::export]]
NumericVector vat_medfilt1_cpp(NumericVector x, int n) {
  return vec2nv(vat::medfilt1(as<arma::vec>(x), n));
}

// [[Rcpp::export]]
NumericVector vat_interp1_linear_cpp(NumericVector x, NumericVector y, NumericVector xq) {
  return vec2nv(vat::interp1_linear(as<arma::vec>(x), as<arma::vec>(y), as<arma::vec>(xq)));
}

// [[Rcpp::export]]
NumericVector vat_interp1_spline_cpp(NumericVector x, NumericVector y, NumericVector xq) {
  return vec2nv(vat::interp1_spline(as<arma::vec>(x), as<arma::vec>(y), as<arma::vec>(xq)));
}

// [[Rcpp::export]]
List vat_findpeaks_cpp(NumericVector x, double min_peak_height = -1e300,
                       int min_peak_distance = 1) {
  vat::PeakResult r = vat::findpeaks(as<arma::vec>(x), min_peak_height,
                                     (arma::uword)std::max(1, min_peak_distance));
  IntegerVector locs(r.locs.n_elem);
  for (arma::uword i = 0; i < r.locs.n_elem; ++i) locs[i] = int(r.locs(i)) + 1;  // 1-indexed
  return List::create(_["locs"] = locs, _["pks"] = wrap(r.pks));
}

// [[Rcpp::export]]
NumericVector vat_resample_cpp(NumericVector x, int p, int q, double beta = 5.0) {
  return vec2nv(vat::resample(as<arma::vec>(x), p, q, beta));
}

// [[Rcpp::export]]
ComplexVector vat_fft_cpp(NumericVector x, int nfft = -1) {
  arma::cx_vec X;
  if (nfft <= 0) X = vat::fft(as<arma::vec>(x));
  else           X = vat::fft(as<arma::vec>(x), (arma::uword)nfft);
  ComplexVector out(X.n_elem);
  for (arma::uword i = 0; i < X.n_elem; ++i) {
    out[i].r = X(i).real();
    out[i].i = X(i).imag();
  }
  return out;
}

// ---------- LPC ----------

// [[Rcpp::export]]
List vat_lpcauto_cpp(NumericVector s, int p) {
  arma::vec ar; double e;
  vat::lpcauto(as<arma::vec>(s), p, ar, e);
  return List::create(_["ar"] = wrap(ar), _["e"] = e);
}

// [[Rcpp::export]]
List vat_burg_cpp(NumericVector s, int p) {
  arma::vec ar; double e;
  vat::burg(as<arma::vec>(s), p, ar, e);
  return List::create(_["ar"] = wrap(ar), _["e"] = e);
}

// [[Rcpp::export]]
NumericVector vat_lpcar2rf_cpp(NumericVector ar) {
  return vec2nv(vat::lpcar2rf(as<arma::vec>(ar)));
}

// [[Rcpp::export]]
NumericVector vat_lpcar2ra_cpp(NumericVector ar) {
  return vec2nv(vat::lpcar2ra(as<arma::vec>(ar)));
}

// [[Rcpp::export]]
NumericVector vat_lpcrf2rr_cpp(NumericVector k) {
  return vec2nv(vat::lpcrf2rr(as<arma::vec>(k)));
}

// [[Rcpp::export]]
double vat_distitar_cpp(NumericVector ar1, NumericVector ar2, bool symmetric = false) {
  return vat::distitar(as<arma::vec>(ar1), as<arma::vec>(ar2), symmetric);
}
