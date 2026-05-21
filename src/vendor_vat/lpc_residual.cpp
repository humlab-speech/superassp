// lpc_residual.cpp
// RcppArmadillo: LP residual computation via overlap-add framing.
// Translates GetLPCresidual.m (Thomas Drugman / John Kane).
//
// Each frame: Hanning window -> Levinson-Durbin LPC -> FIR inverse filter
// -> energy normalisation -> overlap-add into output.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// Levinson-Durbin autocorrelation LPC.
// Returns AR vector (length p+1, first element = 1) and residual energy.
static void levinson_durbin(const arma::vec& s_win, int p,
                             arma::vec& ar, double& e)
{
  int n = s_win.n_elem;
  ar.zeros(p + 1);
  ar(0) = 1.0;

  // Autocorrelation
  arma::vec r(p + 1, arma::fill::zeros);
  for (int k = 0; k <= p; k++) {
    for (int i = 0; i < n - k; i++)
      r(k) += s_win(i) * s_win(i + k);
  }

  if (r(0) == 0.0) { e = 0.0; return; }

  arma::vec a(p, arma::fill::zeros);
  e = r(0);

  for (int m = 0; m < p; m++) {
    double num = r(m + 1);
    for (int j = 0; j < m; j++)
      num += a(j) * r(m - j);
    double k = -num / e;

    // Update a
    arma::vec a_new = a;
    for (int j = 0; j < m; j++)
      a_new(j) += k * a(m - 1 - j);
    a_new(m) = k;
    a = a_new;

    e *= (1.0 - k * k);
    if (e <= 0.0) { e = 1e-10; break; }
  }

  for (int j = 0; j < p; j++)
    ar(j + 1) = a(j);
}

// FIR filter: y[n] = sum_k b[k] * x[n-k]  (no feedback)
static arma::vec fir_filter(const arma::vec& b, const arma::vec& x)
{
  int nx = x.n_elem;
  int nb = b.n_elem;
  arma::vec y(nx, arma::fill::zeros);
  for (int n = 0; n < nx; n++)
    for (int k = 0; k < nb && k <= n; k++)
      y(n) += b(k) * x(n - k);
  return y;
}

// [[Rcpp::export]]
NumericVector vat_lpc_residual_cpp(NumericVector wave_r,
                                int window_len,
                                int shift,
                                int order)
{
  arma::vec wave = Rcpp::as<arma::vec>(wave_r);
  int N = wave.n_elem;

  arma::vec res(N, arma::fill::zeros);
  arma::vec hann_win = 0.5 - 0.5 * arma::cos(
      2.0 * arma::datum::pi * arma::linspace<arma::vec>(0, window_len - 1, window_len)
      / (window_len - 1));

  int start = 0;
  int stop  = start + window_len;

  while (stop < N) {
    arma::vec segment = wave.subvec(start, stop - 1) % hann_win;

    arma::vec ar;
    double e;
    levinson_durbin(segment, order, ar, e);

    arma::vec inv = fir_filter(ar, segment);

    // Energy normalisation
    double e_seg = arma::dot(segment, segment);
    double e_inv = arma::dot(inv, inv);
    if (e_inv > 0.0)
      inv *= std::sqrt(e_seg / e_inv);

    res.subvec(start, stop - 1) += inv;

    start += shift;
    stop   = start + window_len;
  }

  // Normalize
  double mx = arma::max(arma::abs(res));
  if (mx > 0.0) res /= mx;

  return Rcpp::wrap(res);
}
