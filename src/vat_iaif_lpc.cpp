// iaif_lpc.cpp
// RcppArmadillo: Per-GCI LPC analysis and inverse filtering for IAIF.
// Translates calc_residual.m (John Kane, 2012).
//
// For each GCI: extract 2*T0 window, Hamming window, Levinson-Durbin LPC,
// FIR inverse filter with initial conditions, overlap-add residual.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

static void levinson_durbin_iaif(const arma::vec& s_win, int p,
                                  arma::vec& ar, double& e)
{
  int n = s_win.n_elem;
  ar.zeros(p + 1);
  ar(0) = 1.0;

  arma::vec r(p + 1, arma::fill::zeros);
  for (int k = 0; k <= p; k++)
    for (int i = 0; i < n - k; i++)
      r(k) += s_win(i) * s_win(i + k);

  if (r(0) == 0.0) { e = 0.0; return; }

  arma::vec a(p, arma::fill::zeros);
  e = r(0);
  for (int m = 0; m < p; m++) {
    double num = r(m + 1);
    for (int j = 0; j < m; j++)
      num += a(j) * r(m - j);
    double k = -num / e;
    arma::vec a_new = a;
    for (int j = 0; j < m; j++)
      a_new(j) += k * a(m - 1 - j);
    a_new(m) = k;
    a = a_new;
    e *= (1.0 - k * k);
    if (e <= 0.0) { e = 1e-10; break; }
  }
  for (int j = 0; j < p; j++) ar(j + 1) = a(j);
}

// FIR filter with initial state (state = last p samples of previous output/input)
static arma::vec fir_filter_ic(const arma::vec& b, const arma::vec& x,
                                 const arma::vec& zi)
{
  int nx = x.n_elem;
  int nb = b.n_elem;
  arma::vec y(nx, arma::fill::zeros);
  for (int n = 0; n < nx; n++) {
    double val = b(0) * x(n);
    for (int k = 1; k < nb; k++) {
      int xi = n - k;
      if (xi >= 0)
        val += b(k) * x(xi);
      else if ((int)zi.n_elem > -(xi) - 1)
        val += b(k) * zi(zi.n_elem + xi);
    }
    y(n) = val;
  }
  return y;
}

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List vat_calc_residual_cpp(NumericVector x_r,
                        NumericVector x_lpc_r,
                        int ord_lpc,
                        IntegerVector gci)
{
  arma::vec x     = Rcpp::as<arma::vec>(x_r);
  arma::vec x_lpc = Rcpp::as<arma::vec>(x_lpc_r);
  int N      = x.n_elem;
  int n_gci  = gci.size();

  arma::vec vector_res(N, arma::fill::zeros);
  arma::mat ar_lpc(ord_lpc + 1, n_gci, arma::fill::zeros);
  arma::vec e_lpc(n_gci, arma::fill::zeros);

  arma::vec ze_lpc(ord_lpc, arma::fill::zeros);  // initial conditions
  bool has_prev = false;
  arma::vec prev_frame_res;
  arma::vec prev_residual;

  for (int i = 0; i < n_gci; i++) {
    int gci_pos = gci[i] - 1;  // 0-indexed

    int T0_cur;
    if (i > 0)
      T0_cur = gci[i] - gci[i - 1];
    else if (i < n_gci - 1)
      T0_cur = gci[i + 1] - gci[i];
    else
      continue;

    int start = gci_pos - T0_cur;
    int stop  = gci_pos + T0_cur;

    if (start < 0 || stop >= N) continue;
    if (stop - start + 1 <= (int)(ord_lpc * 1.5)) continue;

    // LPC analysis on x_lpc window
    int wlen = stop - start + 1;
    arma::vec hann_win = 0.5 - 0.5 * arma::cos(
        2.0 * arma::datum::pi * arma::linspace<arma::vec>(0, wlen - 1, wlen)
        / (wlen - 1));
    arma::vec frame_lpc = x_lpc.subvec(start, stop) % hann_win;

    arma::vec ar;
    double e;
    levinson_durbin_iaif(frame_lpc, ord_lpc, ar, e);
    ar = arma::real(ar);

    ar_lpc.col(i) = ar;
    e_lpc(i) = e;

    // Inverse filtering with initial conditions
    arma::vec frame_res = x.subvec(start, stop);
    arma::vec residual;
    if (has_prev) {
      residual = fir_filter_ic(ar, frame_res, ze_lpc);
    } else {
      // No initial conditions: plain FIR
      int nx = frame_res.n_elem;
      int nb = ar.n_elem;
      residual.zeros(nx);
      for (int n = 0; n < nx; n++)
        for (int k = 0; k < nb && k <= n; k++)
          residual(n) += ar(k) * frame_res(n - k);
    }

    // Scale by 1/sqrt(e) (MATLAB: filter(ar, sqrt(e), ...))
    if (e > 0.0) residual /= std::sqrt(e);

    // Hamming window the residual
    arma::vec ham_win = 0.54 - 0.46 * arma::cos(
        2.0 * arma::datum::pi * arma::linspace<arma::vec>(0, wlen - 1, wlen)
        / (wlen - 1));
    arma::vec residual_win = residual % ham_win;

    vector_res.subvec(start, stop) += residual_win;

    // Update initial conditions for next frame (last ord_lpc samples)
    if (wlen >= ord_lpc) {
      ze_lpc = frame_res.subvec(wlen - ord_lpc, wlen - 1);
    }
    prev_frame_res = frame_res;
    prev_residual  = residual;
    has_prev = true;
  }

  return List::create(
    Named("vector_res") = Rcpp::wrap(vector_res),
    Named("ar_lpc")     = Rcpp::wrap(ar_lpc),
    Named("e_lpc")      = Rcpp::wrap(e_lpc)
  );
}
