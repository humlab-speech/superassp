// iaif_full.cpp — Top-level IAIF orchestration (Alku 1992).
// Port of IAIF.m. Uses existing vat_calc_residual_cpp from iaif_lpc.cpp.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "vat_dsp.h"

using namespace Rcpp;

// Forward declaration to existing iaif_lpc.cpp helper.
List vat_calc_residual_cpp(NumericVector x_r, NumericVector x_lpc_r,
                        int ord_lpc, IntegerVector gci);

// [[Rcpp::export]]
List vat_iaif_full_cpp(NumericVector x_r, double fs, IntegerVector gci, int p = -1) {
  arma::vec x = as<arma::vec>(x_r);
  if (p <= 0) p = (int)std::round(fs / 1000.0) + 2;

  // ---- HPF: FIR order nc=704, fc=50 Hz, Hamming-windowed sinc ----
  int nc = 704;
  double fc = 50.0;
  double nfc = fc / (fs / 2.0);
  arma::vec b_hp = vat::fir1(nc, nfc, "high");
  arma::vec x_hp = vat::filter(b_hp, arma::vec({1.0}), x);

  // Align: x_filt = [x_hp(nc/2+1 : end), x(end - nc/2 + 1 : end)]
  int half = nc / 2;
  int N = x.n_elem;
  arma::vec x_filt(N);
  // First part: x_hp[half ... N-1] (N - half samples)
  for (int i = 0; i < N - half; ++i) x_filt(i) = x_hp(i + half);
  // Last part: x[N - half ... N - 1] (half samples)
  for (int i = 0; i < half; ++i) x_filt(N - half + i) = x(N - half + i);

  NumericVector xf_nv(N);
  for (int i = 0; i < N; ++i) xf_nv[i] = x_filt(i);

  // ---- Stage 1: emphasise high freqs (order 1) ----
  List st1 = vat_calc_residual_cpp(xf_nv, xf_nv, 1, gci);
  NumericVector x_emph = as<NumericVector>(st1["vector_res"]);

  // ---- Stage 2: first glottal derivative estimate (order p) ----
  List st2 = vat_calc_residual_cpp(xf_nv, x_emph, p, gci);
  NumericVector ug1 = as<NumericVector>(st2["vector_res"]);

  // ---- Stage 3: eliminate source from speech (order 4) ----
  List st3 = vat_calc_residual_cpp(xf_nv, ug1, 4, gci);
  NumericVector vt = as<NumericVector>(st3["vector_res"]);

  // ---- Stage 4: final glottal derivative (order p) ----
  List st4 = vat_calc_residual_cpp(xf_nv, vt, p, gci);
  NumericVector g_iaif = as<NumericVector>(st4["vector_res"]);

  // Integrated glottal flow (trapezoidal cumulative)
  arma::vec dg = as<arma::vec>(g_iaif);
  arma::vec g = vat::cumtrapz(dg, 1.0);

  NumericVector g_nv(g.n_elem);
  for (arma::uword i = 0; i < g.n_elem; ++i) g_nv[i] = g(i);

  return List::create(
    _["g"]      = g_nv,
    _["dg"]     = g_iaif,
    _["ar_lpc"] = st4["ar_lpc"],
    _["e_lpc"]  = st4["e_lpc"]
  );
}
