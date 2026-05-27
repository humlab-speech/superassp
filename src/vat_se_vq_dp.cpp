// se_vq_dp.cpp
// RcppArmadillo: Dynamic programming for SE-VQ GCI detection.
// Translates RESON_dyProg_mat.m (John Kane, 2011).
//
// Algorithm: Ney (1989) Viterbi-style path selection.
// Transition cost: 1 - |corr(pulse_cur, pulse_prev)| * trans_wgt
// Local cost: GCI_relAmp * rel_amp_wgt

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
IntegerVector vat_reson_dynprog_cpp(
    NumericMatrix gci_rel_amp,  // ncands x nframe
    IntegerMatrix gci_n,        // ncands x nframe (1-indexed sample positions)
    double f0_mean,
    NumericVector x,
    double fs,
    double trans_wgt,
    double rel_amp_wgt)
{
  int ncands = gci_n.nrow();
  int nframe = gci_n.ncol();
  int nx     = x.size();

  // cost(n, c) = accumulated cost at frame n for candidate c
  // Initialise with local costs = gci_rel_amp * rel_amp_wgt, transposed to nframe x ncands
  arma::mat cost(nframe, ncands);
  for (int n = 0; n < nframe; n++)
    for (int c = 0; c < ncands; c++)
      cost(n, c) = gci_rel_amp(c, n) * rel_amp_wgt;

  arma::imat prev(nframe, ncands, arma::fill::zeros);

  int pulse_half = (int)std::round(fs / f0_mean / 2.0);

  // Helper: safe slice of x (0-indexed, returns empty if degenerate)
  auto get_pulse = [&](int pos_1idx) -> arma::vec {
    int lo = pos_1idx - 1 - pulse_half;
    int hi = pos_1idx - 1 + pulse_half;
    if (lo < 0)   lo = 0;
    if (hi >= nx) hi = nx - 1;
    if (hi < lo)  return arma::vec();
    return arma::vec(x.begin() + lo, hi - lo + 1);
  };

  for (int n = 1; n < nframe; n++) {
    arma::mat costm(ncands, ncands);  // (p=row, c=col)

    for (int c = 0; c < ncands; c++) {
      arma::vec pulse_cur = get_pulse(gci_n(c, n));

      for (int p = 0; p < ncands; p++) {
        arma::vec pulse_prev = get_pulse(gci_n(p, n - 1));

        double cor = 0.0;
        if (pulse_cur.n_elem > 1 && pulse_cur.n_elem == pulse_prev.n_elem) {
          arma::mat cc = arma::cor(pulse_cur, pulse_prev);
          double v = cc(0, 0);
          if (std::isfinite(v)) cor = v;
        }
        costm(p, c) = (1.0 - std::abs(cor)) * trans_wgt;
      }
    }

    // Add cumulative costs from previous frame: costm(p,c) += cost(n-1, p)
    for (int c = 0; c < ncands; c++)
      for (int p = 0; p < ncands; p++)
        costm(p, c) += cost(n - 1, p);

    // min over rows (previous candidates) for each column (current candidate)
    for (int c = 0; c < ncands; c++) {
      arma::vec col_c = costm.col(c);
      arma::uword best_p = col_c.index_min();
      cost(n, c) += col_c(best_p);
      prev(n, c) = (int)best_p;
    }
  }

  // Traceback
  arma::ivec best(nframe);
  best(nframe - 1) = (int)cost.row(nframe - 1).index_min();
  for (int i = nframe - 1; i > 0; i--)
    best(i - 1) = prev(i, best(i));

  // Extract optimal GCI sample positions (keep 1-indexed)
  IntegerVector gci_opt(nframe);
  for (int n = 0; n < nframe; n++)
    gci_opt[n] = gci_n(best(n), n);

  return gci_opt;
}
