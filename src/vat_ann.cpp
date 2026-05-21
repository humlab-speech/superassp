// vat_ann.cpp — Generic feed-forward ANN forward pass.
// Replicates MATLAB Neural Network Toolbox semantics:
//   layer_in : x_norm = 2*(x - min)/(max - min) - 1   (mapminmax)
//   hidden   : a = tansig(IW * x_norm + b_h)
//   output   : y_norm = purelin(LW * a + b_o)
// Output denormalization is mirror (when applicable). The creak network
// trained by Kane et al. uses tansig hidden + purelin output with sigmoid-like
// output range in [-1, 1] and a decision threshold at 0.3 after median filter.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "vat_dsp.h"

using namespace Rcpp;

static inline arma::vec tansig(const arma::vec& n) {
  arma::vec out(n.n_elem);
  for (arma::uword i = 0; i < n.n_elem; ++i)
    out(i) = 2.0 / (1.0 + std::exp(-2.0 * n(i))) - 1.0;
  return out;
}

// [[Rcpp::export]]
NumericVector vat_ann_forward_cpp(NumericMatrix X_r,
                                NumericMatrix IW_r, NumericVector b_h_r,
                                NumericMatrix LW_r, NumericVector b_o_r,
                                NumericVector mini_r, NumericVector maxi_r,
                                std::string out_act = "logsig") {
  arma::mat  X  = as<arma::mat>(X_r);    // n_features x n_frames
  arma::mat  IW = as<arma::mat>(IW_r);
  arma::vec  bh = as<arma::vec>(b_h_r);
  arma::mat  LW = as<arma::mat>(LW_r);
  arma::vec  bo = as<arma::vec>(b_o_r);
  arma::vec  mini = as<arma::vec>(mini_r);
  arma::vec  maxi = as<arma::vec>(maxi_r);

  arma::uword nF = X.n_rows;
  arma::uword nT = X.n_cols;
  if (nF != mini.n_elem || nF != maxi.n_elem)
    Rcpp::stop("ann_forward: feature dim mismatch with mini/maxi");

  arma::mat Xn(X);
  for (arma::uword k = 0; k < nF; ++k) {
    double mn = mini(k), mx = maxi(k);
    double range = mx - mn;
    for (arma::uword t = 0; t < nT; ++t) {
      double v = Xn(k, t);
      if (!std::isfinite(v)) v = mn;
      Xn(k, t) = (range > 0.0) ? (-1.0 + (v - mn) / range * 2.0) : 0.0;
    }
  }

  // Hidden: a = tansig(IW * Xn + bh)
  arma::mat A(IW.n_rows, nT);
  for (arma::uword t = 0; t < nT; ++t)
    A.col(t) = tansig(IW * Xn.col(t) + bh);

  // Output: y = act(LW * A + bo)  (single output assumed)
  arma::vec Y(nT);
  for (arma::uword t = 0; t < nT; ++t) {
    double raw = arma::as_scalar(LW * A.col(t) + bo);
    if      (out_act == "logsig")  Y(t) = 1.0 / (1.0 + std::exp(-raw));
    else if (out_act == "tansig")  Y(t) = 2.0 / (1.0 + std::exp(-2.0 * raw)) - 1.0;
    else                           Y(t) = raw;  // purelin
  }

  // Median filter 3 (matches MATLAB CreakyDetection_DoClassification)
  arma::vec Yf = vat::medfilt1(Y, 3);

  NumericVector out(nT);
  for (arma::uword t = 0; t < nT; ++t) out[t] = Yf(t);
  return out;
}
