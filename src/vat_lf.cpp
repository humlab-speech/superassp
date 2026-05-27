// vat_lf.cpp — LF model synthesis + Newton-Raphson area-balance solver.
// Ports lf_Area_newton.c (MEX), Rd2R.m, lf_cont.m. dyProg_LF.m (DP fit) is
// a future addition.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "vat_dsp.h"

using namespace Rcpp;

// Newton-Raphson solver for LF model alpha and epsilon parameters.
// Direct port of lf_Area_newton.c (Kane 2012).
static void lf_area_newton(double Tc, double fs, double Tp, double Te,
                           double Ta, double EE,
                           double& alpha_out, double& epsi_out) {
  const double TolFun = 1e-7;
  const int MaxIter = 100;
  double Tb = Tc - Te;
  double omega_g = arma::datum::pi / Tp;

  // Solve epsilon
  double eps0 = 1.0 / Ta;
  double change = 1.0;
  int count = 1;
  while (count <= MaxIter && std::fabs(change) > TolFun) {
    double f_eps = eps0 * Ta - 1.0 + std::exp(-eps0 * Tb);
    double f_eps_prime = Ta - Tb * std::exp(-eps0 * Tb);
    if (f_eps_prime == 0.0) break;
    change = f_eps / f_eps_prime;
    eps0 -= change;
    // Second update — matches MATLAB MEX, gives faster convergence
    eps0 -= (eps0 * Ta - 1.0 + std::exp(-eps0 * Tb)) /
            (Ta - Tb * std::exp(-eps0 * Tb));
    ++count;
  }
  epsi_out = eps0;

  // Solve alpha
  double A2 = (-EE / (eps0 * eps0 * Ta)) *
              (1.0 - std::exp(-eps0 * Tb) * (1.0 + eps0 * Tb));
  double a0 = 0.0;
  double change_a = 1.0;
  int count_a = 1;
  while (count_a <= MaxIter && std::fabs(change_a) > TolFun) {
    double r  = std::sqrt(a0 * a0 + omega_g * omega_g);
    double pA = 2.0 * std::atan((r - a0) / omega_g);
    double s1 = std::sin(omega_g * Te - pA);
    double e1 = std::exp(-a0 * Te);
    double num = r * s1 +
                 (omega_g * e1 -
                  (A2 / EE) * (a0 * a0 + omega_g * omega_g) * std::sin(omega_g * Te));
    double den = std::sin(omega_g * Te) * (1.0 - 2.0 * a0 * A2 / EE) -
                 omega_g * Te * e1;
    if (den == 0.0) break;
    change_a = num / den;
    a0 -= change_a;
    ++count_a;
  }
  alpha_out = a0;
}

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List vat_lf_area_newton_cpp(double Tc, double fs, double Tp, double Te,
                          double Ta, double EE) {
  double a, e;
  lf_area_newton(Tc, fs, Tp, Te, Ta, EE, a, e);
  return List::create(_["alpha"] = a, _["epsi"] = e);
}

// Rd -> Ra, Rk, Rg (LF "R parameter" derivation).
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List vat_rd2r_cpp(double Rd, double EE, double F0) {
  double Ra = (-1.0 + 4.8 * Rd) / 100.0;
  double Rk = (22.4 + 11.8 * Rd) / 100.0;
  double EI = (arma::datum::pi * Rk * EE) / 2.0;
  double UP = (Rd * EE) / (10.0 * F0);
  double Rg = EI / (F0 * UP * arma::datum::pi);
  return List::create(_["Ra"] = Ra, _["Rk"] = Rk, _["Rg"] = Rg);
}

// LF model pulse synthesis.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
NumericVector vat_lf_cont_cpp(double F0, double fs, double Ra, double Rk,
                           double Rg, double EE) {
  double F0min = 20, F0max = 500;
  if (F0 < F0min || F0 > F0max) return NumericVector(0);

  double T0 = 1.0 / F0;
  double Ta = Ra * T0;
  double Te = ((1.0 + Rk) / (2.0 * Rg)) * T0;
  double Tp = Te / (Rk + 1.0);
  double Tb = (1.0 - (Rk + 1.0) / (2.0 * Rg)) / F0;
  double Tc = Tb + Te;

  double alpha, epsi;
  lf_area_newton(Tc, fs, Tp, Te, Ta, EE, alpha, epsi);
  double omega = arma::datum::pi / Tp;
  double E0 = -std::fabs(EE) / (std::exp(alpha * Te) * std::sin(omega * Te));

  double dt = 1.0 / fs;
  int n_op = (int)std::floor(Te / dt);
  int n_cl = (int)std::floor(Tc / dt) - n_op;
  if (n_op <= 0 || n_cl <= 0) return NumericVector(0);

  arma::vec g(n_op + n_cl);
  for (int i = 0; i < n_op; ++i) {
    double t = (i + 1) * dt;
    g(i) = E0 * std::exp(alpha * t) * std::sin(omega * t);
  }
  for (int i = 0; i < n_cl; ++i) {
    double t = (n_op + i + 1) * dt;
    g(n_op + i) = (-EE / (epsi * Ta)) *
                   (std::exp(-epsi * (t - Te)) - std::exp(-epsi * Tb));
  }
  NumericVector out(g.n_elem);
  for (arma::uword i = 0; i < g.n_elem; ++i) out[i] = g(i);
  return out;
}
