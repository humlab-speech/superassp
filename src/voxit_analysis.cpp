#include <RcppArmadillo.h>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
NumericVector compute_f0_stats_simple_cpp(NumericVector f0, NumericVector vuv) {
  // Simple F0 statistics: mean, range, entropy placeholder

  NumericVector f0_voiced;
  for (int i = 0; i < f0.size(); ++i) {
    if (vuv[i] > 0.5) {
      f0_voiced.push_back(f0[i]);
    }
  }

  NumericVector stats(3);

  if (f0_voiced.size() == 0) {
    stats[0] = NA_REAL;
    stats[1] = NA_REAL;
    stats[2] = NA_REAL;
    return stats;
  }

  // Geometric mean
  double log2_sum = 0;
  for (double f : f0_voiced) {
    log2_sum += log(f) / log(2.0);
  }
  stats[0] = pow(2.0, log2_sum / f0_voiced.size());

  // Range in octaves
  double f0_min = *std::min_element(f0_voiced.begin(), f0_voiced.end());
  double f0_max = *std::max_element(f0_voiced.begin(), f0_voiced.end());
  stats[1] = log(f0_max / f0_min) / log(2.0);

  // Voicing percent
  stats[2] = 100.0 * f0_voiced.size() / f0.size();

  return stats;
}

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double compute_pause_count_cpp(NumericVector lin_power, NumericVector vuv,
                               double thresh_db = 10.0) {
  // Simple pause counter

  // Get 90th percentile
  NumericVector sorted = clone(lin_power);
  std::sort(sorted.begin(), sorted.end());
  int idx_90 = (int)(0.9 * sorted.size());
  double p90 = sorted[idx_90];

  double log_p90 = 10.0 * log10(p90 + 1e-10);
  double thresh = log_p90 - thresh_db;

  int pause_count = 0;
  for (int i = 0; i < lin_power.size(); ++i) {
    double log_power = 10.0 * log10(lin_power[i] + 1e-10);
    if (log_power < thresh && vuv[i] < 0.5) {
      pause_count++;
    }
  }

  return pause_count;
}

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double compute_intensity_mean_cpp(NumericVector lin_power) {
  // Mean log power (dB)

  double sum_log = 0;
  for (double p : lin_power) {
    sum_log += 10.0 * log10(p + 1e-10);
  }

  return sum_log / lin_power.size();
}
