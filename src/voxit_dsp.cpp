#include <RcppArmadillo.h>
#include <cmath>
#include <algorithm>

using namespace Rcpp;
using namespace arma;

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
NumericVector sgolay_filter_cpp(NumericVector x, int order = 2, int window = 7) {
  // Savitzky-Golay filter: degree=order, window size
  // Fixed coefficients for degree=2, window=7: [-2, 3, 6, 7, 6, 3, -2] / 21

  if (order != 2 || window != 7) {
    stop("Only degree=2, window=7 supported");
  }

  int n = x.size();
  NumericVector h = NumericVector::create(-2, 3, 6, 7, 6, 3, -2) / 21.0;
  NumericVector result(n);
  int half_window = 3;  // (7-1)/2

  // Causal/symmetric boundary extension
  for (int i = 0; i < n; ++i) {
    double sum = 0;
    for (int j = 0; j < 7; ++j) {
      int idx = i + (j - half_window);
      // Reflect at boundaries
      if (idx < 0) idx = -idx;
      if (idx >= n) idx = 2 * n - 2 - idx;
      sum += h[j] * x[idx];
    }
    result[i] = sum;
  }

  return result;
}

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List contiguous_runs_cpp(IntegerVector x, int target_value = 1) {
  // Find contiguous runs of target_value in x (1-indexed output like MATLAB)
  int n = x.size();
  IntegerVector starts, stops, lengths;

  int i = 0;
  while (i < n) {
    if (x[i] == target_value) {
      int start = i;
      while (i < n && x[i] == target_value) {
        i++;
      }
      starts.push_back(start + 1);  // 1-indexed
      stops.push_back(i);            // 1-indexed (end inclusive)
      lengths.push_back(i - start);
    } else {
      i++;
    }
  }

  return List::create(
    Named("starts") = starts,
    Named("stops") = stops,
    Named("lengths") = lengths
  );
}

// NOTE: quantile_cpp commented out — superassp already has equivalent via VAT.
// Use .vat_quantile_r() or call vat_quantile_cpp() if available.
// This is reserved for backward compat with Rvoxit if needed.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double quantile_cpp(NumericVector x, double p) {
  // Type 7 quantile (linear interpolation, matches R and MATLAB)
  NumericVector sorted_x = clone(x);
  std::sort(sorted_x.begin(), sorted_x.end());

  int n = sorted_x.size();
  double h = (n - 1) * p;
  int lower = (int)floor(h);
  int upper = (int)ceil(h);

  if (lower == upper) {
    return sorted_x[lower];
  }

  double fraction = h - lower;
  return sorted_x[lower] * (1 - fraction) + sorted_x[upper] * fraction;
}

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
IntegerVector histcounts_cpp(NumericVector x, int nbins, double lo, double hi) {
  // Fixed-width histogram bins in [lo, hi], nbins bins
  // Values exactly on right edge go to right bin (MATLAB convention)
  IntegerVector counts(nbins, 0);

  double bin_width = (hi - lo) / nbins;

  for (int i = 0; i < x.size(); ++i) {
    if (x[i] >= lo && x[i] <= hi) {
      int bin = (int)((x[i] - lo) / bin_width);
      // Handle edge case: x[i] == hi should go to last bin
      if (bin >= nbins) bin = nbins - 1;
      counts[bin]++;
    }
  }

  return counts;
}

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
NumericVector interp1_linear_cpp(NumericVector xp, NumericVector fp, NumericVector xi) {
  // Linear interpolation: match values fp at points xp, evaluate at points xi
  // (like MATLAB interp1(xp, fp, xi, 'linear'))

  int n_interp = xi.size();
  NumericVector yi(n_interp);

  for (int i = 0; i < n_interp; ++i) {
    double x_val = xi[i];

    // Find bracketing points
    int idx = 0;
    while (idx < xp.size() - 1 && xp[idx + 1] <= x_val) {
      idx++;
    }

    if (idx == xp.size() - 1) {
      // Beyond range - extrapolate with last slope
      yi[i] = fp[xp.size() - 1];
    } else if (idx == 0 && x_val < xp[0]) {
      // Before range
      yi[i] = fp[0];
    } else {
      // Linear interpolation
      double x0 = xp[idx];
      double x1 = xp[idx + 1];
      double y0 = fp[idx];
      double y1 = fp[idx + 1];

      double fraction = (x_val - x0) / (x1 - x0);
      yi[i] = y0 + fraction * (y1 - y0);
    }
  }

  return yi;
}
