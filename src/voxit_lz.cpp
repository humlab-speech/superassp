#include <RcppArmadillo.h>
#include <cmath>
#include <string>
#include <vector>

using namespace Rcpp;

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double lz_complexity_cpp(LogicalVector S, std::string type = "exhaustive", bool normalize = true) {
  // Lempel-Ziv (1976) complexity of a binary sequence: the number of
  // components in the exhaustive production history of S. This is what
  // Voxit computes with Quang Thai's calc_lz_complexity(S, 'exhaustive', 1)
  // (MATLAB File Exchange 38211, embedded in voxitAnalysis.m), which derives
  // the history from the LZ76 eigenfunction. The count is obtained here with
  // the Kaspar & Schuster (1987) scan, which yields the same number in
  // O(n * c) instead of rebuilding substrings (pinned against a transcription
  // of the MATLAB code in test-voxit-helpers.R). NA counts as TRUE.
  //
  // Voxit reports 100 * C / (n / log2(n)); normalize = TRUE returns that.
  if (type != "exhaustive") stop("Unknown type: only 'exhaustive' is supported");

  const int n = S.size();
  if (n < 1) return 0.0;

  std::vector<unsigned char> b(n);
  for (int j = 0; j < n; ++j) b[j] = (S[j] != 0);   // NA_LOGICAL is non-zero

  int c = 1;
  if (n > 1) {
    int l = 1, i = 0, k = 1, kmax = 1;
    while (true) {
      if (b[i + k - 1] == b[l + k - 1]) {
        ++k;
        if (l + k > n) { ++c; break; }
      } else {
        if (k > kmax) kmax = k;
        ++i;
        if (i == l) {             // no earlier start reproduces more: new component
          ++c;
          l += kmax;
          if (l + 1 > n) break;
          i = 0; k = 1; kmax = 1;
        } else {
          k = 1;
        }
      }
    }
  }

  if (normalize) {
    double denom = n / log2(n);     // Inf for n = 1, giving 0 as in MATLAB
    return 100.0 * c / denom;
  }
  return c;
}
