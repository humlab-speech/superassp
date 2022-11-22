#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double vecmin(NumericVector x) {
  // Rcpp supports STL-style iterators
  NumericVector::iterator it = std::min_element(x.begin(), x.end());
  // we want the value so dereference 
  return *it;
}

