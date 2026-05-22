#include <RcppArmadillo.h>
#include <unordered_set>
#include <unordered_map>
#include <cmath>
#include <string>

using namespace Rcpp;

// [[Rcpp::export]]
double lz_complexity_cpp(LogicalVector S, std::string type = "exhaustive", bool normalize = true) {
  // Lempel-Ziv complexity (LZ76) for binary sequences
  // type: "exhaustive" (reproducible substring) or "primitive" (new substring)
  // normalize: divide by n/log2(n)

  int n = S.size();
  if (n < 1) return 0.0;

  std::unordered_set<std::string> dictionary;
  double c = 0.0;  // complexity count

  if (type == "exhaustive") {
    // Exhaustive: find longest reproducible substring at each position
    int i = 0;
    while (i < n) {
      // Find longest match in dictionary
      int maxlen = 0;
      for (int len = 1; len <= n - i && len <= i; ++len) {
        std::string substr = "";
        for (int j = 0; j < len; ++j) {
          substr += (S[i + j] ? '1' : '0');
        }
        if (dictionary.find(substr) != dictionary.end()) {
          maxlen = len;
        }
      }

      if (maxlen > 0) {
        // Copy: output the found substring pointer
        c += 1.0;
        for (int j = 0; j < maxlen; ++j) {
          std::string substr = "";
          for (int k = 0; k <= j; ++k) {
            substr += (S[i + k] ? '1' : '0');
          }
          dictionary.insert(substr);
        }
        i += maxlen;
      } else {
        // New symbol
        c += 1.0;
        std::string single = "";
        single += (S[i] ? '1' : '0');
        dictionary.insert(single);
        i += 1;
      }
    }
  } else if (type == "primitive") {
    // Primitive: new/repeated pattern words
    int i = 0;
    while (i < n) {
      // Find longest new substring
      int maxlen = 0;
      bool found = false;
      for (int len = 1; len <= n - i; ++len) {
        std::string substr = "";
        for (int j = 0; j < len; ++j) {
          substr += (S[i + j] ? '1' : '0');
        }
        if (dictionary.find(substr) == dictionary.end()) {
          maxlen = len;
          found = true;
        }
      }

      if (found && maxlen > 0) {
        c += 1.0;
        std::string substr = "";
        for (int j = 0; j < maxlen; ++j) {
          substr += (S[i + j] ? '1' : '0');
        }
        dictionary.insert(substr);
        i += maxlen;
      } else {
        // Already in dictionary, use it
        i += 1;
      }
    }
  } else {
    stop("Unknown type: use 'exhaustive' or 'primitive'");
  }

  if (normalize) {
    // Normalize: c / (n / log2(n))
    double denom = n / log2(n);
    return 100.0 * c / denom;  // Scale to 0-100 range
  }

  return c;
}
