#include <RcppArmadillo.h>
#include <unordered_set>
#include <unordered_map>
#include <cmath>
#include <string>
#include <vector>
#include <array>
#include <algorithm>

using namespace Rcpp;

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double lz_complexity_cpp(LogicalVector S, std::string type = "exhaustive", bool normalize = true) {
  // Lempel-Ziv complexity (LZ76) for binary sequences
  // type: "exhaustive" (reproducible substring) or "primitive" (new substring)
  // normalize: divide by n/log2(n)

  int n = S.size();
  if (n < 1) return 0.0;

  double c = 0.0;  // complexity count

  if (type == "exhaustive") {
    // Exhaustive: find longest reproducible substring at each position.
    //
    // Every string ever inserted is accompanied by all of its prefixes
    // (copy step inserts S[i..i+k] for k = 0..maxlen-1; new-symbol step
    // inserts a single symbol), so the dictionary is prefix-closed. Hence
    // "longest len <= min(n-i, i) with S[i..i+len-1] in dictionary" equals
    // the depth reached by walking a binary trie along S[i..], capped at
    // min(n-i, i). This replaces the former O(n^3) substring rebuilding
    // with an O(len) walk and gives identical results.
    std::vector<std::array<int, 2> > trie;   // node 0 = root (empty string)
    trie.reserve(64);
    trie.push_back({{-1, -1}});

    int i = 0;
    while (i < n) {
      int limit = std::min(n - i, i);
      int node = 0;
      int maxlen = 0;
      while (maxlen < limit) {
        int child = trie[node][S[i + maxlen] ? 1 : 0];
        if (child < 0) break;
        node = child;
        ++maxlen;
      }

      // Both branches add one phrase. Copy step: the prefixes of
      // S[i..i+maxlen-1] are already in the trie (just walked).
      c += 1.0;
      if (maxlen > 0) {
        i += maxlen;
      } else {
        // New symbol: insert single symbol S[i]
        int b = S[i] ? 1 : 0;
        if (trie[0][b] < 0) {
          trie[0][b] = (int)trie.size();
          trie.push_back({{-1, -1}});
        }
        i += 1;
      }
    }
  } else if (type == "primitive") {
    // Primitive: new/repeated pattern words.
    // maxlen is the LARGEST len in 1..n-i whose substring is not in the
    // dictionary, so scan len downward from n-i and stop at the first miss
    // (identical result to the former upward scan over all lengths).
    std::unordered_set<std::string> dictionary;
    std::string suffix;
    suffix.reserve(n);
    for (int j = 0; j < n; ++j) suffix += (S[j] ? '1' : '0');

    int i = 0;
    while (i < n) {
      int maxlen = 0;
      for (int len = n - i; len >= 1; --len) {
        if (dictionary.find(suffix.substr(i, len)) == dictionary.end()) {
          maxlen = len;
          break;
        }
      }

      if (maxlen > 0) {
        c += 1.0;
        dictionary.insert(suffix.substr(i, maxlen));
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
