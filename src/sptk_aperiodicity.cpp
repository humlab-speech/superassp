// SPTK Aperiodicity Estimation Wrapper for R
// Provides C++ implementation of D4C algorithm from WORLD vocoder

#include <Rcpp.h>
#include <vector>
#include <string>
#include <cmath>

// SPTK headers
// Note: D4C implementation requires full WORLD vocoder headers which are not
// currently available in the SPTK submodule. This function is deprecated.
// Users should use the Python-based aperiodicities() function instead.
#include "SPTK/analysis/pitch_extraction_by_world.h"
// #include "SPTK/third_party/WORLD/world/d4c.h"
// #include "SPTK/third_party/WORLD/world/common.h"

using namespace Rcpp;

//' D4C Aperiodicity Estimation (C++ Implementation)
//'
//' @description **DEPRECATED**: This function is currently unavailable due to
//' missing WORLD vocoder headers in the SPTK submodule. Please use the
//' Python-based `trk_aperiodicities()` function instead.
//'
//' @param audio_obj An AsspDataObj containing audio data
//' @param minF Minimum F0 in Hz (default: 60)
//' @param maxF Maximum F0 in Hz (default: 400)
//' @param windowShift Frame shift in milliseconds (default: 5)
//' @param voicing_threshold Voicing threshold for F0 detection (default: 0.85)
//' @param threshold D4C threshold parameter (default: 0.85)
//' @param verbose Print processing information (default: FALSE)
//' @return List with aperiodicity (matrix), times (vector), f0 (vector), sample_rate, n_frames, fft_size
//' @export
// [[Rcpp::export]]
List d4c_cpp(SEXP audio_obj,
             double minF = 60.0,
             double maxF = 400.0,
             double windowShift = 5.0,
             double voicing_threshold = 0.85,
             double threshold = 0.85,
             bool verbose = false) {
  
  stop("d4c_cpp() is currently unavailable. Please use trk_aperiodicities() instead.");
  
  return List::create(); // Never reached but needed for compilation
}
