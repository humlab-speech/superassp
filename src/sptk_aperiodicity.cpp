// SPTK Aperiodicity Estimation Wrapper for R
// Provides C++ implementation of D4C algorithm from WORLD vocoder

#include <Rcpp.h>
#include <vector>
#include <string>
#include <cmath>

// SPTK headers
#include "SPTK/analysis/pitch_extraction_by_dio.h"
#include "SPTK/third_party/WORLD/world/d4c.h"
#include "SPTK/third_party/WORLD/world/common.h"

using namespace Rcpp;

//' D4C Aperiodicity Estimation (C++ Implementation)
//'
//' Estimates band aperiodicity using D4C algorithm from WORLD vocoder.
//' This implementation calls the SPTK C++ library directly for optimal performance.
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
  
  // Extract audio data from AsspDataObj
  if (!Rf_inherits(audio_obj, "AsspDataObj")) {
    stop("Input must be an AsspDataObj");
  }
  
  List audio_list(audio_obj);
  if (!audio_list.containsElementNamed("audio")) {
    stop("AsspDataObj must contain 'audio' track");
  }
  
  NumericMatrix audio_matrix = audio_list["audio"];
  int sample_rate = as<int>(audio_list.attr("sampleRate"));
  int n_samples = audio_matrix.nrow();
  
  if (verbose) {
    Rcout << "Processing audio: " << n_samples << " samples at " << sample_rate << " Hz\n";
    Rcout << "F0 range: " << minF << " - " << maxF << " Hz\n";
    Rcout << "Window shift: " << windowShift << " ms\n";
  }
  
  // Convert audio to std::vector<double> for pitch extraction
  std::vector<double> waveform(n_samples);
  for (int i = 0; i < n_samples; i++) {
    waveform[i] = audio_matrix(i, 0);  // First channel
  }
  
  // Step 1: Extract F0 using DIO
  int frame_shift_samples = static_cast<int>(windowShift * sample_rate / 1000.0);
  
  sptk::PitchExtractionByDio dio(
    frame_shift_samples,
    static_cast<double>(sample_rate),
    minF,
    maxF,
    voicing_threshold
  );
  
  if (!dio.IsValid()) {
    stop("Failed to initialize DIO pitch extractor");
  }
  
  std::vector<double> f0;
  std::vector<double> epochs;  // Not used
  sptk::PitchExtractionInterface::Polarity polarity;
  
  if (!dio.Get(waveform, &f0, &epochs, &polarity)) {
    stop("DIO pitch extraction failed");
  }
  
  int n_frames = f0.size();
  
  if (verbose) {
    Rcout << "Extracted " << n_frames << " F0 frames\n";
  }
  
  // Step 2: Prepare temporal positions
  std::vector<double> temporal_positions(n_frames);
  for (int i = 0; i < n_frames; i++) {
    temporal_positions[i] = i * windowShift / 1000.0;  // Convert ms to seconds
  }
  
  // Step 3: Calculate FFT size
  // Use a suitable FFT size based on sample rate
  int fft_size = static_cast<int>(std::pow(2.0, 1.0 + 
    static_cast<int>(std::log(3.0 * sample_rate / 71.0 + 1) / std::log(2.0))));
  
  if (verbose) {
    Rcout << "FFT size: " << fft_size << "\n";
  }
  
  // Step 4: Initialize D4C option
  sptk::world::D4COption option;
  sptk::world::InitializeD4COption(&option);
  option.threshold = threshold;
  
  // Step 5: Allocate aperiodicity matrix
  int n_bins = fft_size / 2 + 1;
  double **aperiodicity = new double*[n_frames];
  for (int i = 0; i < n_frames; i++) {
    aperiodicity[i] = new double[n_bins];
  }
  
  // Step 6: Run D4C
  // Convert waveform to raw pointer for C API
  const double* x_ptr = waveform.data();
  const double* f0_ptr = f0.data();
  const double* time_ptr = temporal_positions.data();
  
  sptk::world::D4C(x_ptr, n_samples, sample_rate,
                   time_ptr, f0_ptr, n_frames,
                   fft_size, &option, aperiodicity);
  
  if (verbose) {
    Rcout << "D4C aperiodicity estimation complete\n";
  }
  
  // Step 7: Convert to R matrix
  NumericMatrix ap_matrix(n_frames, n_bins);
  for (int i = 0; i < n_frames; i++) {
    for (int j = 0; j < n_bins; j++) {
      ap_matrix(i, j) = aperiodicity[i][j];
    }
  }
  
  // Step 8: Create time and f0 vectors for output
  NumericVector times(n_frames);
  NumericVector f0_vec(n_frames);
  for (int i = 0; i < n_frames; i++) {
    times[i] = temporal_positions[i];
    f0_vec[i] = f0[i];
  }
  
  // Clean up
  for (int i = 0; i < n_frames; i++) {
    delete[] aperiodicity[i];
  }
  delete[] aperiodicity;
  
  return List::create(
    Named("aperiodicity") = ap_matrix,
    Named("times") = times,
    Named("f0") = f0_vec,
    Named("sample_rate") = sample_rate,
    Named("n_frames") = n_frames,
    Named("fft_size") = fft_size
  );
}
