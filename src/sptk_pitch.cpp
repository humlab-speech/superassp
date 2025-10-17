// SPTK Pitch Extraction Wrappers for R
// Provides C++ implementations of RAPT, SWIPE, REAPER, and WORLD (DIO) algorithms

#include <Rcpp.h>
#include <vector>
#include <string>

// SPTK headers
#include "SPTK/analysis/pitch_extraction_by_rapt.h"
#include "SPTK/analysis/pitch_extraction_by_swipe.h"
#include "SPTK/analysis/pitch_extraction_by_reaper.h"
#include "SPTK/analysis/pitch_extraction_by_dio.h"

using namespace Rcpp;

//' RAPT Pitch Extraction (C++ Implementation)
//'
//' Robust Algorithm for Pitch Tracking using SPTK library.
//'
//' @param audio_obj An AsspDataObj containing audio data
//' @param minF Minimum F0 in Hz (default: 60)
//' @param maxF Maximum F0 in Hz (default: 400)
//' @param windowShift Frame shift in milliseconds (default: 10)
//' @param voicing_threshold Voicing threshold (default: 0.9)
//' @param verbose Print processing information (default: FALSE)
//' @return List with f0 (matrix), times (vector), sample_rate, n_frames
//' @export
// [[Rcpp::export]]
List rapt_cpp(SEXP audio_obj,
              double minF = 60.0,
              double maxF = 400.0,
              double windowShift = 10.0,
              double voicing_threshold = 0.9,
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
  
  // Convert audio to std::vector<double>
  std::vector<double> waveform(n_samples);
  for (int i = 0; i < n_samples; i++) {
    waveform[i] = audio_matrix(i, 0);  // First channel
  }
  
  // Calculate frame shift in samples
  int frame_shift_samples = static_cast<int>(windowShift * sample_rate / 1000.0);
  
  // Create RAPT extractor
  sptk::PitchExtractionByRapt rapt(
    frame_shift_samples,
    static_cast<double>(sample_rate),
    minF,
    maxF,
    voicing_threshold
  );
  
  if (!rapt.IsValid()) {
    stop("Failed to initialize RAPT pitch extractor");
  }
  
  // Extract pitch
  std::vector<double> f0;
  std::vector<double> epochs;  // Not used by RAPT
  sptk::PitchExtractionInterface::Polarity polarity;
  
  if (!rapt.Get(waveform, &f0, &epochs, &polarity)) {
    stop("RAPT pitch extraction failed");
  }
  
  int n_frames = f0.size();
  
  if (verbose) {
    Rcout << "Extracted " << n_frames << " F0 frames\n";
  }
  
  // Create output matrix (1 column for F0)
  NumericMatrix f0_matrix(n_frames, 1);
  NumericVector times(n_frames);
  
  for (int i = 0; i < n_frames; i++) {
    f0_matrix(i, 0) = f0[i];
    times[i] = i * windowShift / 1000.0;  // Convert ms to seconds
  }
  
  return List::create(
    Named("f0") = f0_matrix,
    Named("times") = times,
    Named("sample_rate") = sample_rate,
    Named("n_frames") = n_frames
  );
}

//' SWIPE Pitch Extraction (C++ Implementation)
//'
//' Sawtooth Waveform Inspired Pitch Estimator using SPTK library.
//'
//' @param audio_obj An AsspDataObj containing audio data
//' @param minF Minimum F0 in Hz (default: 60)
//' @param maxF Maximum F0 in Hz (default: 400)
//' @param windowShift Frame shift in milliseconds (default: 10)
//' @param voicing_threshold Voicing threshold (default: 0.3)
//' @param verbose Print processing information (default: FALSE)
//' @return List with f0 (matrix), times (vector), sample_rate, n_frames
//' @export
// [[Rcpp::export]]
List swipe_cpp(SEXP audio_obj,
               double minF = 60.0,
               double maxF = 400.0,
               double windowShift = 10.0,
               double voicing_threshold = 0.3,
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
  
  // Convert audio to std::vector<double>
  std::vector<double> waveform(n_samples);
  for (int i = 0; i < n_samples; i++) {
    waveform[i] = audio_matrix(i, 0);
  }
  
  // Calculate frame shift in samples
  int frame_shift_samples = static_cast<int>(windowShift * sample_rate / 1000.0);
  
  // Create SWIPE extractor
  sptk::PitchExtractionBySwipe swipe(
    frame_shift_samples,
    static_cast<double>(sample_rate),
    minF,
    maxF,
    voicing_threshold
  );
  
  if (!swipe.IsValid()) {
    stop("Failed to initialize SWIPE pitch extractor");
  }
  
  // Extract pitch
  std::vector<double> f0;
  std::vector<double> epochs;
  sptk::PitchExtractionInterface::Polarity polarity;
  
  if (!swipe.Get(waveform, &f0, &epochs, &polarity)) {
    stop("SWIPE pitch extraction failed");
  }
  
  int n_frames = f0.size();
  
  if (verbose) {
    Rcout << "Extracted " << n_frames << " F0 frames\n";
  }
  
  // Create output matrix
  NumericMatrix f0_matrix(n_frames, 1);
  NumericVector times(n_frames);
  
  for (int i = 0; i < n_frames; i++) {
    f0_matrix(i, 0) = f0[i];
    times[i] = i * windowShift / 1000.0;
  }
  
  return List::create(
    Named("f0") = f0_matrix,
    Named("times") = times,
    Named("sample_rate") = sample_rate,
    Named("n_frames") = n_frames
  );
}

//' REAPER Pitch Extraction (C++ Implementation)
//'
//' Robust Epoch And Pitch EstimatoR using SPTK library.
//'
//' @param audio_obj An AsspDataObj containing audio data
//' @param minF Minimum F0 in Hz (default: 60)
//' @param maxF Maximum F0 in Hz (default: 400)
//' @param windowShift Frame shift in milliseconds (default: 10)
//' @param voicing_threshold Voicing threshold (default: 0.9)
//' @param verbose Print processing information (default: FALSE)
//' @return List with f0 (matrix), times (vector), sample_rate, n_frames, epochs, polarity
//' @export
// [[Rcpp::export]]
List reaper_cpp(SEXP audio_obj,
                double minF = 60.0,
                double maxF = 400.0,
                double windowShift = 10.0,
                double voicing_threshold = 0.9,
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
  
  // Convert audio to std::vector<double>
  std::vector<double> waveform(n_samples);
  for (int i = 0; i < n_samples; i++) {
    waveform[i] = audio_matrix(i, 0);
  }
  
  // Calculate frame shift in samples
  int frame_shift_samples = static_cast<int>(windowShift * sample_rate / 1000.0);
  
  // Create REAPER extractor
  sptk::PitchExtractionByReaper reaper(
    frame_shift_samples,
    static_cast<double>(sample_rate),
    minF,
    maxF,
    voicing_threshold
  );
  
  if (!reaper.IsValid()) {
    stop("Failed to initialize REAPER pitch extractor");
  }
  
  // Extract pitch and epochs
  std::vector<double> f0;
  std::vector<double> epochs;
  sptk::PitchExtractionInterface::Polarity polarity;
  
  if (!reaper.Get(waveform, &f0, &epochs, &polarity)) {
    stop("REAPER pitch extraction failed");
  }
  
  int n_frames = f0.size();
  int n_epochs = epochs.size();
  
  if (verbose) {
    Rcout << "Extracted " << n_frames << " F0 frames and " << n_epochs << " epochs\n";
  }
  
  // Create output matrix
  NumericMatrix f0_matrix(n_frames, 1);
  NumericVector times(n_frames);
  NumericVector epoch_times(n_epochs);
  
  for (int i = 0; i < n_frames; i++) {
    f0_matrix(i, 0) = f0[i];
    times[i] = i * windowShift / 1000.0;
  }
  
  for (int i = 0; i < n_epochs; i++) {
    epoch_times[i] = epochs[i] / sample_rate;  // Convert samples to seconds
  }
  
  std::string polarity_str;
  switch (polarity) {
    case sptk::PitchExtractionInterface::kPositive:
      polarity_str = "positive";
      break;
    case sptk::PitchExtractionInterface::kNegative:
      polarity_str = "negative";
      break;
    default:
      polarity_str = "unknown";
  }
  
  return List::create(
    Named("f0") = f0_matrix,
    Named("times") = times,
    Named("sample_rate") = sample_rate,
    Named("n_frames") = n_frames,
    Named("epochs") = epoch_times,
    Named("n_epochs") = n_epochs,
    Named("polarity") = polarity_str
  );
}

//' DIO (WORLD) Pitch Extraction (C++ Implementation)
//'
//' DIO pitch extraction algorithm from WORLD vocoder using SPTK library.
//'
//' @param audio_obj An AsspDataObj containing audio data
//' @param minF Minimum F0 in Hz (default: 60)
//' @param maxF Maximum F0 in Hz (default: 400)
//' @param windowShift Frame shift in milliseconds (default: 10)
//' @param voicing_threshold Voicing threshold (default: 0.85)
//' @param verbose Print processing information (default: FALSE)
//' @return List with f0 (matrix), times (vector), sample_rate, n_frames
//' @export
// [[Rcpp::export]]
List dio_cpp(SEXP audio_obj,
             double minF = 60.0,
             double maxF = 400.0,
             double windowShift = 10.0,
             double voicing_threshold = 0.85,
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
  
  // Convert audio to std::vector<double>
  std::vector<double> waveform(n_samples);
  for (int i = 0; i < n_samples; i++) {
    waveform[i] = audio_matrix(i, 0);
  }
  
  // Calculate frame shift in samples
  int frame_shift_samples = static_cast<int>(windowShift * sample_rate / 1000.0);
  
  // Create DIO extractor
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
  
  // Extract pitch
  std::vector<double> f0;
  std::vector<double> epochs;
  sptk::PitchExtractionInterface::Polarity polarity;
  
  if (!dio.Get(waveform, &f0, &epochs, &polarity)) {
    stop("DIO pitch extraction failed");
  }
  
  int n_frames = f0.size();
  
  if (verbose) {
    Rcout << "Extracted " << n_frames << " F0 frames\n";
  }
  
  // Create output matrix
  NumericMatrix f0_matrix(n_frames, 1);
  NumericVector times(n_frames);
  
  for (int i = 0; i < n_frames; i++) {
    f0_matrix(i, 0) = f0[i];
    times[i] = i * windowShift / 1000.0;
  }
  
  return List::create(
    Named("f0") = f0_matrix,
    Named("times") = times,
    Named("sample_rate") = sample_rate,
    Named("n_frames") = n_frames
  );
}
