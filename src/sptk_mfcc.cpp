// SPTK MFCC Extraction Wrapper for R
// Provides C++ implementation of MFCC analysis using SPTK library

#include <Rcpp.h>
#include <vector>
#include <cmath>

// SPTK headers
#include "SPTK/analysis/mel_frequency_cepstral_coefficients_analysis.h"
#include "SPTK/conversion/waveform_to_spectrum.h"
#include "SPTK/utils/sptk_utils.h"

using namespace Rcpp;

//' SPTK MFCC Extraction (C++ Implementation)
//'
//' Extract Mel-Frequency Cepstral Coefficients using SPTK library.
//' This is a high-performance C++ implementation that is significantly
//' faster than Python-based implementations.
//'
//' @param audio_obj An AsspDataObj containing audio data
//' @param n_mfcc Number of MFCC coefficients (default: 13)
//' @param n_mels Number of mel filterbanks (default: 40)
//' @param windowShift Frame shift in milliseconds (default: 10.0)
//' @param windowSize Window size in milliseconds (default: 25.0)
//' @param fmin Minimum frequency in Hz (default: 0.0)
//' @param fmax Maximum frequency in Hz (default: sample_rate/2)
//' @param lifter Liftering coefficient (default: 22)
//' @param floor Floor value for mel filterbank output (default: 1.0)
//' @param verbose Print processing information (default: FALSE)
//' @return List with mfcc (matrix), times (vector), sample_rate, n_frames
//' @export
// [[Rcpp::export]]
List sptk_mfcc_cpp(SEXP audio_obj,
                   int n_mfcc = 13,
                   int n_mels = 40,
                   double windowShift = 10.0,
                   double windowSize = 25.0,
                   double fmin = 0.0,
                   double fmax = 0.0,
                   int lifter = 22,
                   double floor = 1.0,
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
    Rcout << "Extracting " << n_mfcc << " MFCCs with " << n_mels << " mel channels\n";
    Rcout << "Window size: " << windowSize << " ms, shift: " << windowShift << " ms\n";
  }
  
  // Validate parameters
  if (n_mfcc <= 0 || n_mfcc >= n_mels) {
    stop("n_mfcc must be positive and less than n_mels");
  }
  
  if (n_mels <= 0) {
    stop("n_mels must be positive");
  }
  
  // Set fmax to Nyquist frequency if not specified
  if (fmax <= 0.0 || fmax > sample_rate / 2.0) {
    fmax = sample_rate / 2.0;
  }
  
  if (fmin < 0.0 || fmin >= fmax) {
    stop("fmin must be non-negative and less than fmax");
  }
  
  // Calculate FFT length (must be power of 2)
  int window_length_samples = static_cast<int>(windowSize * sample_rate / 1000.0);
  int fft_length = 1;
  while (fft_length < window_length_samples) {
    fft_length *= 2;
  }
  
  // Calculate frame shift in samples
  int frame_shift_samples = static_cast<int>(windowShift * sample_rate / 1000.0);
  
  // Calculate number of frames
  int n_frames = (n_samples - window_length_samples) / frame_shift_samples + 1;
  if (n_frames <= 0) {
    stop("Audio is too short for the given window parameters");
  }
  
  if (verbose) {
    Rcout << "FFT length: " << fft_length << " samples\n";
    Rcout << "Window length: " << window_length_samples << " samples\n";
    Rcout << "Frame shift: " << frame_shift_samples << " samples\n";
    Rcout << "Number of frames: " << n_frames << "\n";
  }
  
  // Create MFCC analyzer
  sptk::MelFrequencyCepstralCoefficientsAnalysis mfcc_analyzer(
    fft_length,
    n_mels,
    n_mfcc,
    lifter,
    static_cast<double>(sample_rate),
    fmin,
    fmax,
    floor
  );
  
  if (!mfcc_analyzer.IsValid()) {
    stop("Failed to initialize MFCC analyzer");
  }
  
  // Create waveform to spectrum converter
  sptk::WaveformToSpectrum waveform_to_spectrum(
    fft_length,
    fft_length,
    sptk::SpectrumToSpectrum::InputOutputFormats::kPowerSpectrum,
    0.0,
    -DBL_MAX
  );
  
  if (!waveform_to_spectrum.IsValid()) {
    stop("Failed to initialize spectrum converter");
  }
  
  // Prepare buffers
  sptk::WaveformToSpectrum::Buffer spectrum_buffer;
  sptk::MelFrequencyCepstralCoefficientsAnalysis::Buffer mfcc_buffer;
  
  // Convert audio to std::vector<double>
  std::vector<double> waveform(n_samples);
  for (int i = 0; i < n_samples; i++) {
    waveform[i] = audio_matrix(i, 0);  // First channel
  }
  
  // Output matrices
  NumericMatrix mfcc_matrix(n_frames, n_mfcc + 1);  // +1 for c0
  NumericVector times(n_frames);
  
  // Process each frame
  std::vector<double> frame(fft_length, 0.0);
  std::vector<double> power_spectrum(fft_length / 2 + 1);
  std::vector<double> mfcc_coeffs(n_mfcc + 1);
  double energy;
  
  for (int frame_idx = 0; frame_idx < n_frames; frame_idx++) {
    int start_sample = frame_idx * frame_shift_samples;
    
    // Extract frame with zero padding if needed
    for (int i = 0; i < fft_length; i++) {
      if (start_sample + i < n_samples && i < window_length_samples) {
        frame[i] = waveform[start_sample + i];
      } else {
        frame[i] = 0.0;
      }
    }
    
    // Compute power spectrum
    if (!waveform_to_spectrum.Run(frame, &power_spectrum, &spectrum_buffer)) {
      stop("Failed to compute power spectrum at frame " + std::to_string(frame_idx));
    }
    
    // Extract MFCCs
    if (!mfcc_analyzer.Run(power_spectrum, &mfcc_coeffs, &energy, &mfcc_buffer)) {
      stop("Failed to extract MFCCs at frame " + std::to_string(frame_idx));
    }
    
    // Store results (skip c0 at index 0, store c1-cn)
    for (int coef = 0; coef <= n_mfcc; coef++) {
      mfcc_matrix(frame_idx, coef) = mfcc_coeffs[coef];
    }
    
    // Calculate frame time (center of frame)
    times[frame_idx] = (start_sample + window_length_samples / 2.0) / sample_rate;
  }
  
  if (verbose) {
    Rcout << "Successfully extracted " << n_frames << " MFCC frames\n";
  }
  
  return List::create(
    Named("mfcc") = mfcc_matrix,
    Named("times") = times,
    Named("sample_rate") = sample_rate,
    Named("n_frames") = n_frames,
    Named("n_mfcc") = n_mfcc + 1  // Including c0
  );
}
