// SPTK CheapTrick Spectral Envelope Estimation Wrapper for R
// Provides C++ implementation of CheapTrick from the WORLD vocoder.

#include <Rcpp.h>
#include <vector>
#include <cmath>

// WORLD headers (sptk::world:: namespace, extended sources in src/world/)
#include "cheaptrick.h"
#include "world/common.h"
#include "world/constantnumbers.h"

using namespace Rcpp;

//' CheapTrick Spectral Envelope Estimation (C++ Implementation)
//'
//' @description Estimate the spectral envelope for each frame using the
//' CheapTrick algorithm from the WORLD vocoder (Morise 2015). Takes a
//' pre-computed F0 contour and returns a power spectral envelope matrix.
//'
//' @param audio_obj An AsspDataObj containing audio data
//' @param f0 Numeric vector of F0 values in Hz (0 for unvoiced frames)
//' @param temporal_positions Numeric vector of frame times in seconds
//'   (same length as f0)
//' @param q1 CheapTrick regularization parameter (default -0.15)
//' @param f0_floor Lower F0 bound used for FFT size calculation (default 71.0)
//' @param verbose Print processing information (default FALSE)
//' @return List with elements: spectrogram (matrix n_frames x fft_size/2+1),
//'   temporal_positions (numeric vector), sample_rate (int), n_frames (int),
//'   fft_size (int)
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List cheap_trick_cpp(SEXP audio_obj,
                     NumericVector f0,
                     NumericVector temporal_positions,
                     double q1 = -0.15,
                     double f0_floor = 71.0,
                     bool verbose = false) {

  if (!Rf_inherits(audio_obj, "AsspDataObj"))
    stop("Input must be an AsspDataObj");

  List audio_list(audio_obj);
  if (!audio_list.containsElementNamed("audio"))
    stop("AsspDataObj must contain 'audio' track");

  NumericMatrix audio_matrix = audio_list["audio"];
  int sample_rate = as<int>(audio_list.attr("sampleRate"));
  int n_samples = audio_matrix.nrow();

  int f0_length = f0.size();
  if (temporal_positions.size() != f0_length)
    stop("f0 and temporal_positions must be the same length");

  if (f0_length == 0)
    stop("f0 must have at least one element");

  NumericVector waveform_nv = audio_matrix.column(0);
  std::vector<double> waveform(waveform_nv.begin(), waveform_nv.end());
  std::vector<double> f0_vec(f0.begin(), f0.end());
  std::vector<double> tp_vec(temporal_positions.begin(), temporal_positions.end());

  sptk::world::CheapTrickOption ct_option;
  sptk::world::InitializeCheapTrickOption(sample_rate, &ct_option);
  ct_option.q1 = q1;
  ct_option.f0_floor = f0_floor;

  int fft_size = sptk::world::GetFFTSizeForCheapTrick(sample_rate, &ct_option);
  int sp_length = fft_size / 2 + 1;

  if (verbose)
    Rcout << "CheapTrick: " << f0_length << " frames, fft_size="
          << fft_size << ", sp_length=" << sp_length << "\n";

  // Allocate output spectrogram — use vector-of-vector to guarantee lifetime
  std::vector<std::vector<double>> spec_data(f0_length,
                                             std::vector<double>(sp_length, 0.0));
  std::vector<double*> spectrogram_ptrs(f0_length);
  for (int i = 0; i < f0_length; i++)
    spectrogram_ptrs[i] = spec_data[i].data();

  sptk::world::CheapTrick(waveform.data(), n_samples, sample_rate,
                           tp_vec.data(), f0_vec.data(), f0_length,
                           &ct_option, spectrogram_ptrs.data());

  // Flatten to R matrix (n_frames rows x sp_length cols)
  NumericMatrix sp_matrix(f0_length, sp_length);
  for (int i = 0; i < f0_length; i++)
    for (int j = 0; j < sp_length; j++)
      sp_matrix(i, j) = spec_data[i][j];

  NumericVector times_out(f0_length);
  for (int i = 0; i < f0_length; i++)
    times_out[i] = tp_vec[i];

  return List::create(
    Named("spectrogram")         = sp_matrix,
    Named("temporal_positions")  = times_out,
    Named("sample_rate")         = sample_rate,
    Named("n_frames")            = f0_length,
    Named("fft_size")            = fft_size
  );
}
