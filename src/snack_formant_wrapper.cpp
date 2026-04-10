// Snack Formant Tracker — Rcpp wrapper around snack_formant::compute_formants
//
// Calls the adapted Snack/ESPS LPC + DP formant tracking code
// (compiled from snack_formant.cc) to extract formant frequencies
// and bandwidths.

#include <Rcpp.h>
#include <vector>

// Forward declaration of the public API from snack_formant.cc
namespace snack_formant {
int compute_formants(const double *audio, int n_samples, int sample_rate,
                     int n_formants, int lpc_order, double window_dur,
                     double frame_interval, double preemphasis,
                     int lpc_type, int window_type, double ds_freq,
                     double nom_f1,
                     double **out_freqs, double **out_bands,
                     int *out_n_frames);
}

using namespace Rcpp;

//' Snack Formant Extraction (C++ implementation)
//'
//' LPC + dynamic-programming formant tracker from the Snack Sound Toolkit
//' (Talkin / AT&T / KTH).  Returns formant frequencies and bandwidths.
//'
//' @param audio_obj  An AsspDataObj containing audio data
//' @param numFormants Number of formants to track (default 4, max 7)
//' @param lpcOrder   LPC order (default 12)
//' @param windowLength Analysis window duration in seconds (default 0.049)
//' @param windowShift Frame shift in milliseconds (default 10)
//' @param preEmphasis Pre-emphasis factor (default 0.7)
//' @param dsFreq     Downsample target frequency in Hz (default 10000)
//' @param nomF1      Nominal F1 for DP cost (default -10 = use defaults)
//' @param lpcType    LPC method: 0=autocorrelation, 1=stabilized covariance, 2=covariance (default 0)
//' @param windowType Window type: 0=rectangular, 1=Hamming, 2=cos^4, 3=Hanning (default 2)
//' @param verbose    Print processing info (default FALSE)
//' @return List with fm (frequency matrix), bw (bandwidth matrix), times, sample_rate, n_frames
//' @keywords internal
// [[Rcpp::export]]
List snackf_cpp(SEXP audio_obj,
                int numFormants = 4,
                int lpcOrder = 12,
                double windowLength = 0.049,
                double windowShift = 10.0,
                double preEmphasis = 0.7,
                double dsFreq = 10000.0,
                double nomF1 = -10.0,
                int lpcType = 0,
                int windowType = 2,
                bool verbose = false) {

  if (!Rf_inherits(audio_obj, "AsspDataObj"))
    stop("Input must be an AsspDataObj");

  List audio_list(audio_obj);
  if (!audio_list.containsElementNamed("audio"))
    stop("AsspDataObj must contain 'audio' track");

  NumericMatrix audio_matrix = audio_list["audio"];
  int sample_rate = as<int>(audio_list.attr("sampleRate"));
  int n_samples   = audio_matrix.nrow();

  if (verbose) {
    Rcout << "snackf_cpp: " << n_samples << " samples @ "
          << sample_rate << " Hz, " << numFormants << " formants, LPC order "
          << lpcOrder << "\n";
  }

  // Extract first channel as double array
  std::vector<double> waveform(n_samples);
  for (int i = 0; i < n_samples; i++)
    waveform[i] = audio_matrix(i, 0);

  // Convert windowShift from ms to seconds for frame_interval
  double frame_interval = windowShift / 1000.0;

  // Allocate output pointer arrays
  std::vector<double*> freqs(numFormants, nullptr);
  std::vector<double*> bands(numFormants, nullptr);
  int n_frames = 0;

  int ok = snack_formant::compute_formants(
    waveform.data(), n_samples, sample_rate,
    numFormants, lpcOrder, windowLength, frame_interval,
    preEmphasis, lpcType, windowType, dsFreq, nomF1,
    freqs.data(), bands.data(), &n_frames
  );

  if (!ok || n_frames < 1)
    stop("Snack formant analysis failed (input too short or bad parameters)");

  if (verbose)
    Rcout << "snackf_cpp: extracted " << n_frames << " frames\n";

  // Pack into R matrices: fm is n_frames x numFormants, bw is n_frames x numFormants
  NumericMatrix fm_mat(n_frames, numFormants);
  NumericMatrix bw_mat(n_frames, numFormants);
  NumericVector times(n_frames);

  for (int j = 0; j < numFormants; j++) {
    for (int i = 0; i < n_frames; i++) {
      fm_mat(i, j) = freqs[j][i];
      bw_mat(i, j) = bands[j][i];
    }
    // Free the arrays allocated by compute_formants
    std::free(freqs[j]);
    std::free(bands[j]);
  }

  for (int i = 0; i < n_frames; i++)
    times[i] = i * windowShift / 1000.0;

  return List::create(
    Named("fm")          = fm_mat,
    Named("bw")          = bw_mat,
    Named("times")       = times,
    Named("sample_rate") = sample_rate,
    Named("n_frames")    = n_frames
  );
}
