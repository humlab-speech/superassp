// ESTK Pitchmark C++ Implementation
// Lightweight implementation of Edinburgh Speech Tools pitchmark algorithm
// for finding glottal closure instants in laryngograph signals
//
// Based on algorithm by Mike Macon and Paul Taylor (1997)
// Centre for Speech Technology Research, University of Edinburgh
//
// This implementation works directly with R AsspDataObj structures
// for efficient in-memory processing without temporary files

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

// Simple FIR filter coefficients generator
// Creates a windowed sinc filter
std::vector<double> design_fir_filter(double cutoff_freq, int order, double sample_rate, bool highpass = false) {
  std::vector<double> h(order + 1);
  double fc = cutoff_freq / sample_rate;
  int M = order;
  int M2 = M / 2;

  // Generate windowed sinc function
  for (int i = 0; i <= M; i++) {
    int n = i - M2;
    if (n == 0) {
      h[i] = 2.0 * fc;
    } else {
      h[i] = sin(2.0 * M_PI * fc * n) / (M_PI * n);
    }
    // Apply Hamming window
    h[i] *= 0.54 - 0.46 * cos(2.0 * M_PI * i / M);
  }

  // Normalize
  double sum = 0.0;
  for (double coef : h) sum += coef;
  for (double &coef : h) coef /= sum;

  // Convert to highpass if requested
  if (highpass) {
    for (int i = 0; i <= M; i++) {
      h[i] = -h[i];
    }
    h[M2] += 1.0;
  }

  return h;
}

// Apply FIR filter with forward and backward passes (zero phase)
void fir_filter_double(std::vector<double> &signal, const std::vector<double> &filter) {
  int n = signal.size();
  int m = filter.size();
  std::vector<double> filtered(n, 0.0);

  // Forward pass
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      int idx = i - j;
      if (idx >= 0 && idx < n) {
        filtered[i] += signal[idx] * filter[j];
      }
    }
  }

  signal = filtered;
  std::fill(filtered.begin(), filtered.end(), 0.0);

  // Backward pass
  for (int i = n - 1; i >= 0; i--) {
    for (int j = 0; j < m; j++) {
      int idx = i + j;
      if (idx >= 0 && idx < n) {
        filtered[i] += signal[idx] * filter[j];
      }
    }
  }

  signal = filtered;
}

// Compute delta (differentiation) using regression
void compute_delta(const std::vector<double> &signal, std::vector<double> &delta, int regression_length = 4) {
  int n = signal.size();
  delta.resize(n);

  for (int i = 0; i < n; i++) {
    double sum = 0.0;
    int count = 0;

    for (int j = 1; j <= regression_length; j++) {
      if (i + j < n && i - j >= 0) {
        sum += j * (signal[i + j] - signal[i - j]);
        count += 2 * j * j;
      }
    }

    delta[i] = (count > 0) ? (sum / count) : 0.0;
  }
}

// Simple mean smoothing
void mean_smooth(std::vector<double> &signal, int window_size) {
  if (window_size <= 0) return;

  int n = signal.size();
  std::vector<double> smoothed(n);
  int half_window = window_size / 2;

  for (int i = 0; i < n; i++) {
    double sum = 0.0;
    int count = 0;

    for (int j = -half_window; j <= half_window; j++) {
      int idx = i + j;
      if (idx >= 0 && idx < n) {
        sum += signal[idx];
        count++;
      }
    }

    smoothed[i] = sum / count;
  }

  signal = smoothed;
}

// Find negative zero crossings
std::vector<double> find_neg_zero_crossings(const std::vector<double> &signal, double sample_rate) {
  std::vector<double> crossings;
  int n = signal.size();

  for (int i = 1; i < n; i++) {
    if (signal[i - 1] > 0.0 && signal[i] <= 0.0) {
      // Linear interpolation for sub-sample accuracy
      double t = (double)i / sample_rate;
      crossings.push_back(t);
    }
  }

  return crossings;
}

// Fill pitchmarks (interpolate and remove outliers)
std::vector<double> fill_pitchmarks(const std::vector<double> &pm,
                                      double end_time,
                                      double min_period,
                                      double max_period,
                                      double def_period) {
  std::vector<double> filled;
  if (pm.empty()) return filled;

  double last = 0.0;

  for (double current : pm) {
    if (current > end_time) break;

    double period = current - last;

    if (period < min_period) {
      // Drop this pitchmark (too close)
      continue;
    } else if (period > max_period) {
      // Interpolate
      int num = std::floor(period / def_period);
      if (num > 0) {
        double step = period / num;
        for (int i = 1; i <= num; i++) {
          filled.push_back(last + i * step);
        }
      }
    } else {
      // Keep this pitchmark
      filled.push_back(current);
    }

    last = current;
  }

  // Fill to end if needed
  if (end_time - last > max_period) {
    int num = std::floor((end_time - last) / def_period);
    if (num > 0) {
      double step = (end_time - last) / num;
      for (int i = 1; i <= num; i++) {
        filled.push_back(last + i * step);
      }
    }
  }

  return filled;
}

// Convert pitchmarks to F0 track
NumericMatrix pitchmarks_to_f0(const std::vector<double> &pm) {
  int n = pm.size();
  if (n == 0) {
    return NumericMatrix(0, 1);
  }

  NumericMatrix f0(n, 1);
  double prev = 0.0;

  for (int i = 0; i < n; i++) {
    double period = pm[i] - prev;
    f0(i, 0) = (period > 0) ? (1.0 / period) : 0.0;
    prev = pm[i];
  }

  return f0;
}

//' ESTK Pitchmark - Find glottal closure instants (C++ implementation)
//'
//' This function implements the Edinburgh Speech Tools pitchmark algorithm
//' for finding instants of glottal closure in laryngograph (EGG) waveforms.
//' This is a lightweight C++ implementation that works directly with AsspDataObj
//' structures for efficient in-memory processing.
//'
//' @param audio_obj AsspDataObj containing audio data
//' @param lx_low_frequency Low pass cutoff frequency for lx filtering (Hz)
//' @param lx_low_order Order of low pass lx filter
//' @param lx_high_frequency High pass cutoff frequency for lx filtering (Hz)
//' @param lx_high_order Order of high pass lx filter
//' @param df_low_frequency Low pass cutoff for differentiated signal (Hz)
//' @param df_low_order Order of low pass filter for differentiated signal
//' @param median_order Order of mean smoother (not median in this implementation)
//' @param fill Insert/remove pitchmarks based on min/max/def periods
//' @param min_period Minimum allowed pitch period in seconds
//' @param max_period Maximum allowed pitch period in seconds
//' @param def_period Default pitch period for unvoiced sections in seconds
//' @param invert Invert polarity of signal
//' @param to_f0 Convert pitchmarks to F0 track
//' @param verbose Show processing messages
//'
//' @return List with pitchmark times and optionally F0 values
//' @export
// [[Rcpp::export]]
List estk_pitchmark_cpp(SEXP audio_obj,
                        int lx_low_frequency = 400,
                        int lx_low_order = 19,
                        int lx_high_frequency = 40,
                        int lx_high_order = 19,
                        int df_low_frequency = 1000,
                        int df_low_order = 19,
                        int median_order = 19,
                        bool fill = false,
                        double min_period = 0.003,
                        double max_period = 0.02,
                        double def_period = 0.01,
                        bool invert = false,
                        bool to_f0 = false,
                        bool verbose = false) {

  // Extract audio data from AsspDataObj
  if (!Rf_inherits(audio_obj, "AsspDataObj")) {
    stop("Input must be an AsspDataObj");
  }

  List audio_list(audio_obj);

  // Get audio track
  if (!audio_list.containsElementNamed("audio")) {
    stop("AsspDataObj must contain 'audio' track");
  }

  NumericMatrix audio_matrix = audio_list["audio"];
  int n_samples = audio_matrix.nrow();

  // Get sample rate
  double sample_rate = as<double>(audio_list.attr("sampleRate"));
  double end_time = n_samples / sample_rate;

  if (verbose) {
    Rcout << "Processing " << n_samples << " samples at " << sample_rate << " Hz\n";
  }

  // Extract first channel to vector
  std::vector<double> signal(n_samples);
  for (int i = 0; i < n_samples; i++) {
    signal[i] = audio_matrix(i, 0);
  }

  // Invert if requested
  if (invert) {
    for (double &s : signal) s = -s;
  }

  // Step 1: Low-pass filter
  if (verbose) Rcout << "Applying low-pass filter (" << lx_low_frequency << " Hz, order " << lx_low_order << ")\n";
  auto lp_filter = design_fir_filter(lx_low_frequency, lx_low_order, sample_rate, false);
  fir_filter_double(signal, lp_filter);

  // Step 2: High-pass filter
  if (verbose) Rcout << "Applying high-pass filter (" << lx_high_frequency << " Hz, order " << lx_high_order << ")\n";
  auto hp_filter = design_fir_filter(lx_high_frequency, lx_high_order, sample_rate, true);
  fir_filter_double(signal, hp_filter);

  // Step 3: Compute delta (differentiation)
  if (verbose) Rcout << "Computing delta signal\n";
  std::vector<double> delta_signal;
  compute_delta(signal, delta_signal, 4);

  // Step 4: Optional low-pass filter on delta
  if (df_low_order > 0) {
    if (verbose) Rcout << "Applying low-pass filter to delta (" << df_low_frequency << " Hz)\n";
    auto df_lp_filter = design_fir_filter(df_low_frequency, df_low_order, sample_rate, false);
    fir_filter_double(delta_signal, df_lp_filter);
  }

  // Step 5: Mean smoothing
  if (median_order > 0) {
    if (verbose) Rcout << "Applying mean smoothing (window " << median_order << ")\n";
    mean_smooth(delta_signal, median_order);
  }

  // Step 6: Find negative zero crossings
  if (verbose) Rcout << "Finding negative zero crossings\n";
  std::vector<double> pitchmarks = find_neg_zero_crossings(delta_signal, sample_rate);

  if (verbose) Rcout << "Found " << pitchmarks.size() << " initial pitchmarks\n";

  // Step 7: Optional filling
  if (fill) {
    if (verbose) Rcout << "Filling pitchmarks (min=" << min_period << ", max=" << max_period << ", def=" << def_period << ")\n";
    pitchmarks = fill_pitchmarks(pitchmarks, end_time, min_period, max_period, def_period);
    if (verbose) Rcout << "After filling: " << pitchmarks.size() << " pitchmarks\n";
  }

  // Convert to R structures
  NumericVector pm_times = wrap(pitchmarks);

  // Create result list
  List result;
  result["pitchmarks"] = pm_times;
  result["n_pitchmarks"] = pitchmarks.size();
  result["sample_rate"] = sample_rate;
  result["duration"] = end_time;

  // Optional F0 conversion
  if (to_f0) {
    if (verbose) Rcout << "Converting pitchmarks to F0\n";
    NumericMatrix f0_track = pitchmarks_to_f0(pitchmarks);
    result["f0"] = f0_track;
  }

  return result;
}
