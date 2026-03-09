// gfmiaif.cpp — GFM-IAIF (Glottal Flow Model-based IAIF) in C++
//
// Port of inst/python/gfmiaif/gfmiaif.py and lpc.py
//
// Reference: O. Perrotin and I. V. McLoughlin (2019)
//   "A spectral glottal flow model for source-filter separation of speech",
//   IEEE ICASSP 2019, pp. 7160-7164.

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include "dsp_helpers.hpp"

using namespace Rcpp;

// ============================================================
// DSP helpers shared via dsp_helpers.hpp:
//   levinson_durbin, fir_filter, leaky_integrate
// Local helpers below are specific to GFM-IAIF.
// ============================================================

namespace gfm {

// Biased autocorrelation (note: divides by N, unlike iaif.cpp's unbiased version)
static std::vector<double> autocorr(const double* x, int N, int max_lag) {
  std::vector<double> r(max_lag + 1, 0.0);
  for (int k = 0; k <= max_lag && k < N; k++) {
    double sum = 0.0;
    for (int n = 0; n < N - k; n++) sum += x[n] * x[n + k];
    r[k] = sum / N;
  }
  return r;
}

// LPC: biased autocorrelation + Levinson-Durbin
static std::vector<double> lpc(const double* x, int N, int order) {
  auto r = autocorr(x, N, order + 1);
  return levinson_durbin(r, order);
}

// Hanning window
static std::vector<double> hanning(int N) {
  std::vector<double> w(N);
  for (int n = 0; n < N; n++)
    w[n] = 0.5 * (1.0 - std::cos(2.0 * M_PI * n / (N - 1)));
  return w;
}

// Hamming window
static std::vector<double> hamming(int N) {
  std::vector<double> w(N);
  for (int n = 0; n < N; n++)
    w[n] = 0.54 - 0.46 * std::cos(2.0 * M_PI * n / (N - 1));
  return w;
}

// Blackman window
static std::vector<double> blackman(int N) {
  std::vector<double> w(N);
  for (int n = 0; n < N; n++) {
    double t = 2.0 * M_PI * n / (N - 1);
    w[n] = 0.42 - 0.5 * std::cos(t) + 0.08 * std::cos(2.0 * t);
  }
  return w;
}

// Convolve two polynomials (full convolution)
static std::vector<double> convolve(const std::vector<double>& a,
                                     const std::vector<double>& b) {
  int na = (int)a.size(), nb = (int)b.size();
  std::vector<double> c(na + nb - 1, 0.0);
  for (int i = 0; i < na; i++)
    for (int j = 0; j < nb; j++)
      c[i + j] += a[i] * b[j];
  return c;
}

// Element-wise multiply
static void apply_window(double* out, const double* x, const double* w, int N) {
  for (int n = 0; n < N; n++) out[n] = x[n] * w[n];
}

// ============================================================
// Single-frame GFM-IAIF
// ============================================================

struct GfmResult {
  std::vector<double> av; // nv+1
  std::vector<double> ag; // ng+1
  std::vector<double> al; // 2
};

static GfmResult gfmiaif_frame(const double* s_gvl, int N,
                                int nv, int ng, double d,
                                const double* win) {
  GfmResult res;
  int Lpf = nv + 1; // pre-frame length

  // Pre-frame: ramp from -s_gvl[0] to s_gvl[0]
  std::vector<double> x_gvl(Lpf + N);
  for (int i = 0; i < Lpf; i++)
    x_gvl[i] = -s_gvl[0] + 2.0 * s_gvl[0] * i / (Lpf - 1);
  std::copy(s_gvl, s_gvl + N, x_gvl.begin() + Lpf);
  int total = Lpf + N;

  // Lip radiation coefficients
  res.al = {1.0, -d};

  // Cancel lip radiation: s_gv = filter([1], [1,-d], s_gvl)
  std::vector<double> s_gv = leaky_integrate(s_gvl, N, d);
  std::vector<double> x_gv = leaky_integrate(x_gvl.data(), total, d);

  // Windowed signal for LPC
  std::vector<double> windowed(N);
  apply_window(windowed.data(), s_gv.data(), win, N);

  // --- Gross glottis estimation ---
  // First 1st-order LPC
  std::vector<double> ag1 = lpc(windowed.data(), N, 1);

  // Iterate ng-1 times
  for (int iter = 0; iter < ng - 1; iter++) {
    auto x_v1x = fir_filter(ag1, x_gv.data(), total);
    // Extract after pre-frame
    std::vector<double> s_v1x(x_v1x.begin() + Lpf, x_v1x.end());
    std::vector<double> w_v1x(N);
    apply_window(w_v1x.data(), s_v1x.data(), win, N);
    auto ag1x = lpc(w_v1x.data(), N, 1);
    ag1 = convolve(ag1, ag1x);
  }

  // --- Gross vocal tract estimation ---
  auto x_v1_full = fir_filter(ag1, x_gv.data(), total);
  std::vector<double> s_v1(x_v1_full.begin() + Lpf, x_v1_full.end());
  std::vector<double> w_v1(N);
  apply_window(w_v1.data(), s_v1.data(), win, N);
  auto av1 = lpc(w_v1.data(), N, nv);

  // --- Fine glottis estimation ---
  auto x_g1_full = fir_filter(av1, x_gv.data(), total);
  std::vector<double> s_g1(x_g1_full.begin() + Lpf, x_g1_full.end());
  std::vector<double> w_g1(N);
  apply_window(w_g1.data(), s_g1.data(), win, N);
  res.ag = lpc(w_g1.data(), N, ng);

  // --- Fine vocal tract estimation ---
  auto x_v_full = fir_filter(res.ag, x_gv.data(), total);
  std::vector<double> s_v(x_v_full.begin() + Lpf, x_v_full.end());
  std::vector<double> w_v(N);
  apply_window(w_v.data(), s_v.data(), win, N);
  res.av = lpc(w_v.data(), N, nv);

  return res;
}

} // namespace gfm

// ============================================================
// Frame-based GFM-IAIF processing (exported to R)
// ============================================================

// [[Rcpp::export]]
Rcpp::List gfmiaif_cpp(Rcpp::NumericVector audio, int sample_rate,
                        double window_shift_sec = 0.010,
                        double window_size_sec = 0.032,
                        int nv = 48, int ng = 3, double d = 0.99,
                        std::string window_type = "hann") {

  int N = audio.size();
  int frame_length = (int)(window_size_sec * sample_rate);
  int frame_shift = (int)(window_shift_sec * sample_rate);

  if (frame_length < nv + 2) {
    Rcpp::stop("Frame too short for nv=%d (need >= %d samples)", nv, nv + 2);
  }

  int n_frames = (N - frame_length) / frame_shift + 1;
  if (n_frames < 1) {
    Rcpp::stop("Audio too short: %d samples, need >= %d", N, frame_length);
  }

  // Create window
  std::vector<double> win;
  if (window_type == "hann") win = gfm::hanning(frame_length);
  else if (window_type == "hamming") win = gfm::hamming(frame_length);
  else if (window_type == "blackman") win = gfm::blackman(frame_length);
  else Rcpp::stop("Unknown window type: %s", window_type.c_str());

  // Output matrices
  Rcpp::NumericMatrix av_mat(n_frames, nv + 1);
  Rcpp::NumericMatrix ag_mat(n_frames, ng + 1);
  Rcpp::NumericMatrix al_mat(n_frames, 2);
  Rcpp::NumericVector timestamps(n_frames);

  // Process frames
  for (int i = 0; i < n_frames; i++) {
    int start_idx = i * frame_shift;
    const double* frame_ptr = &audio[start_idx];

    auto res = gfm::gfmiaif_frame(frame_ptr, frame_length, nv, ng, d, win.data());

    // Store results
    for (int j = 0; j <= nv; j++) av_mat(i, j) = res.av[j];
    for (int j = 0; j <= ng; j++) ag_mat(i, j) = res.ag[j];
    for (int j = 0; j < 2; j++) al_mat(i, j) = res.al[j];

    // Center-of-frame timestamp
    timestamps[i] = (start_idx + frame_length / 2.0) / sample_rate;
  }

  return Rcpp::List::create(
    Rcpp::Named("av") = av_mat,
    Rcpp::Named("ag") = ag_mat,
    Rcpp::Named("al") = al_mat,
    Rcpp::Named("timestamps") = timestamps,
    Rcpp::Named("sample_rate_hz") = 1.0 / window_shift_sec,
    Rcpp::Named("n_frames") = n_frames,
    Rcpp::Named("frame_shift_sec") = window_shift_sec,
    Rcpp::Named("frame_size_sec") = window_size_sec
  );
}
