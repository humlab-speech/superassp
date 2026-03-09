// iaif.cpp — IAIF (Iterative Adaptive Inverse Filtering) in C++
//
// Port of inst/python/covarep_python/covarep/glottal/iaif_optimized.py
// and inst/python/covarep_python/covarep/glottal/vq_optimized.py
//
// Reference: P. Alku (1992). "Glottal wave analysis with pitch synchronous
// iterative adaptive inverse filtering". Speech Communication, 11(2-3), 109-118.

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include "dsp_helpers.hpp"

using namespace Rcpp;

// ============================================================
// Core DSP helpers local to iaif.cpp
// (levinson_durbin, fir_filter, leaky_integrate shared via dsp_helpers.hpp)
// ============================================================

// Unbiased autocorrelation for lags 0..order
static std::vector<double> autocorrelation(const double* x, int N, int order) {
  std::vector<double> r(order + 1, 0.0);
  for (int k = 0; k <= order && k < N; k++) {
    double sum = 0.0;
    for (int n = 0; n < N - k; n++) {
      sum += x[n] * x[n + k];
    }
    r[k] = sum;
  }
  return r;
}

// LPC analysis: compute autocorrelation, solve via Levinson-Durbin
static std::vector<double> lpc_analysis(const double* x, int N, int order) {
  std::vector<double> r = autocorrelation(x, N, order);
  if (r[0] == 0.0) {
    std::vector<double> a(order + 1, 0.0);
    a[0] = 1.0;
    return a;
  }
  return levinson_durbin(r, order);
}

// Element-wise multiply x by Hanning window of same length
static void apply_hanning(double* out, const double* x, int N) {
  for (int n = 0; n < N; n++) {
    double w = 0.5 * (1.0 - std::cos(2.0 * M_PI * n / (N - 1)));
    out[n] = x[n] * w;
  }
}

// ============================================================
// High-pass FIR filter design (least-squares, simplified)
// ============================================================

// Design a linear-phase high-pass FIR filter using windowed sinc
// Fstop=40Hz, Fpass=70Hz — we use a simple windowed-sinc approach
static std::vector<double> design_hp_fir(int fs, int Nfir) {
  // Cutoff frequency (midpoint between stop and pass)
  double fc = 55.0 / fs; // normalized cutoff
  std::vector<double> h(Nfir + 1, 0.0);
  int M = Nfir / 2;
  for (int n = 0; n <= Nfir; n++) {
    double nm = n - M;
    if (nm == 0) {
      h[n] = 1.0 - 2.0 * fc;
    } else {
      // High-pass = delta - lowpass
      h[n] = -std::sin(2.0 * M_PI * fc * nm) / (M_PI * nm);
    }
    // Hamming window
    h[n] *= 0.54 - 0.46 * std::cos(2.0 * M_PI * n / Nfir);
  }
  return h;
}

// ============================================================
// IAIF core: 4-step iterative adaptive inverse filtering
// ============================================================

// [[Rcpp::export]]
Rcpp::List iaif_cpp(Rcpp::NumericVector x_in, double fs,
                    int p_vt = -1, int p_gl = -1,
                    double d = 0.99, bool hpfilt = true) {

  int N = x_in.size();
  if (N < 10) {
    return Rcpp::List::create(
      Rcpp::Named("glottal_flow") = Rcpp::NumericVector(0),
      Rcpp::Named("glottal_derivative") = Rcpp::NumericVector(0),
      Rcpp::Named("vt_coeffs") = Rcpp::NumericVector(0),
      Rcpp::Named("gl_coeffs") = Rcpp::NumericVector(0)
    );
  }

  // Default orders (matching MATLAB/Python)
  if (p_vt < 0) p_vt = 2 * (int)std::round(fs / 2000.0) + 4;
  if (p_gl < 0) p_gl = 2 * (int)std::round(fs / 4000.0);
  if (p_gl < 1) p_gl = 1;

  int preflt = p_vt + 1;

  // Copy input
  std::vector<double> x(N);
  for (int i = 0; i < N; i++) x[i] = x_in[i];

  // High-pass filter
  if (hpfilt) {
    int Nfir = (int)std::round(300.0 / 16000.0 * fs);
    if (Nfir % 2 == 1) Nfir++;
    std::vector<double> B = design_hp_fir((int)fs, Nfir);
    int npad = (int)B.size() / 2 - 1;
    // Zero-pad
    std::vector<double> x_padded(N + npad, 0.0);
    std::copy(x.begin(), x.end(), x_padded.begin());
    // Filter
    std::vector<double> x_filt = fir_filter(B, x_padded.data(), (int)x_padded.size());
    // Remove group delay
    int offset = (int)B.size() / 2;
    int new_N = (int)x_filt.size() - offset;
    if (new_N > N) new_N = N;
    if (new_N < preflt + 1) new_N = preflt + 1;
    x.resize(new_N);
    for (int i = 0; i < new_N; i++) {
      x[i] = x_filt[i + offset];
    }
    N = new_N;
  }

  if (N <= p_vt) {
    return Rcpp::List::create(
      Rcpp::Named("glottal_flow") = Rcpp::NumericVector(0),
      Rcpp::Named("glottal_derivative") = Rcpp::NumericVector(0),
      Rcpp::Named("vt_coeffs") = Rcpp::NumericVector(0),
      Rcpp::Named("gl_coeffs") = Rcpp::NumericVector(0)
    );
  }

  // Hanning window
  std::vector<double> win_x(N);
  apply_hanning(win_x.data(), x.data(), N);

  // Pre-frame ramp
  std::vector<double> sig_with_ramp(preflt + N);
  for (int i = 0; i < preflt; i++) {
    sig_with_ramp[i] = -x[0] + (2.0 * x[0]) * i / (preflt - 1);
  }
  std::copy(x.begin(), x.end(), sig_with_ramp.begin() + preflt);
  int total_len = (int)sig_with_ramp.size();

  // === Iteration 1: Estimate glottal+radiation (1st-order LPC) ===
  std::vector<double> Hg1 = lpc_analysis(win_x.data(), N, 1);
  std::vector<double> y_full = fir_filter(Hg1, sig_with_ramp.data(), total_len);
  // Extract after pre-frame
  std::vector<double> y(y_full.begin() + preflt, y_full.end());

  // === Iteration 2: Estimate vocal tract (p_vt order LPC) ===
  std::vector<double> win_y(N);
  apply_hanning(win_y.data(), y.data(), N);
  std::vector<double> Hvt1 = lpc_analysis(win_y.data(), N, p_vt);
  std::vector<double> g1_full = fir_filter(Hvt1, sig_with_ramp.data(), total_len);
  // Leaky integrate
  std::vector<double> g1_int = leaky_integrate(g1_full.data(), total_len, d);
  // Extract after pre-frame
  std::vector<double> g1(g1_int.begin() + preflt, g1_int.end());

  // === Iteration 3: Re-estimate glottal (p_gl order LPC) ===
  std::vector<double> win_g1(N);
  apply_hanning(win_g1.data(), g1.data(), N);
  std::vector<double> Hg2 = lpc_analysis(win_g1.data(), N, p_gl);
  std::vector<double> y2_full = fir_filter(Hg2, sig_with_ramp.data(), total_len);
  std::vector<double> y2_int = leaky_integrate(y2_full.data(), total_len, d);
  std::vector<double> y2(y2_int.begin() + preflt, y2_int.end());

  // === Iteration 4: Final vocal tract estimate ===
  std::vector<double> win_y2(N);
  apply_hanning(win_y2.data(), y2.data(), N);
  std::vector<double> Hvt2 = lpc_analysis(win_y2.data(), N, p_vt);
  std::vector<double> dg_full = fir_filter(Hvt2, sig_with_ramp.data(), total_len);
  // Glottal flow = leaky integrate derivative
  std::vector<double> g_full = leaky_integrate(dg_full.data(), total_len, d);

  // Extract results (remove pre-frame)
  Rcpp::NumericVector glottal_flow(N);
  Rcpp::NumericVector glottal_deriv(N);
  for (int i = 0; i < N; i++) {
    glottal_flow[i] = g_full[preflt + i];
    glottal_deriv[i] = dg_full[preflt + i];
  }

  // Return VT and GL coefficients too
  Rcpp::NumericVector vt_coeffs(Hvt2.begin(), Hvt2.end());
  Rcpp::NumericVector gl_coeffs(Hg2.begin(), Hg2.end());

  return Rcpp::List::create(
    Rcpp::Named("glottal_flow") = glottal_flow,
    Rcpp::Named("glottal_derivative") = glottal_deriv,
    Rcpp::Named("vt_coeffs") = vt_coeffs,
    Rcpp::Named("gl_coeffs") = gl_coeffs
  );
}

// ============================================================
// Voice quality parameter extraction from glottal signals
// ============================================================

// Find peak amplitude in spectrum near target_freq ± bandwidth
static double find_harmonic_peak(const std::vector<double>& spectrum,
                                  const std::vector<double>& freqs,
                                  double target_freq, double bandwidth) {
  double peak = 0.0;
  for (size_t i = 0; i < freqs.size(); i++) {
    if (freqs[i] >= target_freq - bandwidth &&
        freqs[i] <= target_freq + bandwidth) {
      if (spectrum[i] > peak) peak = spectrum[i];
    }
  }
  return peak;
}

// In-place radix-2 Cooley-Tukey FFT (N must be power of 2)
static void fft_inplace(std::vector<double>& re, std::vector<double>& im, int N) {
  // Bit-reversal permutation
  for (int i = 1, j = 0; i < N; i++) {
    int bit = N >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) { std::swap(re[i], re[j]); std::swap(im[i], im[j]); }
  }
  // Butterfly
  for (int len = 2; len <= N; len <<= 1) {
    double ang = -2.0 * M_PI / len;
    double wre = std::cos(ang), wim = std::sin(ang);
    for (int i = 0; i < N; i += len) {
      double ure = 1.0, uim = 0.0;
      for (int j = 0; j < len / 2; j++) {
        double tre = re[i + j + len/2] * ure - im[i + j + len/2] * uim;
        double tim = re[i + j + len/2] * uim + im[i + j + len/2] * ure;
        re[i + j + len/2] = re[i + j] - tre;
        im[i + j + len/2] = im[i + j] - tim;
        re[i + j] += tre;
        im[i + j] += tim;
        double new_ure = ure * wre - uim * wim;
        uim = ure * wim + uim * wre;
        ure = new_ure;
      }
    }
  }
}

// Magnitude spectrum via FFT — uses central window of max 4096 samples
static void compute_spectrum(const double* x, int N,
                              std::vector<double>& mag,
                              std::vector<double>& freqs, double fs) {
  // Limit to a manageable window (center of signal)
  int max_win = 4096;
  const double* data = x;
  int data_len = N;
  std::vector<double> windowed;
  if (N > max_win) {
    int start = (N - max_win) / 2;
    data_len = max_win;
    windowed.resize(data_len);
    for (int i = 0; i < data_len; i++) {
      double w = 0.5 * (1.0 - std::cos(2.0 * M_PI * i / (data_len - 1)));
      windowed[i] = x[start + i] * w;
    }
    data = windowed.data();
  }

  // Next power of 2
  int nfft = 1;
  while (nfft < data_len) nfft <<= 1;
  int n_bins = nfft / 2 + 1;

  std::vector<double> re(nfft, 0.0), im(nfft, 0.0);
  for (int i = 0; i < data_len; i++) re[i] = data[i];

  fft_inplace(re, im, nfft);

  mag.resize(n_bins);
  freqs.resize(n_bins);
  for (int k = 0; k < n_bins; k++) {
    freqs[k] = (double)k * fs / nfft;
    mag[k] = std::sqrt(re[k] * re[k] + im[k] * im[k]);
  }
}

// [[Rcpp::export]]
Rcpp::List extract_vq_params_cpp(Rcpp::NumericVector glottal_flow,
                                  Rcpp::NumericVector glottal_derivative,
                                  double fs,
                                  double f0 = -1.0,
                                  Rcpp::Nullable<Rcpp::IntegerVector> gci = R_NilValue) {

  int N = glottal_flow.size();
  int Nd = glottal_derivative.size();

  // Basic amplitude measures
  double flow_max = *std::max_element(glottal_flow.begin(), glottal_flow.end());
  double flow_min = *std::min_element(glottal_flow.begin(), glottal_flow.end());
  double deriv_peak = 0.0;
  for (int i = 0; i < Nd; i++) {
    double v = std::abs(glottal_derivative[i]);
    if (v > deriv_peak) deriv_peak = v;
  }

  // NAQ and QOQ
  double naq = NA_REAL;
  double qoq = NA_REAL;

  if (gci.isNotNull()) {
    Rcpp::IntegerVector gci_v(gci);
    int n_gci = gci_v.size();
    if (n_gci >= 2) {
      std::vector<double> naq_vals, qoq_vals;
      for (int i = 0; i < n_gci - 1; i++) {
        int start = gci_v[i];
        int end = gci_v[i + 1];
        if (end - start < 10 || end >= N) continue;

        // NAQ
        double d_peak_period = 0.0;
        for (int j = start; j < end && j < Nd; j++) {
          double v = std::abs(glottal_derivative[j]);
          if (v > d_peak_period) d_peak_period = v;
        }
        double f_max = glottal_flow[start], f_min = glottal_flow[start];
        for (int j = start; j < end; j++) {
          if (glottal_flow[j] > f_max) f_max = glottal_flow[j];
          if (glottal_flow[j] < f_min) f_min = glottal_flow[j];
        }
        double f_ac = f_max - f_min;
        if (f_ac > 0.0) {
          double t_period = (double)(end - start) / fs;
          naq_vals.push_back(d_peak_period / (f_ac / t_period));
        }

        // QOQ
        int peak_idx = 0;
        double peak_val = 0.0;
        int period_len = end - start;
        for (int j = 0; j < period_len && start + j < Nd; j++) {
          double v = std::abs(glottal_derivative[start + j]);
          if (v > peak_val) { peak_val = v; peak_idx = j; }
        }
        double open_dur = (double)(period_len - peak_idx) / period_len;
        qoq_vals.push_back(open_dur);
      }
      // Median
      if (!naq_vals.empty()) {
        std::sort(naq_vals.begin(), naq_vals.end());
        naq = naq_vals[naq_vals.size() / 2];
      }
      if (!qoq_vals.empty()) {
        std::sort(qoq_vals.begin(), qoq_vals.end());
        qoq = qoq_vals[qoq_vals.size() / 2];
      }
    }
  }

  // Compute spectrum once; reused for H1-H2, HRF, and PSP
  std::vector<double> spectrum, freqs;
  compute_spectrum(glottal_derivative.begin(), Nd, spectrum, freqs, fs);

  // H1-H2
  double h1_h2 = NA_REAL;
  if (f0 > 0.0 && f0 <= fs / 4.0) {
    double h1 = find_harmonic_peak(spectrum, freqs, f0, 50.0);
    double h2 = find_harmonic_peak(spectrum, freqs, 2.0 * f0, 50.0);
    if (h1 > 0.0 && h2 > 0.0) {
      h1_h2 = 20.0 * std::log10(h1 / h2);
    }
  }

  // HRF (Harmonic Richness Factor)
  double hrf = NA_REAL;
  {
    double low_energy = 0.0, high_energy = 0.0;
    for (size_t i = 0; i < freqs.size(); i++) {
      double power = spectrum[i] * spectrum[i];
      if (freqs[i] < 2000.0) low_energy += power;
      else high_energy += power;
    }
    if (low_energy > 0.0) {
      hrf = 10.0 * std::log10((high_energy + 1e-10) / low_energy);
    }
  }

  // PSP (Parabolic Spectral Parameter)
  double psp = NA_REAL;
  {
    // Fit parabola to log spectrum up to 5kHz
    std::vector<double> xx, yy;
    for (size_t i = 0; i < freqs.size() && freqs[i] <= 5000.0; i++) {
      xx.push_back(freqs[i] / 1000.0);
      yy.push_back(20.0 * std::log10(spectrum[i] + 1e-10));
    }
    int n_pts = (int)xx.size();
    if (n_pts >= 3) {
      // Least squares: y = a*x^2 + b*x + c
      // Normal equations: X'X * [a,b,c]' = X'y
      double S0 = n_pts, S1 = 0, S2 = 0, S3 = 0, S4 = 0;
      double Sy = 0, Sxy = 0, Sx2y = 0;
      for (int i = 0; i < n_pts; i++) {
        double xi = xx[i], xi2 = xi * xi, yi = yy[i];
        S1 += xi; S2 += xi2; S3 += xi2 * xi; S4 += xi2 * xi2;
        Sy += yi; Sxy += xi * yi; Sx2y += xi2 * yi;
      }
      // Solve 3x3 system via Cramer's rule
      // [S4 S3 S2] [a]   [Sx2y]
      // [S3 S2 S1] [b] = [Sxy]
      // [S2 S1 S0] [c]   [Sy]
      double det = S4*(S2*S0 - S1*S1) - S3*(S3*S0 - S1*S2) + S2*(S3*S1 - S2*S2);
      if (std::abs(det) > 1e-15) {
        double a_num = Sx2y*(S2*S0 - S1*S1) - S3*(Sxy*S0 - S1*Sy) + S2*(Sxy*S1 - S2*Sy);
        psp = a_num / det;
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("glottal_flow_max") = flow_max,
    Rcpp::Named("glottal_flow_min") = flow_min,
    Rcpp::Named("glottal_derivative_peak") = deriv_peak,
    Rcpp::Named("NAQ") = naq,
    Rcpp::Named("QOQ") = qoq,
    Rcpp::Named("H1_H2") = h1_h2,
    Rcpp::Named("HRF") = hrf,
    Rcpp::Named("PSP") = psp
  );
}
