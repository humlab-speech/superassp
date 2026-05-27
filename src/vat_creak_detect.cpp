// creak_detect.cpp — Bit-faithful port of the Kane-Drugman + Ishi creak feature
// extraction pipeline (Kane, Drugman, Gobl, "Improved automatic detection of
// creak", CSL 2013) feeding the bundled ANN classifier (vat_ann.cpp).
//
// Mirrors:
//   creak_fcns/get_ALL_creak_features.m
//   creak_fcns/get_KD_creak_features.m
//   creak_fcns/sil_unv_features.m
//   creak_fcns/get_short_pow.m
//   creak_fcns/get_creak_H2H1.m
//   creak_fcns/res_peak.m
//   creak_fcns/get_res_peak_prom.m
//   creak_fcns/get_ishi_params_inter.m
//   creak_fcns/getIFP.m
//   creak_fcns/ishi_creak_detection_ORIG.m
//   creak_fcns/getZeroXrate.m
//   creak_fcns/get_delta_mat.m

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "vat_dsp.h"
#include "vat_lpc.h"

using namespace Rcpp;

// ---------------- Common helpers ----------------

// MATLAB xcorr(a, b): returns length 2N-1 vector for equal-length inputs.
// out(k) = sum_n a(n) * b(n - (k - (N-1))) for k = 0..2N-2.
static arma::vec xcorr_full(const arma::vec& a, const arma::vec& b) {
  int N = a.n_elem;
  arma::vec out(2 * N - 1, arma::fill::zeros);
  for (int k = 0; k < 2 * N - 1; ++k) {
    int lag = k - (N - 1);
    int n_lo = std::max(0, lag);
    int n_hi = std::min(N - 1, N - 1 + lag);
    double s = 0.0;
    for (int n = n_lo; n <= n_hi; ++n) s += a(n) * b(n - lag);
    out(k) = s;
  }
  return out;
}

// Equivalent to MATLAB xcorr(x, maxlag, 'coeff'): normalized to lag 0 = 1,
// only returns lags in [-maxlag, +maxlag] (length 2*maxlag+1).
static void xcorr_coeff(const arma::vec& x, int maxlag,
                         arma::vec& acf, arma::vec& lags) {
  int N = x.n_elem;
  double norm0 = arma::dot(x, x);
  if (norm0 <= 0) norm0 = 1.0;
  int L = 2 * maxlag + 1;
  acf.zeros(L);
  lags.set_size(L);
  for (int k = 0; k < L; ++k) {
    int lag = k - maxlag;
    lags(k) = lag;
    int n_lo = std::max(0, lag);
    int n_hi = std::min(N - 1, N - 1 + lag);
    double s = 0.0;
    for (int n = n_lo; n <= n_hi; ++n) s += x(n) * x(n - lag);
    acf(k) = s / norm0;
  }
}

// Moving average smooth (MATLAB smooth(x, n), centered, edges shorten span).
static arma::vec smooth_ma(const arma::vec& x, int span) {
  if (span <= 1) return x;
  int N = x.n_elem;
  arma::vec out(N);
  int half = span / 2;
  for (int i = 0; i < N; ++i) {
    int lo = std::max(0, i - half);
    int hi = std::min(N - 1, i + half);
    int w = std::min(i, N - 1 - i);
    w = std::min(w, half);
    lo = i - w; hi = i + w;
    out(i) = arma::mean(x.subvec(lo, hi));
  }
  return out;
}

// Resample-to-length via linear interpolation.
static arma::vec interp_to_length(const arma::vec& x, int N_out) {
  if ((int)x.n_elem == N_out) return x;
  if (x.n_elem < 2) {
    arma::vec out(N_out, arma::fill::zeros);
    if (x.n_elem == 1) out.fill(x(0));
    return out;
  }
  arma::vec xs = arma::linspace<arma::vec>(0.0, double(N_out - 1), x.n_elem);
  arma::vec xq = arma::regspace<arma::vec>(0, N_out - 1);
  return vat::interp1_linear(xs, x, xq);
}

// ---------------- sil_unv_features ----------------
// Returns Es (log mean-energy per frame after Hanning) and ZCs_ms (zero-cross
// count, scaled like MATLAB code), plus the time-grid pos (centers in samples).
struct SilUnvResult {
  arma::vec Es;       // log energy (frame_len_ms / shift 5 ms)
  arma::vec ZCs_ms;   // zero-cross count / 1000 * fs
  arma::vec pos;      // frame center sample index (1-based equivalent in MATLAB)
  arma::vec Xpos;     // zero-cross frame center positions (for ZCs)
};
static SilUnvResult sil_unv_features(const arma::vec& wave_in, double fs,
                                      double frame_len_ms = 10.0) {
  arma::vec wave = wave_in * std::pow(2.0, 15);
  int frame_len = (int)std::round(frame_len_ms / 1000.0 * fs);
  int shift     = (int)std::round(5.0 / 1000.0 * fs);
  int N = wave.n_elem;
  if (frame_len < 2) frame_len = 2;

  // Zero crossing rate over the same frame_len/shift grid.
  arma::vec zc, xpos;
  {
    std::vector<double> v_zc, v_xpos;
    int start = 0, stop = frame_len - 1;
    while (stop < N) {
      // Count sign changes
      int xings = 0;
      for (int n = start + 1; n <= stop; ++n)
        if ((wave(n) >= 0) != (wave(n - 1) >= 0)) ++xings;
      v_zc.push_back(double(xings) / double(frame_len));
      v_xpos.push_back(0.5 * (start + stop + 1));  // 1-based midpoint
      start += shift;
      stop = start + frame_len - 1;
    }
    zc.set_size(v_zc.size());
    xpos.set_size(v_xpos.size());
    for (size_t i = 0; i < v_zc.size(); ++i) { zc(i) = v_zc[i]; xpos(i) = v_xpos[i]; }
  }
  arma::vec zc_scaled = zc / 1000.0 * fs;  // matches MATLAB ZCs_ms = ZCs/1000*Fs

  // Frame-wise log-energy.
  arma::vec hwin = vat::hanning(frame_len);
  std::vector<double> v_E, v_pos;
  int start = 0, stop = frame_len - 1;
  while (stop < N) {
    arma::vec seg = wave.subvec(start, stop) % hwin;
    v_E.push_back(arma::mean(seg % seg));
    start += shift;
    stop = start + frame_len - 1;
    v_pos.push_back(0.5 * (start + stop + 1));
  }
  SilUnvResult r;
  r.Es.set_size(v_E.size());
  r.pos.set_size(v_pos.size());
  for (size_t i = 0; i < v_E.size(); ++i) {
    r.Es(i)  = std::log(v_E[i] + 1e-30);
    r.pos(i) = v_pos[i];
  }
  r.ZCs_ms = zc_scaled;
  r.Xpos   = xpos;
  return r;
}

// ---------------- get_short_pow ----------------
// Returns (pow_dB length n_short, pow_std length n_short, pow_std_inter length N).
struct ShortPowResult {
  arma::vec pow_dB;
  arma::vec pow_std;
  arma::vec pow_std_inter;
  arma::vec t_pow;  // 1-based MATLAB-style center positions in samples
};
static ShortPowResult get_short_pow(const arma::vec& x, double fs) {
  int veryShort_len   = (int)std::round(4.0 * fs / 1000.0);
  int veryShort_shift = (int)std::round(2.0 * fs / 1000.0);
  int N = x.n_elem;
  std::vector<double> vp, vt;
  int start = 0, finish = veryShort_len - 1;
  while (finish < N) {
    double s = 0.0;
    for (int k = start; k <= finish; ++k) s += x(k) * x(k);
    vp.push_back(s / double(veryShort_len));
    vt.push_back(0.5 * (start + finish + 1) + 1.0);  // 1-based center
    start += veryShort_shift;
    finish = start + veryShort_len - 1;
  }
  int M = vp.size();
  arma::vec pow_dB(M), t_pow(M);
  for (int i = 0; i < M; ++i) {
    pow_dB(i) = 20.0 * std::log10(vp[i] + 1e-30);
    t_pow(i) = vt[i];
  }
  // Replace -Inf with min of finite values (MATLAB Kane quirk)
  std::vector<double> finite_vals;
  for (int i = 0; i < M; ++i)
    if (std::isfinite(pow_dB(i))) finite_vals.push_back(pow_dB(i));
  double minfin = finite_vals.empty() ? -200.0 :
    *std::min_element(finite_vals.begin(), finite_vals.end());
  for (int i = 0; i < M; ++i)
    if (!std::isfinite(pow_dB(i))) pow_dB(i) = minfin;

  arma::vec pow_std(M, arma::fill::zeros);
  int std_len = 16;
  for (int n = std_len; n < M - std_len; ++n)
    pow_std(n) = arma::stddev(pow_dB.subvec(n - std_len, n + std_len));
  pow_std = vat::medfilt1(pow_std, 13);

  ShortPowResult r;
  r.pow_dB = pow_dB;
  r.pow_std = pow_std;
  r.t_pow = t_pow;
  r.pow_std_inter = interp_to_length(pow_std, N);
  return r;
}

// ---------------- get_creak_H2H1 ----------------
struct H2H1Result { arma::vec H2H1; arma::vec creakF0; };
static H2H1Result get_creak_h2h1(const arma::vec& res, double fs, double F0mean) {
  int N = res.n_elem;
  int Stop_init = (int)std::round(50.0 / 1000.0 * fs);
  int Hop = (int)std::round(10.0 / 1000.0 * fs);
  arma::vec win = vat::hanning(Stop_init);
  int moving_average_size = (int)std::round(100.0 / 1000.0 * fs);
  int moving_average_frame = std::max(1, moving_average_size / Hop);

  // Resonator outputs.
  auto reson = [&](double Rho) {
    double Phi = 2.0 * arma::datum::pi * F0mean / fs;
    arma::vec b = {1.0, 0.0, 0.0};
    arma::vec a = {1.0, -2.0 * Rho * std::cos(Phi), Rho * Rho};
    return vat::filter(b, a, res);
  };
  arma::vec rep  = reson(0.97);
  arma::vec rep2 = reson(0.80);
  double mx = arma::max(arma::abs(rep));
  if (mx > 0.0) rep /= mx;

  std::vector<double> f0_vec, delta_vec;
  int start = 0, stop = Stop_init - 1;
  int nfft = (int)fs;
  while (stop < (int)rep.n_elem) {
    // F0 estimate from rep2
    arma::vec Sig = rep2.subvec(start, stop) % win;
    arma::vec corrs = xcorr_full(Sig, Sig);
    arma::vec C1 = corrs.subvec(Sig.n_elem, corrs.n_elem - 1);  // drop first len(Sig)
    int Lc = C1.n_elem;
    // Unbiased correction
    for (int k = 0; k < Lc; ++k)
      C1(k) *= double(Lc) / double(Lc - k);
    arma::vec C2 = C1;
    int k = 0;
    while (k < (int)C2.n_elem && C2(k) > 0) { C2(k) = 0; ++k; }
    arma::uword posi = 0;
    C2.max(posi);
    int posi_i = (int)posi;
    if (posi_i < 5) posi_i = 5;
    int F0 = (int)std::round(fs / double(posi_i));
    f0_vec.push_back(F0);

    // H2-H1 from rep
    arma::vec Sig2 = rep.subvec(start, stop) % win;
    arma::vec corrs_r = xcorr_full(Sig2, Sig2);
    arma::vec Cr = corrs_r.subvec(Sig2.n_elem, corrs_r.n_elem - 1);

    arma::cx_vec Spec_c = vat::fft(Cr, (arma::uword)nfft);
    int half = nfft / 2;
    arma::vec Spec(half);
    for (int i = 0; i < half; ++i) Spec(i) = std::abs(Spec_c(i));
    double nrm = std::sqrt(arma::dot(Spec, Spec));
    if (nrm > 0.0) Spec /= nrm;
    for (int i = 0; i < half; ++i) Spec(i) = 20.0 * std::log10(Spec(i) + 1e-30);

    double delta = 0.0;
    if (2 * F0 < half && F0 < half && F0 >= 1)
      delta = Spec(2 * F0) - Spec(F0);
    delta_vec.push_back(delta);

    start += Hop;
    stop = start + Stop_init - 1;
  }

  H2H1Result out;
  if (delta_vec.empty()) {
    out.H2H1 = arma::zeros<arma::vec>(N);
    out.creakF0 = arma::zeros<arma::vec>(N);
    return out;
  }
  arma::vec delta(delta_vec.size()), f0v(f0_vec.size());
  for (size_t i = 0; i < delta_vec.size(); ++i) { delta(i) = delta_vec[i]; f0v(i) = f0_vec[i]; }
  delta = smooth_ma(delta, moving_average_frame);
  f0v   = vat::medfilt1(f0v, 7);
  f0v   = smooth_ma(f0v, moving_average_frame);

  out.H2H1    = interp_to_length(delta, N);
  out.creakF0 = interp_to_length(f0v, N);
  return out;
}

// ---------------- res_peak / get_res_peak_prom ----------------
// res_peak returns per-sample interpolated peak_inter contour.
static arma::vec res_peak(const arma::vec& x, double fs, double F0mean,
                            const arma::vec& res, const arma::vec& Es) {
  int N = x.n_elem;
  double maxCreakF0;
  if (F0mean > 100)       maxCreakF0 = 90;
  else if (F0mean >= 85)  maxCreakF0 = 65;
  else if (F0mean < 85)   maxCreakF0 = 55;
  else                    maxCreakF0 = 80;

  int winLen = (int)std::round(fs / maxCreakF0) * 2;
  double Phi = 2.0 * arma::datum::pi * F0mean / fs;
  double Rho = 0.8;
  arma::vec b = {1.0, 0.0, 0.0};
  arma::vec a = {1.0, -2.0 * Rho * std::cos(Phi), Rho * Rho};
  arma::vec rep = vat::filter(b, a, res);

  // get_res_peak_prom
  arma::vec Es_inter = interp_to_length(Es, N);
  int ener_tol = (int)std::round(400.0 / 1000.0 * fs);
  double ener_thresh = 4.0;
  int winLenRem = (int)std::round(winLen * 0.2);

  int shift = (int)std::round(20.0 / 1000.0 * fs);
  int Mrep = rep.n_elem;
  std::vector<double> v_peak_prom, v_peak_t;
  int start = 0, stop = winLen - 1;
  while (stop < Mrep) {
    double pt = 0.5 * (start + stop + 1) + 1.0;
    arma::vec frame_init = rep.subvec(start, stop);
    arma::uword maxIdx = arma::abs(frame_init).index_max();
    int gpos = (int)maxIdx + start;
    double prom = 0.0;
    if (gpos - winLen / 2 >= 0 && gpos + winLen / 2 < Mrep) {
      int start2 = gpos - winLen / 2;
      int stop2  = gpos + winLen / 2;
      arma::vec frame2 = rep.subvec(start2, stop2);
      arma::uword m2 = arma::abs(frame2).index_max();
      if (frame2(m2) < 0) frame2 *= -1;
      double fmax = frame2.max();
      if (fmax > 0) frame2 /= fmax;
      int midpoint = winLen / 2;
      int lo = std::max(0, midpoint - winLenRem);
      int hi = std::min((int)frame2.n_elem - 1, midpoint + winLenRem);
      frame2.subvec(lo, hi).zeros();
      frame2.transform([](double v) { return std::max(0.0, v); });
      prom = 1.0 - frame2.max();

      int se = std::max(0, gpos - ener_tol);
      int ee = std::min(N - 1, gpos + ener_tol);
      double max_e = (se <= ee) ? Es_inter.subvec(se, ee).max() : 0.0;
      if (gpos >= 0 && gpos < N && Es_inter(gpos) < max_e - ener_thresh)
        prom = 0.0;
    }
    v_peak_prom.push_back(prom);
    v_peak_t.push_back(pt);
    start += shift;
    stop = start + winLen - 1;
  }

  if (v_peak_t.size() < 2) return arma::zeros<arma::vec>(N);
  arma::vec pp(v_peak_prom.size()), pt(v_peak_t.size());
  for (size_t i = 0; i < v_peak_prom.size(); ++i) { pp(i) = v_peak_prom[i]; pt(i) = v_peak_t[i]; }
  pp = vat::medfilt1(pp, 5);

  // Interp1 onto 1..N samples
  arma::vec xq = arma::regspace<arma::vec>(1, N);
  arma::vec out = vat::interp1_linear(pt, pp, xq);
  out.replace(arma::datum::nan, 0.0);
  return out;
}

// ---------------- getIFP ----------------
// Returns (IFP, t_IFP); applies the successive-frames silencing post-filter.
static void get_ifp(const arma::vec& x_filt, double fs, double IFPthresh,
                     arma::vec& IFP, arma::vec& t_IFP) {
  int N = (int)std::round(32.0 * fs / 1000.0);
  int shift = (int)std::round(10.0 * fs / 1000.0);
  int maxLag = (int)std::round(15.0 * fs / 1000.0);
  int peak_tol = (int)std::round(1.0 * fs / 1000.0);
  int Nx = x_filt.n_elem;
  std::vector<double> v_IFP, v_t;
  int start = 0, finish = N - 1;
  while (finish < Nx) {
    arma::vec frame = x_filt.subvec(start, finish);
    arma::vec acf, lags;
    xcorr_coeff(frame, maxLag, acf, lags);
    // Keep lags > 0 only
    arma::vec acf_pos = acf.subvec(maxLag + 1, acf.n_elem - 1);  // lags 1..maxLag
    // Zero first peak_tol bins
    int zero_n = std::min((int)acf_pos.n_elem, peak_tol);
    if (zero_n > 0) acf_pos.subvec(0, zero_n - 1).zeros();
    // Unbiased normalisation factor
    int La = acf_pos.n_elem;
    for (int k = 0; k < La; ++k)
      acf_pos(k) *= double(N) / double(N - (k + 1));
    // Pick harmonic peaks at m * max_idx
    arma::uword imax;
    double acf_max1 = acf_pos.max(imax);
    std::vector<double> acf_maxes; acf_maxes.push_back(acf_max1);
    int m1 = (int)imax + 1;  // 1-indexed lag
    int mm = 2;
    while (m1 * mm < La) {
      int begin = std::max(0, m1 * mm - peak_tol - 1);
      int stop2 = std::min(La - 1, m1 * mm + peak_tol - 1);
      if (stop2 < begin) break;
      double mx = acf_pos.subvec(begin, stop2).max();
      acf_maxes.push_back(mx);
      ++mm;
    }
    double v = acf_maxes.empty() ? 0.0 :
      *std::min_element(acf_maxes.begin(), acf_maxes.end());
    v_IFP.push_back(v);
    v_t.push_back(0.5 * (start + finish + 1) + 1.0);
    start += shift;
    finish = start + N - 1;
  }
  int M = v_IFP.size();
  IFP.set_size(M); t_IFP.set_size(M);
  for (int i = 0; i < M; ++i) { IFP(i) = v_IFP[i]; t_IFP(i) = v_t[i]; }
  // Successive-frame silencing
  int succ = 3;
  for (int n = 0; n + succ < M; ++n)
    if (IFP(n + 1) < IFPthresh || IFP(n + 2) < IFPthresh) IFP(n) = 0;
  // Clamp to [0,1]
  IFP.transform([](double v) { return std::max(0.0, std::min(1.0, v)); });
}

// ---------------- ishi feature extraction ----------------
struct IshiResult {
  arma::vec IFP_inter;        // length N (per-sample)
  arma::vec IPS_inter;        // length N
  arma::vec PwP_rise_inter;   // length N
  arma::vec PwP_fall_inter;   // length N
};
static IshiResult get_ishi_params(const arma::vec& x_in, double fs) {
  int N = x_in.n_elem;
  IshiResult r;
  r.IFP_inter      = arma::zeros<arma::vec>(N);
  r.IPS_inter      = arma::zeros<arma::vec>(N);
  r.PwP_rise_inter = arma::zeros<arma::vec>(N);
  r.PwP_fall_inter = arma::zeros<arma::vec>(N);

  // Scale + 100-1500 Hz bandpass (Butterworth order 4)
  arma::vec x = x_in * std::pow(2.0, 15);
  arma::vec b, a;
  arma::vec Wn = {100.0 / (fs / 2.0), 1500.0 / (fs / 2.0)};
  vat::butter(4, Wn, "bandpass", b, a);
  arma::vec x_filt = vat::filtfilt(b, a, x);

  // Very short-term power contour
  int vs_len = (int)std::round(4.0 * fs / 1000.0);
  int vs_shift = (int)std::round(2.0 * fs / 1000.0);
  std::vector<double> v_pow, v_tpow;
  int start = 0, finish = vs_len - 1;
  while (finish < (int)x_filt.n_elem) {
    double s = 0.0;
    for (int k = start; k <= finish; ++k) s += x_filt(k) * x_filt(k);
    v_pow.push_back(s / double(vs_len));
    v_tpow.push_back(0.5 * (start + finish + 1) + 1.0);
    start += vs_shift;
    finish = start + vs_len - 1;
  }
  int M = v_pow.size();
  if (M < 5) return r;
  arma::vec pow_dB(M), t_pow(M);
  for (int i = 0; i < M; ++i) {
    pow_dB(i) = 20.0 * std::log10(v_pow[i] + 1e-30);
    t_pow(i) = v_tpow[i];
  }

  // Peaks of pow_dB
  vat::PeakResult pr = vat::findpeaks(pow_dB);
  if (pr.locs.n_elem < 3) return r;

  // Local-peak filtering (drop peaks within 2 dB of neighbours 3 bins out)
  double localpeak_thresh = 2.0;
  std::vector<arma::uword> kept_loc;
  std::vector<double> kept_pk;
  for (arma::uword i = 0; i < pr.locs.n_elem; ++i) {
    int p = pr.locs(i);
    int st = std::max(0, p - 3);
    int fi = std::min(M - 1, p + 3);
    double pkv = pr.pks(i);
    if (pkv - localpeak_thresh >= pow_dB(st) &&
        pkv - localpeak_thresh >= pow_dB(fi)) {
      kept_loc.push_back(p);
      kept_pk.push_back(pkv);
    }
  }
  // Drop peaks > 30 dB below max
  double maxPow = kept_pk.empty() ? 0.0 :
    *std::max_element(kept_pk.begin(), kept_pk.end());
  double max2remove_dB = 30.0;
  std::vector<arma::uword> kept_loc2;
  std::vector<double> kept_pk2;
  for (size_t i = 0; i < kept_pk.size(); ++i)
    if (kept_pk[i] >= maxPow - max2remove_dB) {
      kept_loc2.push_back(kept_loc[i]);
      kept_pk2.push_back(kept_pk[i]);
    }
  int Np = kept_loc2.size();
  if (Np < 3) return r;

  // PwP rise / fall around each peak
  arma::vec PwP_rise(Np, arma::fill::zeros);
  arma::vec PwP_fall(Np, arma::fill::zeros);
  for (int n = 0; n < Np; ++n) {
    int p = kept_loc2[n];
    int st = std::max(0, p - 5);
    int fi = std::min(M - 1, p + 5);
    if (p > 0 && st <= p - 1) {
      double mx = pow_dB(st);
      for (int k = st; k <= p - 1; ++k) mx = std::min(mx, pow_dB(k));
      PwP_rise(n) = kept_pk2[n] - mx;
    }
    if (p < M - 1 && p + 1 <= fi) {
      double mx = pow_dB(p + 1);
      for (int k = p + 1; k <= fi; ++k) mx = std::min(mx, pow_dB(k));
      PwP_fall(n) = kept_pk2[n] - mx;
    }
  }

  // IFP (frame-based)
  arma::vec IFP, t_IFP;
  get_ifp(x_filt, fs, 0.5, IFP, t_IFP);

  // IPS (peak-synchronous)
  arma::vec IPS(Np, arma::fill::zeros);
  int N_ips = (int)std::round(15.0 * fs / 1000.0);
  int maxLags = N_ips;
  int Tmax = (int)std::round(100.0 * fs / 1000.0);
  // Convert peak locs to sample indices
  std::vector<int> peak_samp(Np);
  for (int n = 0; n < Np; ++n) peak_samp[n] = (int)std::round(t_pow(kept_loc2[n]));
  int Nx = x_filt.n_elem;
  for (int n = 1; n < Np; ++n) {
    if (peak_samp[n] - peak_samp[n - 1] > Tmax) { IPS(n) = 0; continue; }
    int half = N_ips / 2;
    int s1 = peak_samp[n] - half, e1 = peak_samp[n] + half;
    int s2 = peak_samp[n - 1] - half, e2 = peak_samp[n - 1] + half;
    if (s1 < 0 || e1 >= Nx || s2 < 0 || e2 >= Nx) { IPS(n) = 0; continue; }
    arma::vec a1 = x_filt.subvec(s1, e1);
    arma::vec a2 = x_filt.subvec(s2, e2);
    arma::vec acf, lags;
    xcorr_coeff(a1 - a2 + a1, maxLags, acf, lags);  // placeholder — replace
    // Proper: cross-correlation between a1, a2, normalized to lag 0 of each
    // Implemented manually:
    int Lh = 2 * maxLags + 1;
    arma::vec cc(Lh);
    double n1 = std::sqrt(arma::dot(a1, a1) * arma::dot(a2, a2));
    if (n1 <= 0) { IPS(n) = 0; continue; }
    int La = a1.n_elem;
    for (int k = 0; k < Lh; ++k) {
      int lag = k - maxLags;
      int n_lo = std::max(0, lag);
      int n_hi = std::min(La - 1, La - 1 + lag);
      double s = 0.0;
      for (int nn = n_lo; nn <= n_hi; ++nn) s += a1(nn) * a2(nn - lag);
      cc(k) = s / n1;
    }
    IPS(n) = cc.max();
  }

  // Interp Ishi outputs to per-sample
  if (t_IFP.n_elem >= 2 && IFP.n_elem >= 2) {
    arma::vec xq = arma::regspace<arma::vec>(1, N);
    arma::vec ifp_i = vat::interp1_linear(t_IFP, IFP, xq);
    ifp_i.replace(arma::datum::nan, 0.0);
    r.IFP_inter = ifp_i;
  }

  // PwP / IPS scattering via nearest-peak rule (35 ms tolerance)
  int maxLen = (int)std::round(35.0 / 1000.0 * fs);
  for (int n = 0; n < N; ++n) {
    int best_idx = -1, best_dist = maxLen + 1;
    for (int p = 0; p < Np; ++p) {
      if (peak_samp[p] < n) continue;
      int d = peak_samp[p] - n;
      if (d < best_dist) { best_dist = d; best_idx = p; }
    }
    if (best_idx >= 0 && best_dist < maxLen) {
      r.PwP_rise_inter(n) = PwP_rise(best_idx);
      r.PwP_fall_inter(n) = PwP_fall(best_idx);
      r.IPS_inter(n)      = IPS(best_idx);
    }
  }
  return r;
}

// ---------------- get_delta_mat ----------------
// MATLAB: delta_mat(2:end, n) = diff(mat(:, n)); first row = 0.
static arma::mat delta_mat(const arma::mat& m) {
  arma::mat d(m.n_rows, m.n_cols, arma::fill::zeros);
  if (m.n_rows < 2) return d;
  for (arma::uword n = 0; n < m.n_cols; ++n)
    for (arma::uword i = 1; i < m.n_rows; ++i)
      d(i, n) = m(i, n) - m(i - 1, n);
  return d;
}

// ---------------- get_KD_creak_features (top-level KD pipeline) ----------------
// Builds the 8-feature matrix used by get_ALL_creak_features:
// [H2H1 res_p ZCR F0 F0mean enerN pow_std creakF0] resampled to 10 ms.

// LP residual via Drugman's GetLPCresidual (mirrors the one in gci_se_vq.cpp).
static arma::vec lpc_residual_drugman(const arma::vec& wave, int L, int shift, int order) {
  int N = wave.n_elem;
  arma::vec res(N, arma::fill::zeros);
  arma::vec hann_win = vat::hanning(L + 1);
  int start = 0, stop = start + L;
  while (stop < N) {
    arma::vec seg = wave.subvec(start, stop) % hann_win;
    arma::vec ar; double e;
    arma::vec r(order + 1, arma::fill::zeros);
    int nseg = seg.n_elem;
    for (int k = 0; k <= order; ++k)
      for (int i = 0; i < nseg - k; ++i) r(k) += seg(i) * seg(i + k);
    vat::levinson(r, order, ar, e);
    arma::vec inv = vat::filter(ar, arma::vec({1.0}), seg);
    double e_seg = arma::dot(seg, seg);
    double e_inv = arma::dot(inv, inv);
    if (e_inv > 0.0) inv *= std::sqrt(e_seg / e_inv);
    res.subvec(start, stop) += inv;
    start += shift;
    stop = start + L;
  }
  double mx = arma::max(arma::abs(res));
  if (mx > 0.0) res /= mx;
  return res;
}

// Forward decl to vat_srh_pitch_cpp from srh_pitch.cpp (returns List with f0, VUV, fs).
List vat_srh_pitch_cpp(NumericVector wave_r, double Fs, double F0min, double F0max);

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
NumericMatrix vat_creak_features_cpp(NumericVector x_r, double fs) {
  arma::vec x = as<arma::vec>(x_r);
  int N = x.n_elem;
  int winShift = (int)std::round(10.0 / 1000.0 * fs);

  // 1. LP residual + SRH pitch
  int lpcOrd = (int)std::round(fs / 1000.0) + 2;
  int L  = (int)std::round(25.0 / 1000.0 * fs);
  int sh = (int)std::round(5.0  / 1000.0 * fs);
  arma::vec res = lpc_residual_drugman(x, L, sh, lpcOrd);

  NumericVector x_nv(x.n_elem);
  for (arma::uword i = 0; i < x.n_elem; ++i) x_nv[i] = x(i);
  List ps = vat_srh_pitch_cpp(x_nv, fs, 20.0, 500.0);
  arma::vec f0  = as<arma::vec>(ps["f0"]);
  arma::ivec vuv = as<arma::ivec>(ps["VUV"]);
  double f0min = 20.0, f0max = 500.0;
  // F0mean = median of voiced f0 in valid range
  std::vector<double> vf;
  for (arma::uword i = 0; i < f0.n_elem; ++i)
    if (f0(i) > f0min && f0(i) < f0max && vuv(i) == 1) vf.push_back(f0(i));
  double F0mean = vf.empty() ? 100.0 : [&]{
    arma::vec v(vf.size()); for (size_t i=0;i<vf.size();++i) v(i)=vf[i]; return arma::median(v);
  }();

  // 2. sil_unv_features (32 ms frames as in MATLAB caller)
  SilUnvResult su = sil_unv_features(x, fs, 32.0);

  // 3. Short power std
  ShortPowResult sp = get_short_pow(x, fs);

  // 4. H2H1 + creakF0 (per-sample contours)
  H2H1Result h = get_creak_h2h1(res, fs, F0mean);

  // 5. Residual peak prominence (per-sample contour)
  arma::vec peak_inter = res_peak(x, fs, F0mean, res, su.Es);

  // Resample feature contours to 10 ms time grid
  // time = winShift:winShift:N (1-based) → in C++ 0-based: winShift-1, 2*winShift-1, ...
  std::vector<int> tgrid;
  for (int t = winShift; t <= N; t += winShift) tgrid.push_back(t - 1);
  int nf = tgrid.size();

  // enerN: Es_inter (per-sample) sampled at tgrid; ener_norm = Es - max(Es)
  arma::vec Es_norm = su.Es - arma::max(su.Es);
  arma::vec Es_inter = interp_to_length(Es_norm, N);
  // ZC interp
  arma::vec ZC_inter;
  if (su.Xpos.n_elem >= 2) {
    arma::vec xq = arma::regspace<arma::vec>(1, N);
    ZC_inter = vat::interp1_linear(su.Xpos, su.ZCs_ms, xq);
  } else {
    ZC_inter = arma::zeros<arma::vec>(N);
  }
  ZC_inter.replace(arma::datum::nan, 0.0);

  // F0 per-frame contour from SRH (interp directly to tgrid length, matching
  // MATLAB: F0 = interp1(linspace(1, n_frames, length(f0)), f0, 1:n_frames))
  arma::vec f0_per_frame;
  {
    arma::vec ts = arma::linspace<arma::vec>(1.0, double(nf), f0.n_elem);
    arma::vec xq = arma::regspace<arma::vec>(1, nf);
    f0_per_frame = vat::interp1_linear(ts, f0, xq);
  }

  // 6. Ishi features
  IshiResult ishi = get_ishi_params(x, fs);

  // ---- Assemble 12-base features at tgrid ----
  arma::mat F_static(nf, 12, arma::fill::zeros);
  for (int m = 0; m < nf; ++m) {
    int t = tgrid[m];
    int ti = std::min(N - 1, std::max(0, t));
    F_static(m, 0) = h.H2H1(ti);                  // col 1: H2H1
    F_static(m, 1) = peak_inter(ti);              // col 2: res_p
    F_static(m, 2) = ZC_inter(ti);                // col 3: ZCR
    F_static(m, 3) = ishi.IFP_inter(ti);          // col 4: IFP
    F_static(m, 4) = ishi.IPS_inter(ti);          // col 5: IPS
    F_static(m, 5) = ishi.PwP_fall_inter(ti);     // col 6: PwP.fall
    F_static(m, 6) = ishi.PwP_rise_inter(ti);     // col 7: PwP.rise
    F_static(m, 7) = f0_per_frame(m);             // col 8: F0
    F_static(m, 8) = F0mean;                      // col 9: F0mean
    F_static(m, 9) = Es_inter(ti);                // col 10: enerN
    F_static(m,10) = sp.pow_std_inter(ti);        // col 11: pow_std
    F_static(m,11) = h.creakF0(ti);               // col 12: creakF0
  }

  // Replace any NaN/Inf with 0 to keep ANN inputs finite
  F_static.replace(arma::datum::nan, 0.0);
  F_static.replace(arma::datum::inf, 0.0);
  F_static.replace(-arma::datum::inf, 0.0);

  arma::mat F_d  = delta_mat(F_static);
  arma::mat F_dd = delta_mat(F_d);
  arma::mat F = arma::join_horiz(F_static, F_d);
  F = arma::join_horiz(F, F_dd);  // n_frames x 36

  NumericMatrix out(F.n_rows, F.n_cols);
  for (arma::uword i = 0; i < F.n_rows; ++i)
    for (arma::uword j = 0; j < F.n_cols; ++j)
      out(i, j) = F(i, j);
  return out;
}

// Forward decl to vat_ann.cpp.
NumericVector vat_ann_forward_cpp(NumericMatrix X_r,
                                NumericMatrix IW_r, NumericVector b_h_r,
                                NumericMatrix LW_r, NumericVector b_o_r,
                                NumericVector mini_r, NumericVector maxi_r,
                                std::string out_act);

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List vat_creak_detect_cpp(NumericVector x_r, double fs,
                       NumericMatrix IW, NumericVector b_h,
                       NumericMatrix LW, NumericVector b_o,
                       NumericVector mini, NumericVector maxi,
                       double threshold = 0.3) {
  NumericMatrix feats = vat_creak_features_cpp(x_r, fs);
  // ann_forward expects n_features x n_frames
  NumericMatrix X(feats.ncol(), feats.nrow());
  for (int i = 0; i < feats.nrow(); ++i)
    for (int j = 0; j < feats.ncol(); ++j)
      X(j, i) = feats(i, j);
  NumericVector post = vat_ann_forward_cpp(X, IW, b_h, LW, b_o, mini, maxi, "logsig");
  IntegerVector dec(post.size());
  for (int i = 0; i < post.size(); ++i) dec[i] = (post[i] > threshold) ? 1 : 0;
  NumericVector tvec(post.size());
  for (int i = 0; i < post.size(); ++i) tvec[i] = (i + 1) * 0.010;
  return List::create(_["features"] = feats, _["posterior"] = post,
                      _["decision"] = dec, _["time"] = tvec);
}
