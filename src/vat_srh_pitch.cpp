// srh_pitch.cpp — SRH (Summation of Residual Harmonics) pitch tracking.
// Port of SRH_PitchTracking.m + SRH_EstimatePitch (Drugman & Alwan 2011).
//
// Pipeline: (optional) resample to 16 kHz -> LP residual -> 2-iter SRH search
// with median-F0-driven range refinement -> SRH-threshold VUV decisions.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "vat_dsp.h"
#include "vat_lpc.h"

using namespace Rcpp;

// Drugman's GetLPCresidual: 25 ms / 5 ms framing, Hanning window, autocorr LPC,
// FIR inverse filter, energy normalize per frame, overlap-add, scale by max|.|.
static arma::vec get_lpc_residual(const arma::vec& wave, int L, int shift, int order) {
  int N = wave.n_elem;
  arma::vec res(N, arma::fill::zeros);
  arma::vec hann_win = vat::hanning(L + 1);

  int start = 0, stop = start + L;
  while (stop < N) {
    arma::vec seg = wave.subvec(start, stop) % hann_win;
    arma::vec ar; double e;
    // MATLAB lpc() = autocorrelation Levinson (no Hamming windowing of autocorr).
    arma::vec r(order + 1, arma::fill::zeros);
    int n_seg = seg.n_elem;
    for (int k = 0; k <= order; ++k)
      for (int i = 0; i < n_seg - k; ++i) r(k) += seg(i) * seg(i + k);
    vat::levinson(r, order, ar, e);
    // FIR inverse: filter(A, 1, seg)
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

// One SRH iteration. sig: residual (column). Returns f0 (1-Hz resolution, padded
// with 5 zeros each side as in MATLAB), SRH amplitude per frame, time-start.
static void srh_estimate_pitch(const arma::vec& sig, double Fs,
                                int F0min, int F0max,
                                arma::vec& f0, arma::vec& SRHVal,
                                arma::ivec& tstart) {
  int N = sig.n_elem;
  int start = 0;
  int win  = (int)std::round(0.100 * Fs);   // 100 ms
  int shift = (int)std::round(0.010 * Fs);  // 10 ms
  int stop = win - 1;
  if (stop >= N) {
    f0.zeros(10); SRHVal.zeros(10); tstart.zeros(0);
    return;
  }
  int nframes = (N - win) / shift + 1;
  arma::vec frame_f0(nframes, arma::fill::zeros);
  arma::vec frame_sr(nframes, arma::fill::zeros);
  arma::ivec frame_t(nframes, arma::fill::zeros);
  arma::vec bw = vat::blackman(win);

  int nfft = (int)Fs;  // bin = 1 Hz
  arma::uword half = nfft / 2;

  for (int idx = 0; idx < nframes; ++idx) {
    int s0 = start + idx * shift;
    int s1 = s0 + win - 1;
    if (s1 >= N) break;
    arma::vec seg = sig.subvec(s0, s1) % bw;
    seg -= arma::mean(seg);

    arma::cx_vec X = vat::fft(seg, (arma::uword)nfft);
    arma::vec Spec(half);
    for (arma::uword i = 0; i < half; ++i) Spec(i) = std::abs(X(i));
    double nrm = std::sqrt(arma::dot(Spec, Spec));
    if (nrm > 0.0) Spec /= nrm;

    arma::vec srh_curve(F0max + 1, arma::fill::zeros);
    for (int f = F0min; f <= F0max; ++f) {
      auto S = [&](int k) {
        if (k < 0 || k >= (int)half) return 0.0;
        return Spec(k);
      };
      double pos = S(f) + S(2 * f) + S(3 * f) + S(4 * f) + S(5 * f);
      double neg = S((int)std::round(1.5 * f)) + S((int)std::round(2.5 * f)) +
                   S((int)std::round(3.5 * f)) + S((int)std::round(4.5 * f));
      srh_curve(f) = pos - neg;
    }
    arma::uword pos;
    srh_curve.max(pos);
    frame_f0(idx) = double(pos);
    frame_sr(idx) = srh_curve(pos);
    frame_t(idx) = s0 + 1;  // 1-indexed start time
  }

  // Pad with 5 zeros each side (MATLAB behavior)
  int pad = 5;
  int total = nframes + 2 * pad;
  f0.zeros(total);
  SRHVal.zeros(total);
  for (int i = 0; i < nframes; ++i) {
    f0(i + pad) = frame_f0(i);
    SRHVal(i + pad) = frame_sr(i);
  }
  tstart = frame_t;
}

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List vat_srh_pitch_cpp(NumericVector wave_r, double Fs,
                   double F0min, double F0max) {
  arma::vec wave = as<arma::vec>(wave_r);

  // Resample to 16 kHz if needed (rational p/q approx)
  double Fs_work = Fs;
  if (Fs > 16000.0) {
    // Approximate rational ratio: 16000/Fs
    auto gcd_fn = [](int a, int b){ while (b) { int t = b; b = a % b; a = t; } return a; };
    int g = gcd_fn((int)Fs, 16000);
    int q = (int)Fs / g;
    int p = 16000 / g;
    // Fallback for arbitrary fractional ratios
    if (std::abs(Fs * p / q - 16000.0) > 1e-6) {
      // simple decimation/interp: target 16000
      int factor = (int)std::round(Fs / 16000.0);
      arma::vec down((wave.n_elem + factor - 1) / factor);
      for (arma::uword i = 0, j = 0; i < wave.n_elem; i += factor, ++j)
        if (j < down.n_elem) down(j) = wave(i);
      wave = down;
    } else {
      wave = vat::resample(wave, p, q);
    }
    Fs_work = 16000.0;
  }

  int LPCorder = (int)std::round(0.75 * Fs_work / 1000.0);
  int L = (int)std::round(0.025 * Fs_work);
  int sh = (int)std::round(0.005 * Fs_work);
  arma::vec res = get_lpc_residual(wave, L, sh, LPCorder);

  arma::vec f0, SRHVal;
  arma::ivec t;
  int f0min_cur = (int)std::round(F0min);
  int f0max_cur = (int)std::round(F0max);

  for (int iter = 0; iter < 2; ++iter) {
    srh_estimate_pitch(res, Fs_work, f0min_cur, f0max_cur, f0, SRHVal, t);
    // Median F0 over frames where SRH > 0.1
    arma::uvec idx = arma::find(SRHVal > 0.1);
    if (idx.n_elem > 1) {
      arma::vec sub(idx.n_elem);
      for (arma::uword i = 0; i < idx.n_elem; ++i) sub(i) = f0(idx(i));
      double med = arma::median(sub);
      f0min_cur = (int)std::round(0.5 * med);
      f0max_cur = (int)std::round(2.0 * med);
      if (f0min_cur < 20) f0min_cur = 20;
      if (f0max_cur > (int)Fs_work / 2) f0max_cur = (int)Fs_work / 2;
    }
  }

  // VUV decision
  double sigma = arma::stddev(SRHVal);
  double thresh = (sigma > 0.05) ? 0.085 : 0.07;
  IntegerVector vuv(f0.n_elem);
  for (arma::uword i = 0; i < f0.n_elem; ++i) vuv[i] = (SRHVal(i) > thresh) ? 1 : 0;

  // Time vector in seconds (10 ms hop), padded
  arma::vec time(f0.n_elem);
  for (arma::uword i = 0; i < f0.n_elem; ++i) time(i) = double(i) * 0.010;

  NumericVector f0_out(f0.n_elem), srh_out(f0.n_elem), t_out(f0.n_elem);
  for (arma::uword i = 0; i < f0.n_elem; ++i) {
    f0_out[i]  = f0(i);
    srh_out[i] = SRHVal(i);
    t_out[i]   = time(i);
  }

  return List::create(
    _["f0"]     = f0_out,
    _["VUV"]    = vuv,
    _["SRHVal"] = srh_out,
    _["time"]   = t_out,
    _["fs"]     = Fs_work
  );
}
