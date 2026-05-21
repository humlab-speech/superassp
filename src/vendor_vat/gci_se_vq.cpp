// gci_se_vq.cpp — Full SE-VQ GCI detection pipeline (Kane & Gobl 2013).
// Ports SE_VQ.m + SE_VQ_varF0.m and all supporting helpers
// (RCVD_reson_GCI, get_MBS, get_MBS_GCI_intervals, search_res_interval_peaks,
// GCI_creak_postproc, create_continuous_smooth_f0).

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "vat_dsp.h"
#include "vat_lpc.h"

using namespace Rcpp;

// ------- LP residual a la GetLPCresidual.m -------
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

// ------- RCVD_reson_GCI: filtfilt second-order resonator at F0mean -------
static arma::vec rcvd_reson_gci(const arma::vec& res, double fs, double F0mean) {
  double Phi = 2.0 * arma::datum::pi * F0mean / fs;
  double Rho = 0.9;
  arma::vec b = {1.0, 0.0, 0.0};
  arma::vec a = {1.0, -2.0 * Rho * std::cos(Phi), Rho * Rho};
  arma::vec y = vat::filtfilt(b, a, res);
  double mx = arma::max(arma::abs(y));
  if (mx > 0.0) y /= mx;
  return y;
}

// ------- get_MBS: mean-based signal -------
static arma::vec get_mbs(const arma::vec& x, double fs, const arma::vec& T0mean) {
  int N = x.n_elem;
  arma::vec MBS(N, arma::fill::zeros);
  // Clamp T0 to [fs/500, fs/20] to defend against Inf/NaN from f0_samp.
  double T0_max_clamp = fs / 20.0;
  double T0_min_clamp = fs / 500.0;
  auto clamp_T0 = [&](double t) {
    if (!std::isfinite(t) || t <= 0.0) t = fs / 100.0;
    if (t < T0_min_clamp) t = T0_min_clamp;
    if (t > T0_max_clamp) t = T0_max_clamp;
    return t;
  };
  int halfL = (int)std::round(1.6 * clamp_T0(T0mean(0)) / 2.0);
  int StepExp = 3;
  int Step = 1 << StepExp;

  for (int m = halfL; m < N - halfL; m += Step) {
    int hl;
    if (T0mean.n_elem == 1) hl = (int)std::round(1.7 * clamp_T0(T0mean(0)) / 2.0);
    else                    hl = (int)std::round(1.7 * clamp_T0(T0mean(m)) / 2.0);
    arma::vec bw = vat::blackman(2 * hl + 1);
    int start = m - hl, stop = m + hl;
    if (stop >= N) break;
    if (start < 0) continue;
    arma::vec vec = x.subvec(start, stop) % bw;
    MBS(m) = arma::mean(vec);
  }

  // Interpolate non-zero samples to fill the signal
  std::vector<double> tt, vv;
  for (int i = 0; i < N; ++i)
    if (MBS(i) != 0.0) { tt.push_back(double(i)); vv.push_back(MBS(i)); }
  if (tt.size() >= 2) {
    arma::vec txi = arma::regspace<arma::vec>(0, N - 1);
    arma::vec tvec(tt.size()), vvec(vv.size());
    for (size_t i = 0; i < tt.size(); ++i) { tvec(i) = tt[i]; vvec(i) = vv[i]; }
    MBS = vat::interp1_linear(tvec, vvec, txi);
  }
  MBS.replace(arma::datum::nan, 0.0);

  // zero-phase HPF at 70 Hz with 10 Hz transition (Butterworth)
  arma::vec b, a;
  vat::butter(4, std::min(0.99, 2.0 * 70.0 / fs), "high", b, a);
  MBS = vat::filtfilt(b, a, MBS);
  double mx = MBS.max();
  if (mx != 0.0) MBS /= mx;
  // Smooth with span 7 (MATLAB 'smooth' = moving average)
  arma::vec sm(N);
  int half = 3;
  for (int i = 0; i < N; ++i) {
    int lo = std::max(0, i - half), hi = std::min(N - 1, i + half);
    sm(i) = arma::mean(MBS.subvec(lo, hi));
  }
  return sm;
}

// ------- get_MBS_GCI_intervals: intervals around MBS negative peaks -------
static arma::imat get_mbs_intervals(const arma::vec& MBS, double fs,
                                     const arma::vec& T0mean, double F0max) {
  double F0max2 = F0max * 2.0;
  arma::uword T0max = (arma::uword)std::round(fs / F0max2);
  arma::vec neg = -MBS;
  vat::PeakResult pr = vat::findpeaks(neg, -1e300, std::max<arma::uword>(1, T0max));
  arma::uvec idx = pr.locs;
  int N = idx.n_elem;
  double search_rate = 0.28;
  double search_left = 0.01;
  arma::imat interval(N, 2, arma::fill::zeros);
  int Ms = MBS.n_elem;
  // Clamp T0 to a sane positive range (defends against NaN/Inf from f0_samp).
  double T0_max_clamp = fs / 20.0;
  double T0_min_clamp = fs / 500.0;
  for (int n = 0; n < N; ++n) {
    int i = idx(n);
    double T0 = (T0mean.n_elem > 1 && i < (int)T0mean.n_elem) ? T0mean(i) : T0mean(0);
    if (!std::isfinite(T0) || T0 <= 0.0) T0 = fs / 100.0;
    if (T0 < T0_min_clamp) T0 = T0_min_clamp;
    if (T0 > T0_max_clamp) T0 = T0_max_clamp;
    int start = i - (int)std::round(T0 * search_left);
    int stop  = i + (int)std::round(T0 * search_rate);
    if (start < 0) start = 0;
    if (stop >= Ms && start < Ms) stop = Ms - 1;
    else if (stop >= Ms) { N = n; break; }
    interval(n, 0) = start;
    interval(n, 1) = stop;
  }
  return interval.rows(0, N - 1);
}

// ------- search_res_interval_peaks: top Ncand peaks per interval -------
static void search_res_intervals(const arma::vec& res, const arma::imat& interval,
                                  int Ncand, arma::imat& GCI_N, arma::mat& GCI_relAmp) {
  int N = interval.n_rows;
  GCI_N.zeros(N, Ncand);
  GCI_relAmp.zeros(N, Ncand);
  for (int n = 0; n < N; ++n) {
    int start = interval(n, 0), stop = interval(n, 1);
    if (stop <= start) {
      GCI_N.row(n).fill(-1);
      continue;
    }
    arma::vec seg = res.subvec(start, stop);
    int M = seg.n_elem;
    if (M < Ncand) {
      arma::uword pos; double mx = seg.max(pos);
      (void)mx;
      GCI_N.row(n).fill((int)pos + start);
    } else {
      arma::uvec ord = arma::sort_index(seg, "descend");
      double mx = seg(ord(0));
      for (int c = 0; c < Ncand; ++c) {
        GCI_N(n, c) = (int)ord(c) + start;
        GCI_relAmp(n, c) = (mx > 0.0) ? 1.0 - seg(ord(c)) / mx : 0.0;
      }
    }
  }
}

// ------- GCI_creak_postproc -------
static arma::ivec gci_creak_postproc(const arma::ivec& GCI_in, const arma::ivec& creak,
                                      int search_reg, const arma::vec& rep,
                                      double removeThresh, int repNum) {
  std::vector<int> noncreak, creakGCI;
  for (arma::uword i = 0; i < GCI_in.n_elem; ++i) {
    int g = GCI_in(i);
    if (g >= 0 && g < (int)creak.n_elem && creak(g) == 1) creakGCI.push_back(g);
    else                                                   noncreak.push_back(g);
  }
  for (int m = 0; m < repNum; ++m) {
    size_t n = 1;
    while (n + 1 < creakGCI.size()) {
      int g = creakGCI[n];
      int lo = std::max(0, g - search_reg);
      int hi = std::min((int)rep.n_elem - 1, g + search_reg);
      double cur = std::abs(arma::min(rep.subvec(lo, hi)));
      double nbr = 0.5 * (std::abs(rep(creakGCI[n - 1])) + std::abs(rep(creakGCI[n + 1])));
      if (nbr * removeThresh > cur) {
        creakGCI.erase(creakGCI.begin() + n);
      } else {
        ++n;
      }
    }
  }
  std::vector<int> all = noncreak;
  all.insert(all.end(), creakGCI.begin(), creakGCI.end());
  std::sort(all.begin(), all.end());
  all.erase(std::unique(all.begin(), all.end()), all.end());
  arma::ivec out(all.size());
  for (size_t i = 0; i < all.size(); ++i) out(i) = all[i];
  return out;
}

// ------- create_continuous_smooth_f0 -------
static arma::vec create_continuous_smooth_f0(const arma::vec& f0_in, const arma::ivec& VUV,
                                              int Nsamples) {
  double F0_min = 50, F0_max = 500;
  int med_len = 17, sm_len = 17;

  arma::vec f0 = f0_in;
  // collect voiced f0 in range
  std::vector<double> v;
  for (arma::uword i = 0; i < f0.n_elem; ++i)
    if (VUV(i) == 1 && f0(i) > F0_min && f0(i) < F0_max) v.push_back(f0(i));
  if (v.empty()) {
    arma::vec out(Nsamples); out.fill(150.0);
    return out;
  }
  arma::vec vv(v.size()); for (size_t i = 0; i < v.size(); ++i) vv(i) = v[i];
  arma::vec vv_med = vat::medfilt1(vv, med_len);
  double F0lo = vv_med.min(), F0hi = vv_med.max();
  f0.transform([&](double x) { return std::min(F0hi, std::max(F0lo, x)); });

  // Treat VUV==0 as unvoiced regions; linearly interpolate boundaries
  for (arma::uword i = 0; i < f0.n_elem; ++i)
    if (VUV(i) == 0) f0(i) = arma::datum::nan;
  // Forward/backward fill via linear interp over non-NaN indices
  std::vector<double> xs, ys;
  for (arma::uword i = 0; i < f0.n_elem; ++i)
    if (!std::isnan(f0(i))) { xs.push_back(double(i)); ys.push_back(f0(i)); }
  if (xs.empty()) {
    arma::vec out(Nsamples); out.fill(150.0);
    return out;
  }
  arma::vec xq = arma::regspace<arma::vec>(0, f0.n_elem - 1);
  arma::vec xv(xs.size()), yv(ys.size());
  for (size_t i = 0; i < xs.size(); ++i) { xv(i) = xs[i]; yv(i) = ys[i]; }
  arma::vec f0_inter = vat::interp1_linear(xv, yv, xq);

  // Median + moving-average smoothing
  f0_inter = vat::medfilt1(f0_inter, med_len);
  arma::vec sm(f0_inter.n_elem);
  int half = sm_len / 2;
  for (arma::uword i = 0; i < f0_inter.n_elem; ++i) {
    int lo = std::max(0, int(i) - half);
    int hi = std::min(int(f0_inter.n_elem) - 1, int(i) + half);
    sm(i) = arma::mean(f0_inter.subvec(lo, hi));
  }

  // Resample to length Nsamples
  arma::vec xs2 = arma::linspace<arma::vec>(0.0, double(Nsamples - 1), sm.n_elem);
  arma::vec xq2 = arma::regspace<arma::vec>(0, Nsamples - 1);
  return vat::interp1_linear(xs2, sm, xq2);
}

// Re-use Viterbi DP from se_vq_dp.cpp (declared here).
IntegerVector vat_reson_dynprog_cpp(NumericMatrix gci_rel_amp, IntegerMatrix gci_n,
                                 double f0_mean, NumericVector x, double fs,
                                 double trans_wgt, double rel_amp_wgt);

// ------- Top-level SE-VQ (and SE-VQ_varF0 toggle) -------
// [[Rcpp::export]]
List vat_se_vq_cpp(NumericVector x_r, double fs,
                NumericVector f0_r, IntegerVector vuv_r,
                Nullable<IntegerVector> creak_r = R_NilValue,
                bool var_f0 = false) {
  arma::vec x = as<arma::vec>(x_r);
  arma::vec f0 = as<arma::vec>(f0_r);
  arma::ivec vuv = as<arma::ivec>(vuv_r);

  // Settings
  double F0min_set = 20, F0max_set = 500;
  std::vector<double> vf;
  for (arma::uword i = 0; i < f0.n_elem; ++i)
    if (f0(i) > F0min_set && f0(i) < F0max_set && vuv(i) == 1) vf.push_back(f0(i));
  double F0mean = vf.empty() ? 100.0 :
    [&]{ arma::vec v(vf.size()); for (size_t i = 0; i < vf.size(); ++i) v(i) = vf[i];
         return arma::median(v); }();
  double F0maxObs = vf.empty() ? 500.0 :
    [&]{ arma::vec v(vf.size()); for (size_t i = 0; i < vf.size(); ++i) v(i) = vf[i];
         return arma::max(vat::medfilt1(v, 13)); }();
  if (F0mean < 70.0) F0mean = 80.0;

  double T0mean_scalar = fs / F0mean;
  arma::vec T0mean;
  if (var_f0) {
    arma::vec f0_samp = create_continuous_smooth_f0(f0, vuv, x.n_elem);
    T0mean = fs / f0_samp;
  } else {
    T0mean = arma::vec({T0mean_scalar});
  }

  double winLen_ms = 25, winShift_ms = 5;
  int L = (int)std::round(winLen_ms / 1000.0 * fs);
  int sh = (int)std::round(winShift_ms / 1000.0 * fs);
  int LPC_ord = (int)std::round(fs / 1000.0) + 2;
  int Ncand = 5;
  double trans_wgt = 1.0, relAmp_wgt = 0.3;

  arma::vec res = lpc_residual_drugman(x, L, sh, LPC_ord);
  arma::vec rep = rcvd_reson_gci(res, fs, F0mean);
  arma::vec MBS = get_mbs(x, fs, T0mean);
  arma::imat interval = get_mbs_intervals(MBS, fs, T0mean, F0maxObs);
  arma::imat GCI_N; arma::mat GCI_relAmp;
  search_res_intervals(res, interval, Ncand, GCI_N, GCI_relAmp);

  // To re-use existing Viterbi DP (NumericMatrix/IntegerMatrix), wrap arma -> Rcpp.
  int N = GCI_N.n_rows;
  NumericMatrix amp_m(Ncand, N);
  IntegerMatrix n_m(Ncand, N);
  for (int n = 0; n < N; ++n)
    for (int c = 0; c < Ncand; ++c) {
      amp_m(c, n) = GCI_relAmp(n, c);
      n_m(c, n)   = GCI_N(n, c) + 1;  // DP expects 1-indexed positions
    }
  NumericVector x_nv(x.n_elem);
  for (arma::uword i = 0; i < x.n_elem; ++i) x_nv[i] = x(i);
  IntegerVector gci_opt = vat_reson_dynprog_cpp(amp_m, n_m, F0mean, x_nv, fs,
                                             trans_wgt, relAmp_wgt);

  // Convert to 0-indexed arma::ivec for postproc
  arma::ivec gci0(gci_opt.size());
  for (int i = 0; i < gci_opt.size(); ++i) gci0(i) = gci_opt[i] - 1;

  // Optional creaky postproc
  if (creak_r.isNotNull()) {
    IntegerVector creak_nv(creak_r);
    if ((arma::uword)creak_nv.size() == x.n_elem) {
      arma::ivec creak = as<arma::ivec>(creak_nv);
      int search_reg = (int)std::round(1.3 / 1000.0 * fs);
      gci0 = gci_creak_postproc(gci0, creak, search_reg, rep, 0.4, 2);
    }
  }

  // Sort + unique
  std::vector<int> g;
  for (arma::uword i = 0; i < gci0.n_elem; ++i)
    if (gci0(i) >= 0 && gci0(i) < (int)x.n_elem) g.push_back(gci0(i));
  std::sort(g.begin(), g.end());
  g.erase(std::unique(g.begin(), g.end()), g.end());
  IntegerVector gci_out(g.size());
  for (size_t i = 0; i < g.size(); ++i) gci_out[i] = g[i] + 1;  // 1-indexed for R/MATLAB

  // Wrap residual / resonator / MBS for downstream use
  NumericVector res_nv(res.n_elem), rep_nv(rep.n_elem), mbs_nv(MBS.n_elem);
  for (arma::uword i = 0; i < res.n_elem; ++i) { res_nv[i] = res(i); rep_nv[i] = rep(i); mbs_nv[i] = MBS(i); }

  return List::create(
    _["GCI"] = gci_out,
    _["res"] = res_nv,
    _["rep"] = rep_nv,
    _["MBS"] = mbs_nv,
    _["F0mean"] = F0mean,
    _["F0max"] = F0maxObs
  );
}
