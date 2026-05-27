// voice_quality.cpp — Voice quality parameter extraction.
// Ports get_NAQ_QOQ_H1H2.m + supporting helpers (integrat, findAmid_t).

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "vat_dsp.h"

using namespace Rcpp;

// integrat: y(n) = y(n-1) + x(n)/fs   (forward Euler accumulator)
static arma::vec integrat(const arma::vec& x, double fs) {
  arma::vec y(x.n_elem);
  double Ts = 1.0 / fs;
  y(0) = Ts * x(0);
  for (arma::uword n = 1; n < x.n_elem; ++n) y(n) = Ts * x(n) + y(n - 1);
  return y;
}

// findAmid_t: walk left/right from Tz while signal > Amid.
static void findAmid_t(const arma::vec& g, double Amid, int Tz, int& T1, int& T2) {
  T1 = T2 = 0;
  if (Tz <= 0) return;
  int n = Tz;
  while (n > 3 && g(n) > Amid) --n;
  T1 = n;
  n = Tz;
  int N = (int)g.n_elem;
  while (n < N - 3 && g(n) > Amid) ++n;
  T2 = n;
}

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List vat_naq_qoq_h1h2_cpp(NumericVector glot_r, double fs, IntegerVector gci_r) {
  arma::vec glot = as<arma::vec>(glot_r);
  arma::ivec gci = as<arma::ivec>(gci_r);
  int Ng = gci.n_elem;
  int Nx = glot.n_elem;
  arma::vec NAQ(Ng, arma::fill::zeros);
  arma::vec QOQ(Ng, arma::fill::zeros);
  arma::vec H1H2(Ng, arma::fill::zeros);
  arma::vec HRF(Ng, arma::fill::zeros);

  double F0min = 20, F0max = 500;
  double qoq_level = 0.5;
  int T0_num = 3;

  arma::vec glot_int = integrat(glot, fs);
  int glot_shift = (int)std::round(0.5 / 1000.0 * fs);

  for (int n = 0; n < Ng; ++n) {
    int start, stop, T0;
    if (n == 0) {
      start = 0;
      stop  = gci(0) - 1;
      if (Ng < 2) continue;
      T0 = gci(1) - gci(0);
    } else {
      start = gci(n - 1) - 1;
      stop  = gci(n) - 1;
      T0 = gci(n) - gci(n - 1);
    }
    if (start < 0) start = 0;
    if (stop >= Nx) stop = Nx - 1;
    if (T0 <= 0 || stop <= start) continue;

    double F0 = fs / double(T0);
    if (!std::isfinite(F0) || F0 <= F0min || F0 >= F0max) continue;

    // Linear baseline from glot_int(start) -> glot_int(stop)
    int span = stop - start + 1;
    arma::vec line(span);
    double a = glot_int(start), b = glot_int(stop);
    for (int i = 0; i < span; ++i)
      line(i) = a + (b - a) * double(i) / double(std::max(1, span - 1));
    arma::vec glot_int_cur = glot_int.subvec(start, stop);
    arma::vec glot_int_comp = glot_int_cur - line;

    int stop2 = (stop + glot_shift < Nx) ? stop + glot_shift : stop;
    arma::vec glot_cur = glot.subvec(start, stop2);

    // H1-H2 frame: ±T0*T0_num/2 around current GCI
    int gc = gci(n) - 1;
    int half = (int)std::round(T0 * T0_num / 2.0);
    int f_start = std::max(0, gc - half);
    int f_stop  = std::min(Nx - 1, gc + half);
    if (f_stop <= f_start) continue;

    arma::vec f_frame = glot.subvec(f_start, f_stop);
    arma::vec ham = vat::hamming(f_frame.n_elem);
    arma::vec f_win = f_frame % ham;
    arma::cx_vec X = vat::fft(f_win, (arma::uword)fs);
    int half_nf = (int)fs / 2;
    arma::vec f_spec(half_nf);
    for (int i = 0; i < half_nf; ++i)
      f_spec(i) = 20.0 * std::log10(std::abs(X(i)) + 1e-12);

    // Spectral peaks separated by at least F0/2 (in bins ~= Hz)
    vat::PeakResult pk = vat::findpeaks(f_spec, -1e300,
                                        (arma::uword)std::max(1, int(F0 / 2.0)));
    if (pk.locs.n_elem < 5) continue;

    // Closest peak to F0 and 2*F0
    arma::uword f0_idx = 0, f02_idx = 0;
    double best1 = 1e30, best2 = 1e30;
    for (arma::uword i = 0; i < pk.locs.n_elem; ++i) {
      double d1 = std::abs(double(pk.locs(i)) - F0);
      double d2 = std::abs(double(pk.locs(i)) - 2.0 * F0);
      if (d1 < best1) { best1 = d1; f0_idx = i; }
      if (d2 < best2) { best2 = d2; f02_idx = i; }
    }
    H1H2(n) = pk.pks(f0_idx) - pk.pks(f02_idx);
    if (pk.pks(0) != 0.0) HRF(n) = pk.pks(f0_idx) / pk.pks(0);

    // NAQ, QOQ
    double d_peak = arma::max(arma::abs(glot_cur));
    arma::uword max_idx;
    double f_ac = glot_int_comp.max(max_idx);
    double Amid = f_ac * qoq_level;
    int T1 = 0, T2 = 0;
    findAmid_t(glot_int_comp, Amid, (int)max_idx, T1, T2);
    if (d_peak > 0.0) NAQ(n) = (f_ac / d_peak) * F0;
    QOQ(n) = double(T2 - T1) / (fs / F0);
    if (QOQ(n) < 0 || QOQ(n) > 1) QOQ(n) = 0;
  }

  NumericVector naq_o(Ng), qoq_o(Ng), h1h2_o(Ng), hrf_o(Ng);
  for (int i = 0; i < Ng; ++i) {
    naq_o[i] = NAQ(i);
    qoq_o[i] = QOQ(i);
    h1h2_o[i] = H1H2(i);
    hrf_o[i] = HRF(i);
  }
  return List::create(_["NAQ"]=naq_o, _["QOQ"]=qoq_o,
                      _["H1H2"]=h1h2_o, _["HRF"]=hrf_o);
}

// ---------- Daless wavelet bank (peakSlope + MDQ shared) ----------
// Mother wavelet at scale s_n = 2^i:
//   t = (-1000:1000)/fs
//   h(t) = -cos(2*pi*f_o*(t/s_n)) * exp(-(t/s_n)^2 / (2*tau^2))
// where f_o = fs/2, tau = 1/(2*f_o).
static arma::vec daless_mw(int i_idx, double fs) {
  double s_n = std::pow(2.0, double(i_idx));
  double f_o = fs / 2.0;
  double tau = 1.0 / (2.0 * f_o);
  arma::vec t(2001);
  for (int k = 0; k <= 2000; ++k) t(k) = double(k - 1000) / fs;
  arma::vec u = t / s_n;
  arma::vec h(2001);
  for (arma::uword k = 0; k < 2001; ++k)
    h(k) = -std::cos(2.0 * arma::datum::pi * f_o * u(k)) *
            std::exp(-u(k) * u(k) / (2.0 * tau * tau));
  return h;
}

// Zero-phase convolution: full conv, then take center-aligned N samples.
static arma::vec daless_filt(const arma::vec& s, const arma::vec& h_i) {
  arma::vec full = arma::conv(s, h_i, "full");
  int halfLen = (int)std::ceil(double(h_i.n_elem) / 2.0);
  // MATLAB: y_i = s_h_i(halfLen : halfLen + length(s) - 1)
  return full.subvec(halfLen - 1, halfLen - 1 + s.n_elem - 1);
}

// Full multi-scale decomposition. i_set: vector of integer scale exponents.
static arma::mat daless_decomp(const arma::vec& s, double fs, const arma::ivec& i_set) {
  int K = i_set.n_elem;
  arma::mat y(K, s.n_elem, arma::fill::zeros);
  for (int n = 0; n < K; ++n) {
    arma::vec h_i = daless_mw(i_set(n), fs);
    arma::vec y_i = daless_filt(s, h_i);
    y.row(n) = y_i.t();
  }
  return y;
}

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List vat_daless_decomp_cpp(NumericVector s_r, double fs, IntegerVector i_r) {
  arma::vec s = as<arma::vec>(s_r);
  arma::ivec i_set = as<arma::ivec>(i_r);
  arma::mat y = daless_decomp(s, fs, i_set);
  return List::create(_["y"] = y);
}

// ---------- peakSlope ----------
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
NumericVector vat_peakslope_cpp(NumericVector s_r, double fs) {
  arma::vec s = as<arma::vec>(s_r);
  int frameLen   = (int)std::round(0.040 * fs);
  int frameShift = (int)std::round(0.010 * fs);
  arma::ivec i_set = arma::regspace<arma::ivec>(0, 6);  // 7 scales
  arma::mat y = daless_decomp(s, fs, i_set);

  int K = i_set.n_elem;
  int N = s.n_elem;
  int nframes = (N - frameLen) / frameShift + 1;
  if (nframes <= 0) return NumericVector(0);
  arma::vec ps(nframes, arma::fill::zeros);

  for (int m = 0; m < nframes; ++m) {
    int start = m * frameShift;
    int finish = start + frameLen - 1;
    if (finish >= N) break;
    arma::vec maxima(K);
    for (int k = 0; k < K; ++k)
      maxima(k) = arma::max(arma::abs(y.row(k).subvec(start, finish)));
    // Reverse to follow frequency order (low->high since i=0 is high-freq)
    arma::vec rev = arma::reverse(maxima);
    arma::vec logm = arma::log10(rev + 1e-12);
    // Linear regression slope
    int L = logm.n_elem;
    arma::vec tv(L);
    for (int k = 0; k < L; ++k) tv(k) = double(k + 1);
    double tm = arma::mean(tv), ym = arma::mean(logm);
    double num = 0.0, den = 0.0;
    for (int k = 0; k < L; ++k) {
      num += (tv(k) - tm) * (logm(k) - ym);
      den += (tv(k) - tm) * (tv(k) - tm);
    }
    ps(m) = (den > 0.0) ? num / den : 0.0;
  }
  NumericVector out(ps.n_elem);
  for (arma::uword i = 0; i < ps.n_elem; ++i) out[i] = std::isfinite(ps(i)) ? ps(i) : 0.0;
  return out;
}

// ---------- MDQ ----------
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
NumericVector vat_mdq_cpp(NumericVector res_r, double fs, IntegerVector gci_r) {
  arma::vec res = as<arma::vec>(res_r);
  arma::ivec gci = as<arma::ivec>(gci_r);
  // Deduplicate + sort
  std::vector<int> g;
  for (arma::uword i = 0; i < gci.n_elem; ++i)
    if (gci(i) >= 1 && gci(i) <= (int)res.n_elem) g.push_back(gci(i));
  std::sort(g.begin(), g.end());
  g.erase(std::unique(g.begin(), g.end()), g.end());
  int Ng = g.size();
  arma::vec MDQ(Ng, arma::fill::zeros);
  if (Ng < 2) return NumericVector(0);

  arma::ivec i_set = arma::regspace<arma::ivec>(2, 6);  // 5 scales
  arma::mat y = daless_decomp(res, fs, i_set);
  int K = i_set.n_elem;
  int N = res.n_elem;
  double searchRate = 0.2;

  for (int n = 0; n < Ng; ++n) {
    int T0 = (n == 0) ? g[1] - g[0] : g[n] - g[n - 1];
    if (T0 <= 0) continue;
    int gc0 = g[n] - 1;
    int searchLen = (int)std::round(searchRate * T0);
    int start_ser = gc0 - searchLen;
    int finish_ser = gc0 + searchLen;
    if (gc0 - T0 < 0 || gc0 + T0 >= N) continue;
    if (start_ser < 0 || finish_ser >= N) continue;

    int midpoint = searchLen;
    arma::vec dist_cur(K);
    for (int m = 0; m < K; ++m) {
      arma::rowvec row = y.row(m).subvec(start_ser, finish_ser);
      arma::uword maxIdx = row.index_max();
      dist_cur(m) = std::abs(double(midpoint) - double(maxIdx));
    }
    MDQ(n) = arma::mean(dist_cur) / double(T0);
  }

  for (arma::uword i = 0; i < MDQ.n_elem; ++i)
    if (!std::isfinite(MDQ(i))) MDQ(i) = 0.0;

  NumericVector out(MDQ.n_elem);
  for (arma::uword i = 0; i < MDQ.n_elem; ++i) out[i] = MDQ(i);
  return out;
}
