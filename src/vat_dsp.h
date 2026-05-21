// vat_dsp.h — Shared DSP primitives for the Voice Analysis Toolkit Rcpp port.
// Header-only. MATLAB-faithful semantics.

#ifndef VAT_DSP_H
#define VAT_DSP_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <algorithm>
#include <complex>

namespace vat {

// ---------- Windows ----------
inline arma::vec hamming(arma::uword N, bool symmetric = true) {
  if (N == 0) return arma::vec();
  if (N == 1) return arma::vec({1.0});
  arma::vec w(N);
  double denom = symmetric ? double(N - 1) : double(N);
  for (arma::uword n = 0; n < N; ++n)
    w(n) = 0.54 - 0.46 * std::cos(2.0 * arma::datum::pi * double(n) / denom);
  return w;
}

inline arma::vec hanning(arma::uword N, bool symmetric = true) {
  if (N == 0) return arma::vec();
  if (N == 1) return arma::vec({1.0});
  arma::vec w(N);
  double denom = symmetric ? double(N - 1) : double(N);
  for (arma::uword n = 0; n < N; ++n)
    w(n) = 0.5 - 0.5 * std::cos(2.0 * arma::datum::pi * double(n) / denom);
  return w;
}

inline arma::vec blackman(arma::uword N, bool symmetric = true) {
  if (N == 0) return arma::vec();
  if (N == 1) return arma::vec({1.0});
  arma::vec w(N);
  double denom = symmetric ? double(N - 1) : double(N);
  for (arma::uword n = 0; n < N; ++n) {
    double t = 2.0 * arma::datum::pi * double(n) / denom;
    w(n) = 0.42 - 0.5 * std::cos(t) + 0.08 * std::cos(2.0 * t);
  }
  return w;
}

// MATLAB-compatible Kaiser window (used by resample).
inline double bessel_i0(double x) {
  double ax = std::fabs(x);
  if (ax < 3.75) {
    double y = x / 3.75; y = y * y;
    return 1.0 + y*(3.5156229 + y*(3.0899424 + y*(1.2067492
              + y*(0.2659732 + y*(0.0360768 + y*0.0045813)))));
  }
  double y = 3.75 / ax;
  return (std::exp(ax)/std::sqrt(ax))*(0.39894228 + y*(0.01328592
        + y*(0.00225319 + y*(-0.00157565 + y*(0.00916281
        + y*(-0.02057706 + y*(0.02635537 + y*(-0.01647633 + y*0.00392377))))))));
}

inline arma::vec kaiser(arma::uword N, double beta) {
  if (N == 0) return arma::vec();
  if (N == 1) return arma::vec({1.0});
  arma::vec w(N);
  double denom = bessel_i0(beta);
  double half = double(N - 1) / 2.0;
  for (arma::uword n = 0; n < N; ++n) {
    double r = (double(n) - half) / half;
    w(n) = bessel_i0(beta * std::sqrt(std::max(0.0, 1.0 - r*r))) / denom;
  }
  return w;
}

// ---------- IIR filter (Direct-Form-II Transposed) ----------
// MATLAB-equivalent y = filter(b, a, x [, zi]).
// zi length max(len(b), len(a)) - 1.
inline arma::vec filter(arma::vec b, arma::vec a, const arma::vec& x,
                        const arma::vec& zi = arma::vec()) {
  if (a.n_elem == 0 || a(0) == 0.0)
    Rcpp::stop("filter: a(0) must be nonzero");
  double a0 = a(0);
  if (a0 != 1.0) { b /= a0; a /= a0; }

  arma::uword nb = b.n_elem, na = a.n_elem;
  arma::uword nfilt = std::max(nb, na);
  if (nb < nfilt) b.resize(nfilt);
  if (na < nfilt) a.resize(nfilt);

  arma::vec z(nfilt - 1, arma::fill::zeros);
  if (zi.n_elem >= nfilt - 1)
    z = zi.subvec(0, nfilt - 2);

  arma::uword nx = x.n_elem;
  arma::vec y(nx);
  for (arma::uword n = 0; n < nx; ++n) {
    double xn = x(n);
    double yn = b(0) * xn + (nfilt > 1 ? z(0) : 0.0);
    y(n) = yn;
    for (arma::uword k = 1; k < nfilt - 1; ++k)
      z(k - 1) = b(k) * xn + z(k) - a(k) * yn;
    if (nfilt > 1)
      z(nfilt - 2) = b(nfilt - 1) * xn - a(nfilt - 1) * yn;
  }
  return y;
}

// Steady-state initial conditions for filtfilt; solves (I - A) zi = B - a(1)*b(0:end)
// where A and B are the companion-form matrices. Matches MATLAB filtfilt.m.
inline arma::vec filtfilt_initial_state(const arma::vec& b, const arma::vec& a) {
  arma::uword n = std::max(b.n_elem, a.n_elem);
  arma::vec bb = b, aa = a;
  if (bb.n_elem < n) bb.resize(n);
  if (aa.n_elem < n) aa.resize(n);
  if (n < 2) return arma::vec();

  arma::mat M = arma::eye<arma::mat>(n - 1, n - 1);
  // Subtract A_companion: first column -a(1..end); identity sub-diagonal
  for (arma::uword i = 0; i < n - 1; ++i)
    M(i, 0) += aa(i + 1);
  for (arma::uword i = 0; i < n - 2; ++i)
    M(i, i + 1) -= 1.0;

  arma::vec rhs(n - 1);
  for (arma::uword i = 0; i < n - 1; ++i)
    rhs(i) = bb(i + 1) - aa(i + 1) * bb(0);

  arma::vec zi;
  if (!arma::solve(zi, M, rhs))
    zi = arma::zeros<arma::vec>(n - 1);
  return zi;
}

// Zero-phase filtering. MATLAB filtfilt(b, a, x).
inline arma::vec filtfilt(const arma::vec& b, const arma::vec& a, const arma::vec& x) {
  arma::uword nfilt = std::max(b.n_elem, a.n_elem);
  if (nfilt < 2) return x;
  arma::uword nfact = 3 * (nfilt - 1);
  if (x.n_elem <= nfact)
    Rcpp::stop("filtfilt: input must be > 3*(filter_order)");

  // Edge reflection: 2*x(0) - x(nfact:-1:1) and 2*x(end) - x(end-1:-1:end-nfact)
  arma::vec pre(nfact), post(nfact);
  for (arma::uword i = 0; i < nfact; ++i) {
    pre(i)  = 2.0 * x(0)            - x(nfact - i);
    post(i) = 2.0 * x(x.n_elem - 1) - x(x.n_elem - 2 - i);
  }
  arma::vec ext = arma::join_vert(arma::join_vert(pre, x), post);

  arma::vec zi = filtfilt_initial_state(b, a);

  arma::vec y1 = filter(b, a, ext, zi * ext(0));
  arma::vec y1r = arma::reverse(y1);
  arma::vec y2 = filter(b, a, y1r, zi * y1r(0));
  arma::vec y  = arma::reverse(y2);

  return y.subvec(nfact, y.n_elem - nfact - 1);
}

// ---------- FIR design (fir1) ----------
// Hamming-windowed-sinc, length n+1, symmetric. type: "low"/"high"/"stop"/"pass".
// Wn in [0,1] where 1 = Nyquist.
inline arma::vec fir1(arma::uword n, arma::vec Wn, std::string type = "low") {
  if (n < 1) Rcpp::stop("fir1: order must be >= 1");
  arma::uword L = n + 1;
  bool odd = (n % 2 == 1);
  // Wn sorted ascending
  Wn = arma::sort(Wn);

  // Bands: alternate stop/pass starting from "type"
  std::vector<double> edges{0.0};
  for (arma::uword i = 0; i < Wn.n_elem; ++i) edges.push_back(Wn(i));
  edges.push_back(1.0);

  std::vector<int> mags;
  bool start_pass;
  if (type == "low" || type == "bandpass" || type == "pass")
    start_pass = (Wn.n_elem == 1 || Wn.n_elem == 2);
  else
    start_pass = false;
  if (type == "low")  start_pass = true;
  if (type == "high") start_pass = false;
  if (type == "stop") start_pass = true;
  if (type == "bandpass" || type == "pass") start_pass = false;

  int cur = start_pass ? 1 : 0;
  for (size_t i = 0; i + 1 < edges.size(); ++i) {
    mags.push_back(cur);
    cur = 1 - cur;
  }

  // MATLAB fir1: highpass/stop require Type I (odd length L => even order n).
  if (odd && (type == "high" || type == "stop"))
    Rcpp::stop("fir1: highpass/stop requires even order (odd length)");

  arma::vec h(L, arma::fill::zeros);
  double half = double(n) / 2.0;
  for (size_t k = 0; k < mags.size(); ++k) {
    if (mags[k] == 0) continue;
    double f1 = edges[k], f2 = edges[k + 1];
    for (arma::uword m = 0; m < L; ++m) {
      double t = double(m) - half;
      if (std::fabs(t) < 1e-12)
        h(m) += (f2 - f1);
      else
        h(m) += (std::sin(arma::datum::pi * f2 * t)
               - std::sin(arma::datum::pi * f1 * t)) / (arma::datum::pi * t);
    }
  }
  arma::vec w = hamming(L);
  h %= w;

  // Normalize so passband gain = 1 at center of first passband
  double f_norm = 0.0;
  for (size_t k = 0; k < mags.size(); ++k) {
    if (mags[k] == 1) { f_norm = 0.5 * (edges[k] + edges[k + 1]); break; }
  }
  arma::cx_double resp(0.0, 0.0);
  for (arma::uword m = 0; m < L; ++m)
    resp += h(m) * std::exp(arma::cx_double(0.0, -arma::datum::pi * f_norm * double(m)));
  double g = std::abs(resp);
  if (g > 0.0) h /= g;
  return h;
}

inline arma::vec fir1(arma::uword n, double Wn, std::string type = "low") {
  arma::vec v({Wn});
  return fir1(n, v, type);
}

// ---------- Butterworth (analog prototype + bilinear) ----------
// Returns [b, a] as 2-row matrix (row 0 = b, row 1 = a). type: low/high/bandpass/stop.
inline void butter(int n, arma::vec Wn, std::string type,
                   arma::vec& b_out, arma::vec& a_out) {
  Wn = arma::sort(Wn);
  // Pre-warp
  arma::vec u(Wn.n_elem);
  for (arma::uword i = 0; i < Wn.n_elem; ++i)
    u(i) = 2.0 * std::tan(arma::datum::pi * Wn(i) / 2.0);

  // Analog low-pass prototype poles
  std::vector<arma::cx_double> pa(n);
  for (int k = 0; k < n; ++k) {
    double th = arma::datum::pi * (2.0 * (k + 1) + n - 1) / (2.0 * n);
    pa[k] = arma::cx_double(std::cos(th), std::sin(th));
  }
  arma::cx_double ka(1.0, 0.0);
  for (int k = 0; k < n; ++k) ka *= -pa[k];

  // Frequency transform
  std::vector<arma::cx_double> pa_t, za_t;
  arma::cx_double k_t = ka;
  if (type == "low") {
    double Wn1 = u(0);
    for (auto& p : pa) pa_t.push_back(p * Wn1);
    k_t *= std::pow(Wn1, n);
  } else if (type == "high") {
    double Wn1 = u(0);
    for (auto& p : pa) {
      pa_t.push_back(arma::cx_double(Wn1, 0.0) / p);
      za_t.push_back(arma::cx_double(0.0, 0.0));
    }
    for (auto& p : pa) k_t /= -p;
  } else if (type == "bandpass" || type == "pass") {
    double Bw = u(1) - u(0);
    double Wo = std::sqrt(u(0) * u(1));
    for (auto& p : pa) {
      arma::cx_double pb = p * Bw / 2.0;
      arma::cx_double disc = std::sqrt(pb * pb - Wo * Wo);
      pa_t.push_back(pb + disc);
      pa_t.push_back(pb - disc);
      za_t.push_back(arma::cx_double(0.0, 0.0));
    }
    k_t *= std::pow(Bw, n);
  } else if (type == "stop") {
    double Bw = u(1) - u(0);
    double Wo = std::sqrt(u(0) * u(1));
    for (auto& p : pa) {
      arma::cx_double pb = (Bw / 2.0) / p;
      arma::cx_double disc = std::sqrt(pb * pb - Wo * Wo);
      pa_t.push_back(pb + disc);
      pa_t.push_back(pb - disc);
      za_t.push_back(arma::cx_double(0.0,  Wo));
      za_t.push_back(arma::cx_double(0.0, -Wo));
    }
    for (auto& p : pa) k_t /= -p;
  } else {
    Rcpp::stop("butter: unknown type");
  }

  // Bilinear transform: z = (2+s)/(2-s)
  std::vector<arma::cx_double> pd, zd;
  arma::cx_double kd = k_t;
  for (auto& p : pa_t) {
    pd.push_back((2.0 + p) / (2.0 - p));
    kd /= (2.0 - p);
  }
  arma::uword n_pa = pa_t.size(), n_za = za_t.size();
  for (auto& z : za_t) {
    zd.push_back((2.0 + z) / (2.0 - z));
    kd *= (2.0 - z);
  }
  while (zd.size() < pd.size()) { zd.push_back(arma::cx_double(-1.0, 0.0)); }
  (void)n_pa; (void)n_za;

  // Polynomials from roots
  auto poly = [](const std::vector<arma::cx_double>& r) {
    arma::cx_vec p(1); p(0) = arma::cx_double(1.0, 0.0);
    for (auto& ri : r) {
      arma::cx_vec p2(p.n_elem + 1, arma::fill::zeros);
      for (arma::uword i = 0; i < p.n_elem; ++i) {
        p2(i)     += p(i);
        p2(i + 1) -= ri * p(i);
      }
      p = p2;
    }
    return p;
  };
  arma::cx_vec az = poly(pd);
  arma::cx_vec bz = poly(zd);
  bz *= kd;

  b_out.set_size(bz.n_elem); a_out.set_size(az.n_elem);
  for (arma::uword i = 0; i < bz.n_elem; ++i) b_out(i) = bz(i).real();
  for (arma::uword i = 0; i < az.n_elem; ++i) a_out(i) = az(i).real();
}

inline void butter(int n, double Wn, std::string type,
                   arma::vec& b_out, arma::vec& a_out) {
  arma::vec v({Wn});
  butter(n, v, type, b_out, a_out);
}

// ---------- Median filter (medfilt1, zero-pad) ----------
inline arma::vec medfilt1(const arma::vec& x, arma::uword n) {
  if (n < 2) return x;
  arma::uword N = x.n_elem;
  arma::vec y(N);
  arma::uword half = n / 2;
  std::vector<double> buf(n);
  for (arma::uword i = 0; i < N; ++i) {
    for (arma::uword k = 0; k < n; ++k) {
      long idx = long(i) + long(k) - long(half);
      buf[k] = (idx < 0 || idx >= long(N)) ? 0.0 : x(idx);
    }
    std::sort(buf.begin(), buf.end());
    if (n % 2 == 1) y(i) = buf[half];
    else            y(i) = 0.5 * (buf[half - 1] + buf[half]);
  }
  return y;
}

// ---------- Linear interpolation (MATLAB interp1, "linear", extrap allowed) ----------
inline arma::vec interp1_linear(const arma::vec& x, const arma::vec& y,
                                const arma::vec& xq, double extrap = 0.0) {
  arma::uword n = x.n_elem, nq = xq.n_elem;
  arma::vec yq(nq);
  for (arma::uword i = 0; i < nq; ++i) {
    double xi = xq(i);
    if (xi <= x(0)) {
      if (std::isnan(extrap)) yq(i) = std::nan("");
      else {
        double dx = x(1) - x(0);
        yq(i) = (dx != 0.0) ? y(0) + (xi - x(0)) * (y(1) - y(0)) / dx : y(0);
      }
      continue;
    }
    if (xi >= x(n - 1)) {
      if (std::isnan(extrap)) yq(i) = std::nan("");
      else {
        double dx = x(n - 1) - x(n - 2);
        yq(i) = (dx != 0.0)
              ? y(n - 1) + (xi - x(n - 1)) * (y(n - 1) - y(n - 2)) / dx
              : y(n - 1);
      }
      continue;
    }
    // Binary search
    arma::uword lo = 0, hi = n - 1;
    while (hi - lo > 1) {
      arma::uword m = (lo + hi) / 2;
      if (x(m) <= xi) lo = m; else hi = m;
    }
    double dx = x(hi) - x(lo);
    yq(i) = (dx != 0.0) ? y(lo) + (xi - x(lo)) * (y(hi) - y(lo)) / dx : y(lo);
  }
  return yq;
}

// Natural cubic spline interpolation (MATLAB "spline" mode).
inline arma::vec interp1_spline(const arma::vec& x, const arma::vec& y,
                                const arma::vec& xq) {
  arma::uword n = x.n_elem;
  arma::vec h(n - 1), alpha(n - 1, arma::fill::zeros);
  for (arma::uword i = 0; i < n - 1; ++i) h(i) = x(i + 1) - x(i);
  for (arma::uword i = 1; i < n - 1; ++i)
    alpha(i) = 3.0 * ((y(i + 1) - y(i)) / h(i) - (y(i) - y(i - 1)) / h(i - 1));

  arma::vec l(n), mu(n, arma::fill::zeros), z(n, arma::fill::zeros);
  l(0) = 1.0;
  for (arma::uword i = 1; i < n - 1; ++i) {
    l(i)  = 2.0 * (x(i + 1) - x(i - 1)) - h(i - 1) * mu(i - 1);
    mu(i) = h(i) / l(i);
    z(i)  = (alpha(i) - h(i - 1) * z(i - 1)) / l(i);
  }
  l(n - 1) = 1.0;
  arma::vec c(n, arma::fill::zeros), b(n - 1), d(n - 1);
  for (long j = long(n) - 2; j >= 0; --j) {
    c(j) = z(j) - mu(j) * c(j + 1);
    b(j) = (y(j + 1) - y(j)) / h(j) - h(j) * (c(j + 1) + 2.0 * c(j)) / 3.0;
    d(j) = (c(j + 1) - c(j)) / (3.0 * h(j));
  }

  arma::vec yq(xq.n_elem);
  for (arma::uword i = 0; i < xq.n_elem; ++i) {
    double xi = xq(i);
    arma::uword lo = 0, hi = n - 1;
    if (xi <= x(0)) lo = 0;
    else if (xi >= x(n - 1)) lo = n - 2;
    else {
      while (hi - lo > 1) {
        arma::uword m = (lo + hi) / 2;
        if (x(m) <= xi) lo = m; else hi = m;
      }
    }
    double dx = xi - x(lo);
    yq(i) = y(lo) + b(lo) * dx + c(lo) * dx * dx + d(lo) * dx * dx * dx;
  }
  return yq;
}

// ---------- Find peaks (subset of MATLAB findpeaks) ----------
struct PeakResult { arma::uvec locs; arma::vec pks; };
inline PeakResult findpeaks(const arma::vec& x,
                             double min_peak_height = -arma::datum::inf,
                             arma::uword min_peak_distance = 1) {
  arma::uword N = x.n_elem;
  std::vector<arma::uword> raw;
  for (arma::uword i = 1; i + 1 < N; ++i)
    if (x(i) > x(i - 1) && x(i) > x(i + 1) && x(i) >= min_peak_height)
      raw.push_back(i);
  // Sort by peak height descending, then suppress within distance
  std::sort(raw.begin(), raw.end(),
            [&](arma::uword a, arma::uword b) { return x(a) > x(b); });
  std::vector<bool> kept(raw.size(), true);
  for (arma::uword i = 0; i < raw.size(); ++i) {
    if (!kept[i]) continue;
    for (arma::uword j = i + 1; j < raw.size(); ++j) {
      if (kept[j] &&
          (raw[j] > raw[i] ? raw[j] - raw[i] : raw[i] - raw[j]) < min_peak_distance)
        kept[j] = false;
    }
  }
  std::vector<arma::uword> survivors;
  for (arma::uword i = 0; i < raw.size(); ++i)
    if (kept[i]) survivors.push_back(raw[i]);
  std::sort(survivors.begin(), survivors.end());
  PeakResult r;
  r.locs.set_size(survivors.size());
  r.pks.set_size(survivors.size());
  for (arma::uword i = 0; i < survivors.size(); ++i) {
    r.locs(i) = survivors[i];
    r.pks(i)  = x(survivors[i]);
  }
  return r;
}

// ---------- Resample (rational p/q via polyphase + Kaiser FIR) ----------
inline arma::vec resample(const arma::vec& x, int p, int q, double beta = 5.0) {
  if (p <= 0 || q <= 0) Rcpp::stop("resample: p,q must be positive");
  auto gcd_fn = [](int a, int b){ while (b) { int t = b; b = a % b; a = t; } return a; };
  int g = gcd_fn(p, q);
  p /= g; q /= g;
  if (p == 1 && q == 1) return x;
  int Lx = x.n_elem;
  int N = 10;  // half-length factor (MATLAB default)
  double fc = 1.0 / double(std::max(p, q));
  int L = 2 * N * std::max(p, q) + 1;
  arma::vec h = fir1(L - 1, fc, "low");
  // Scale so polyphase upsample gain = p
  h *= double(p);
  // Apply Kaiser window
  arma::vec kw = kaiser(L, beta);
  h %= kw;

  // Upsample by p, filter, downsample by q
  int Ly = int(std::ceil(double(Lx) * p / double(q)));
  arma::vec y(Ly, arma::fill::zeros);
  int half = (L - 1) / 2;
  for (int n = 0; n < Ly; ++n) {
    long m = long(n) * q;
    long lo = m - half, hi = m + half;
    double acc = 0.0;
    for (long k = lo; k <= hi; ++k) {
      if (k % p != 0 && k >= 0) continue;
      long xi = k / p;
      if (xi >= 0 && xi < Lx) acc += h(k - lo) * x(xi);
    }
    y(n) = acc;
  }
  return y;
}

// ---------- FFT helpers (MATLAB-compatible scaling) ----------
inline arma::cx_vec fft(const arma::vec& x) { return arma::fft(arma::cx_vec(x, arma::zeros(x.n_elem))); }
inline arma::cx_vec fft(const arma::vec& x, arma::uword nfft) {
  arma::cx_vec xc(nfft, arma::fill::zeros);
  arma::uword n = std::min<arma::uword>(x.n_elem, nfft);
  for (arma::uword i = 0; i < n; ++i) xc(i) = x(i);
  return arma::fft(xc);
}
inline arma::cx_vec ifft(const arma::cx_vec& X) { return arma::ifft(X); }

// ---------- Misc ----------
inline arma::vec zero_phase_filter(const arma::vec& b, const arma::vec& a, const arma::vec& x) {
  return filtfilt(b, a, x);
}

inline arma::vec cumtrapz(const arma::vec& y, double dt = 1.0) {
  arma::uword N = y.n_elem;
  arma::vec z(N, arma::fill::zeros);
  for (arma::uword i = 1; i < N; ++i)
    z(i) = z(i - 1) + 0.5 * (y(i - 1) + y(i)) * dt;
  return z;
}

}  // namespace vat
#endif
