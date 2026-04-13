#include <RcppArmadillo.h>

using namespace Rcpp;

namespace {

arma::vec levinson_durbin(const arma::vec& x, const int order) {
  arma::vec r(order + 1, arma::fill::zeros);
  for (int lag = 0; lag <= order; ++lag) {
    double value = 0.0;
    for (arma::uword i = 0; i + lag < x.n_elem; ++i) {
      value += x[i] * x[i + lag];
    }
    r[lag] = value;
  }

  arma::vec a(order + 1, arma::fill::zeros);
  a[0] = 1.0;

  double e = r[0];
  if (e == 0.0) {
    return a;
  }

  for (int i = 1; i <= order; ++i) {
    double lambda = r[i];
    for (int j = 1; j < i; ++j) {
      lambda -= a[j] * r[i - j];
    }

    lambda /= e;

    arma::vec a_new = a;
    a_new[i] = lambda;
    for (int j = 1; j < i; ++j) {
      a_new[j] = a[j] - lambda * a[i - j];
    }

    a = a_new;
    e *= (1.0 - lambda * lambda);
    if (e <= 0.0) {
      break;
    }
  }

  return a;
}

arma::vec apply_fir(const arma::vec& b, const arma::vec& x) {
  arma::vec y(x.n_elem, arma::fill::zeros);
  for (arma::uword n = 0; n < x.n_elem; ++n) {
    double acc = 0.0;
    for (arma::uword k = 0; k < b.n_elem && k <= n; ++k) {
      acc += b[k] * x[n - k];
    }
    y[n] = acc;
  }
  return y;
}

arma::vec hann_window(const int n) {
  arma::vec out(n, arma::fill::zeros);
  if (n <= 1) {
    if (n == 1) {
      out[0] = 1.0;
    }
    return out;
  }

  for (int i = 0; i < n; ++i) {
    out[i] = 0.5 - 0.5 * std::cos((2.0 * M_PI * i) / static_cast<double>(n - 1));
  }
  return out;
}

arma::mat reshape_column_major(const arma::vec& x, const int rows, const int cols) {
  arma::mat out(rows, cols, arma::fill::zeros);
  for (int col = 0; col < cols; ++col) {
    for (int row = 0; row < rows; ++row) {
      out(row, col) = x[col * rows + row];
    }
  }
  return out;
}

}  // namespace

// [[Rcpp::export]]
arma::mat pre_emphasis_cpp(const arma::vec& x, const double preemp) {
  arma::mat out(x.n_elem, 1, arma::fill::zeros);
  if (x.n_elem == 0) {
    return out;
  }

  out(0, 0) = x[0];
  for (arma::uword i = 1; i < x.n_elem; ++i) {
    out(i, 0) = x[i] - preemp * x[i - 1];
  }
  return out;
}

// [[Rcpp::export]]
arma::mat tvlp_l2_cpp(const arma::vec& x, const int p, const int q) {
  const int n_samples = static_cast<int>(x.n_elem);
  const int n_equations = n_samples - p;
  const int n_coeffs = p * (q + 1);

  arma::mat ypu(n_equations, n_coeffs, arma::fill::zeros);
  arma::vec yn(n_equations, arma::fill::zeros);

  int m = 0;
  for (int n = p; n < n_samples; ++n) {
    yn[m] = x[n];

    for (int i = 1; i <= p; ++i) {
      for (int j = 0; j <= q; ++j) {
        const int col_idx = (i - 1) * (q + 1) + j;
        ypu(m, col_idx) = std::pow(static_cast<double>(n - p), j) * x[n - i];
      }
    }
    ++m;
  }

  arma::vec x_l2;
  arma::solve(x_l2, ypu, yn);
  return reshape_column_major(x_l2, q + 1, p);
}

// [[Rcpp::export]]
arma::mat tvwlp_l2_cpp(const arma::vec& x, const int p, const int q, const arma::vec& w) {
  const int n_samples = static_cast<int>(x.n_elem);
  const int n_equations = n_samples - p;
  const int n_coeffs = p * (q + 1);

  arma::mat ypu(n_equations, n_coeffs, arma::fill::zeros);
  arma::vec yn(n_equations, arma::fill::zeros);

  int m = 0;
  for (int n = p; n < n_samples; ++n) {
    yn[m] = w[n] * x[n];

    for (int i = 1; i <= p; ++i) {
      for (int j = 0; j <= q; ++j) {
        const int col_idx = (i - 1) * (q + 1) + j;
        ypu(m, col_idx) = w[n] * std::pow(static_cast<double>(n - p), j) * x[n - i];
      }
    }
    ++m;
  }

  arma::vec x_l2;
  arma::solve(x_l2, ypu, yn);
  return reshape_column_major(x_l2, q + 1, p);
}

// [[Rcpp::export]]
List tvlptoformants_akitofi_cpp(const arma::mat& aki, const int nx, const int npeaks, const double fs) {
  const int q = static_cast<int>(aki.n_rows) - 1;
  const int p = static_cast<int>(aki.n_cols);

  arma::rowvec tn = arma::linspace<arma::rowvec>(0.0, static_cast<double>(nx - 1), nx);
  arma::mat akn(p, nx, arma::fill::zeros);

  for (int k = 0; k < p; ++k) {
    for (int i = 0; i <= q; ++i) {
      akn.row(k) += aki(i, k) * arma::pow(tn, i);
    }
  }

  arma::mat ak(p + 1, nx, arma::fill::zeros);
  ak.row(0).ones();
  ak.rows(1, p) = -akn;

  arma::mat fi(nx, npeaks, arma::fill::zeros);
  arma::mat bw(nx, npeaks, arma::fill::zeros);
  for (int idx = 0; idx < nx; ++idx) {
    arma::cx_vec roots_i = arma::roots(ak.col(idx));
    std::vector<std::pair<double, double>> freq_bw;
    freq_bw.reserve(roots_i.n_elem);

    for (arma::uword root_idx = 0; root_idx < roots_i.n_elem; ++root_idx) {
      const double freq = std::arg(roots_i[root_idx]) * fs / (2.0 * M_PI);
      if (freq > 0.0) {
        const double bandwidth = -std::log(std::abs(roots_i[root_idx])) * fs / M_PI;
        freq_bw.push_back(std::make_pair(freq, bandwidth));
      }
    }

    std::sort(freq_bw.begin(), freq_bw.end());
    const int n_available = std::min(npeaks, static_cast<int>(freq_bw.size()));
    for (int j = 0; j < n_available; ++j) {
      fi(idx, j) = freq_bw[j].first;
      bw(idx, j) = freq_bw[j].second;
    }
  }

  return List::create(
    _["fi"] = fi,
    _["bw"] = bw,
    _["ak"] = ak
  );
}

// [[Rcpp::export]]
arma::vec get_lpc_residual_cpp(const arma::vec& wave, const int l, const int shift, const int order) {
  const arma::uword n_samples = wave.n_elem;
  arma::vec res(n_samples, arma::fill::zeros);
  const arma::vec hann_win = hann_window(l + 1);

  int start = 0;
  int stop = start + l;

  while (stop < static_cast<int>(n_samples)) {
    arma::vec segment = wave.subvec(start, stop);
    segment %= hann_win;

    const arma::vec a = levinson_durbin(segment, order);
    arma::vec inv = apply_fir(a, segment);

    const double seg_energy = arma::dot(segment, segment);
    const double inv_energy = arma::dot(inv, inv);
    if (inv_energy > 0.0) {
      inv *= std::sqrt(seg_energy / inv_energy);
    }

    res.subvec(start, stop) += inv;
    start += shift;
    stop += shift;
  }

  const double res_max = arma::abs(res).max();
  if (res_max > 0.0) {
    res /= res_max;
  }

  return res;
}
