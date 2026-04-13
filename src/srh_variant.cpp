#include <RcppArmadillo.h>

namespace {

constexpr double kMatlabEps = 2.220446049250313e-16;

struct CoreResult {
  Rcpp::NumericMatrix f0_candidate;
  Rcpp::NumericVector srh_val;
};

double finite_abs_sum(const std::vector<double>& x) {
  double out = 0.0;
  for (double value : x) {
    if (std::isfinite(value)) {
      out += std::abs(value);
    }
  }
  return out;
}

double max_finite_value(const std::vector<double>& x) {
  double best = -std::numeric_limits<double>::infinity();
  for (double value : x) {
    if (std::isfinite(value) && value > best) {
      best = value;
    }
  }
  return std::isfinite(best) ? best : 0.0;
}

double mean_stddev(const Rcpp::NumericVector& x) {
  const int n = x.size();
  if (n <= 1) {
    return 0.0;
  }

  double mean = 0.0;
  for (double value : x) {
    mean += value;
  }
  mean /= static_cast<double>(n);

  double accum = 0.0;
  for (double value : x) {
    const double delta = value - mean;
    accum += delta * delta;
  }
  return std::sqrt(accum / static_cast<double>(n - 1));
}

double vector_median(std::vector<double> x) {
  if (x.empty()) {
    return NA_REAL;
  }

  std::sort(x.begin(), x.end());
  const std::size_t n = x.size();
  if (n % 2 == 1) {
    return x[n / 2];
  }
  return 0.5 * (x[n / 2 - 1] + x[n / 2]);
}

Rcpp::NumericVector symmetric_hann(const int n) {
  Rcpp::NumericVector out(n);
  if (n == 1) {
    out[0] = 1.0;
    return out;
  }

  const double denom = static_cast<double>(n - 1);
  for (int i = 0; i < n; ++i) {
    out[i] = 0.5 - 0.5 * std::cos((2.0 * M_PI * static_cast<double>(i)) / denom);
  }
  return out;
}

Rcpp::NumericVector symmetric_blackman(const int n) {
  Rcpp::NumericVector out(n);
  if (n == 1) {
    out[0] = 1.0;
    return out;
  }

  const double denom = static_cast<double>(n - 1);
  for (int i = 0; i < n; ++i) {
    const double angle = (2.0 * M_PI * static_cast<double>(i)) / denom;
    out[i] = 0.42 - 0.5 * std::cos(angle) + 0.08 * std::cos(2.0 * angle);
  }
  return out;
}

std::vector<double> autocorrelation(const std::vector<double>& x, const int order) {
  const int n = static_cast<int>(x.size());
  std::vector<double> r(order + 1, 0.0);

  for (int lag = 0; lag <= order; ++lag) {
    double accum = 0.0;
    for (int i = lag; i < n; ++i) {
      accum += x[i] * x[i - lag];
    }
    r[lag] = accum;
  }

  return r;
}

std::vector<double> levinson_durbin(const std::vector<double>& r, const int order) {
  std::vector<double> a(order + 1, 0.0);
  a[0] = 1.0;

  if (r.empty() || r[0] <= 0.0) {
    return a;
  }

  double e = r[0];

  for (int i = 1; i <= order; ++i) {
    double acc = r[i];
    for (int j = 1; j < i; ++j) {
      acc += a[j] * r[i - j];
    }

    const double k = -acc / (e + kMatlabEps);
    std::vector<double> updated = a;
    for (int j = 1; j < i; ++j) {
      updated[j] = a[j] + k * a[i - j];
    }
    updated[i] = k;
    a.swap(updated);

    e *= (1.0 - k * k);
    if (e <= kMatlabEps) {
      e = kMatlabEps;
    }
  }

  return a;
}

Rcpp::NumericVector fir_filter(const std::vector<double>& b, const Rcpp::NumericVector& x) {
  const int n = x.size();
  const int m = static_cast<int>(b.size());
  Rcpp::NumericVector y(n);

  for (int i = 0; i < n; ++i) {
    double accum = 0.0;
    for (int k = 0; k < m; ++k) {
      const int idx = i - k;
      if (idx < 0) {
        break;
      }
      accum += b[k] * x[idx];
    }
    y[i] = accum;
  }

  return y;
}

Rcpp::NumericVector cal_lpc(const Rcpp::NumericVector& wave,
                           const int frame_length,
                           const int hop_size,
                           const int lpc_order) {
  const int wave_len = wave.size();
  Rcpp::NumericVector hann_win = symmetric_hann(frame_length + 1);
  Rcpp::NumericVector residual(wave_len);

  int start = 0;
  int stop = start + frame_length;

  while (stop < wave_len) {
    Rcpp::NumericVector sig_frame(frame_length + 1);
    for (int i = 0; i <= frame_length; ++i) {
      sig_frame[i] = wave[start + i] * hann_win[i];
    }

    std::vector<double> frame_vec(sig_frame.begin(), sig_frame.end());
    const std::vector<double> acf = autocorrelation(frame_vec, lpc_order);
    const std::vector<double> a = levinson_durbin(acf, lpc_order);

    Rcpp::NumericVector inv = fir_filter(a, sig_frame);

    double sig_energy = 0.0;
    double inv_energy = 0.0;
    for (int i = 0; i <= frame_length; ++i) {
      sig_energy += sig_frame[i] * sig_frame[i];
      inv_energy += inv[i] * inv[i];
    }
    const double scale = std::sqrt(sig_energy / (inv_energy + kMatlabEps));

    for (int i = 0; i <= frame_length; ++i) {
      residual[start + i] += inv[i] * scale;
    }

    start += hop_size;
    stop += hop_size;
  }

  double max_abs = 0.0;
  for (double value : residual) {
    max_abs = std::max(max_abs, std::abs(value));
  }

  if (max_abs > 0.0) {
    residual = residual / max_abs;
  }

  return residual;
}

std::vector<int> get_peak_centers(const std::vector<double>& x,
                                  const int mpd,
                                  const double mhd,
                                  const int edge_low,
                                  const int edge_high,
                                  const int n_peak) {
  std::vector<double> pxx = x;
  const int n = static_cast<int>(pxx.size());
  const double mph = 0.0001 * max_finite_value(pxx);
  const double nan = std::numeric_limits<double>::quiet_NaN();
  std::vector<int> centers;
  std::vector<double> peaks;

  auto pos = [&](const int idx_1_based) -> double {
    return pxx[idx_1_based - 1];
  };

  auto is_nan_pos = [&](const int idx_1_based) -> bool {
    return !std::isfinite(pos(idx_1_based));
  };

  auto set_nan_range = [&](int left_1_based, int right_1_based) {
    left_1_based = std::max(1, left_1_based);
    right_1_based = std::min(n, right_1_based);
    for (int idx = left_1_based; idx <= right_1_based; ++idx) {
      pxx[idx - 1] = nan;
    }
  };

  if (n >= edge_low + 1 && pos(edge_low) >= pos(edge_low + 1)) {
    int idx_right = edge_low;
    while (idx_right <= n - 1 && pos(idx_right) >= pos(idx_right + 1)) {
      ++idx_right;
    }
    --idx_right;
    set_nan_range(1, idx_right);
  }

  if (n >= 2 && pos(n) >= pos(n - 1)) {
    int idx_left = n;
    while (idx_left > 1 && pos(idx_left) >= pos(idx_left - 1)) {
      --idx_left;
    }
    ++idx_left;
    set_nan_range(idx_left, n);
  }

  while (static_cast<int>(centers.size()) < n_peak && finite_abs_sum(pxx) > 0.0) {
    int temp_loc = 1;
    double temp_peak = -std::numeric_limits<double>::infinity();
    for (int idx = 1; idx <= n; ++idx) {
      const double value = pos(idx);
      if (std::isfinite(value) && value > temp_peak) {
        temp_peak = value;
        temp_loc = idx;
      }
    }

    if (temp_loc == 1 || temp_loc == n) {
      set_nan_range(temp_loc, temp_loc);
      continue;
    }

    if (is_nan_pos(temp_loc - 1)) {
      int idx_right = temp_loc + 1;
      while (idx_right <= n && pos(idx_right - 1) >= pos(idx_right)) {
        ++idx_right;
      }
      --idx_right;
      set_nan_range(temp_loc, idx_right);
      continue;
    }

    if (is_nan_pos(temp_loc + 1)) {
      int idx_left = temp_loc - 1;
      while (idx_left > 0 && pos(idx_left) <= pos(idx_left + 1)) {
        --idx_left;
      }
      ++idx_left;
      set_nan_range(idx_left, temp_loc);
      continue;
    }

    if (temp_peak > mph && pos(temp_loc) >= pos(temp_loc - 1) && pos(temp_loc) >= pos(temp_loc + 1)) {
      int idx_left = temp_loc - 1;
      while (idx_left > 0 && pos(idx_left) <= pos(idx_left + 1)) {
        --idx_left;
      }
      ++idx_left;

      int idx_right = temp_loc + 1;
      while (idx_right <= n && pos(idx_right - 1) >= pos(idx_right)) {
        ++idx_right;
      }
      --idx_right;

      if (pos(temp_loc) > pos(idx_left) && pos(temp_loc) > pos(idx_right)) {
        bool accept = true;
        if (!centers.empty()) {
          int min_dist = std::numeric_limits<int>::max();
          double min_height = std::numeric_limits<double>::infinity();
          for (std::size_t i = 0; i < centers.size(); ++i) {
            min_dist = std::min(min_dist, std::abs(temp_loc - centers[i]));
            min_height = std::min(min_height, std::abs(temp_peak - peaks[i]));
          }
          accept = (min_dist > mpd) || (min_height > mhd);
        }

        if (accept) {
          centers.push_back(temp_loc);
          peaks.push_back(temp_peak);
        }
      }

      set_nan_range(idx_left, idx_right);
    } else {
      break;
    }
  }

  return centers;
}

CoreResult srh_core(const Rcpp::NumericVector& sig,
                    const int fs,
                    const int edge_low,
                    const int edge_high,
                    const int n_candidate) {
  const int stop = 2048;
  const int hop_size = 200;
  const int delta_f = 2;
  const int frame_num = (sig.size() >= stop) ? ((sig.size() - stop) / hop_size + 1) : 0;

  Rcpp::NumericMatrix f0_candidate(n_candidate, frame_num);
  Rcpp::NumericVector srh_val(frame_num);
  if (frame_num == 0) {
    return {f0_candidate, srh_val};
  }

  Rcpp::NumericVector black_win = symmetric_blackman(stop);
  int start = 0;

  for (int index = 0; index < frame_num; ++index) {
    arma::vec sig_frame(stop);
    double mean = 0.0;
    for (int i = 0; i < stop; ++i) {
      const double value = sig[start + i] * black_win[i];
      sig_frame[i] = value;
      mean += value;
    }
    mean /= static_cast<double>(stop);
    sig_frame -= mean;

    arma::cx_vec spec_c = arma::fft(sig_frame, fs);
    arma::vec spec = arma::abs(spec_c.head(fs / 2));
    const double norm = std::sqrt(arma::accu(arma::square(spec)));
    if (norm > 0.0) {
      spec /= norm;
    }

    std::vector<double> srhs(edge_high, 0.0);
    for (int freq = edge_low; freq <= edge_high; ++freq) {
      auto harmonic_max = [&](const int center_1_based) {
        double best = -std::numeric_limits<double>::infinity();
        for (int idx = center_1_based - delta_f; idx <= center_1_based + delta_f; ++idx) {
          const int cpp_idx = idx - 1;
          if (cpp_idx >= 0 && cpp_idx < spec.n_elem) {
            best = std::max(best, spec[cpp_idx]);
          }
        }
        return std::isfinite(best) ? best : 0.0;
      };

      const double harmonic_sum =
        harmonic_max(1 * freq) +
        harmonic_max(2 * freq) +
        harmonic_max(3 * freq) +
        harmonic_max(4 * freq) +
        harmonic_max(5 * freq);

      const int idx15 = std::max(1, static_cast<int>(std::lround(1.5 * static_cast<double>(freq))));
      const int idx25 = std::max(1, static_cast<int>(std::lround(2.5 * static_cast<double>(freq))));
      const int idx35 = std::max(1, static_cast<int>(std::lround(3.5 * static_cast<double>(freq))));
      const int idx45 = std::max(1, static_cast<int>(std::lround(4.5 * static_cast<double>(freq))));

      const double anti_harmonic =
        spec[idx15 - 1] +
        spec[idx25 - 1] +
        spec[idx35 - 1] +
        spec[idx45 - 1];

      srhs[freq - 1] = harmonic_sum - anti_harmonic;
    }

    const double min_srh = *std::min_element(srhs.begin(), srhs.end());
    for (int idx = 0; idx < edge_low - 1; ++idx) {
      srhs[idx] = min_srh;
    }
    for (double& value : srhs) {
      value = value - min_srh + kMatlabEps;
    }

    const int peak_span = std::min(static_cast<int>(srhs.size()), 400);
    std::vector<double> srhs_peaks(srhs.begin(), srhs.begin() + peak_span);
    std::vector<int> peak_centers = get_peak_centers(srhs_peaks, 1, 0.0, edge_low, edge_high, n_candidate);

    if (peak_centers.empty()) {
      peak_centers.push_back(edge_low);
    }
    if (n_candidate == 2 && peak_centers.size() == 1) {
      peak_centers.push_back(peak_centers.front());
    }

    for (int row = 0; row < n_candidate; ++row) {
      const int center = peak_centers[std::min(row, static_cast<int>(peak_centers.size()) - 1)];
      f0_candidate(row, index) = center;
    }
    srh_val[index] = srhs[peak_centers.front() - 1];

    start += hop_size;
  }

  return {f0_candidate, srh_val};
}

Rcpp::NumericVector movmedian3(const Rcpp::NumericVector& x) {
  const int n = x.size();
  Rcpp::NumericVector out(n);

  for (int i = 0; i < n; ++i) {
    std::vector<double> window;
    for (int j = std::max(0, i - 1); j <= std::min(n - 1, i + 1); ++j) {
      window.push_back(x[j]);
    }
    out[i] = vector_median(window);
  }

  return out;
}

Rcpp::List srh_variant_impl(const Rcpp::NumericVector& wave,
                            const int fs,
                            const Rcpp::IntegerVector& edge) {
  const double delta_threshold = 0.11;
  const int length_threshold = 5;
  const int lpc_order = static_cast<int>(std::llround((3.0 / 4.0) * static_cast<double>(fs) / 1000.0));
  const int frame_length = static_cast<int>(std::llround((25.0 / 1000.0) * static_cast<double>(fs)));
  const int hop_size = static_cast<int>(std::llround((5.0 / 1000.0) * static_cast<double>(fs)));

  Rcpp::NumericVector residual = cal_lpc(wave, frame_length, hop_size, lpc_order);

  CoreResult pass1 = srh_core(residual, fs, edge[0], edge[1], 1);

  std::vector<double> voiced_candidates;
  for (int i = 0; i < pass1.srh_val.size(); ++i) {
    if (pass1.srh_val[i] > 0.1) {
      voiced_candidates.push_back(pass1.f0_candidate(0, i));
    }
  }
  const double f0_median = vector_median(voiced_candidates);

  Rcpp::IntegerVector edge_refined = Rcpp::clone(edge);
  if (!Rcpp::NumericVector::is_na(f0_median)) {
    edge_refined[0] = std::max(static_cast<int>(std::llround(0.5 * f0_median)), edge_refined[0]);
    edge_refined[1] = std::min(static_cast<int>(std::llround(2.0 * f0_median)), edge_refined[1]);
  }

  CoreResult pass2 = srh_core(residual, fs, edge_refined[0], edge_refined[1], 2);

  const int f0_num = pass2.f0_candidate.ncol();
  Rcpp::NumericVector f0(f0_num);
  for (int i = 0; i < f0_num; ++i) {
    f0[i] = pass2.f0_candidate(0, i);
  }
  for (int i = 2; i < f0_num - 1; ++i) {
    if (f0[i - 1] == f0[i - 2] &&
        f0[i] == f0[i - 1] + 1.0 &&
        std::abs(f0[i + 1] - f0[i]) > 20.0) {
      f0[i] = f0[i - 1];
    }
  }
  Rcpp::NumericVector f0_initial = Rcpp::clone(f0);

  int idx_stop = 1;
  int ii = 1;
  while (ii < f0_num - 1) {
    int idx_start = ii;
    int idx_end = ii + 1;
    while (idx_end <= f0_num &&
           (std::abs(f0[idx_end - 1] - f0[idx_end - 2]) / (f0[idx_end - 2] + kMatlabEps) < delta_threshold)) {
      ++idx_end;
    }
    --idx_end;

    if (idx_end - idx_start > length_threshold) {
      idx_start = idx_start - 1;
      while (idx_start > idx_stop) {
        bool any_match = false;
        double best_shift = std::numeric_limits<double>::infinity();
        int best_row = 0;
        for (int row = 0; row < pass2.f0_candidate.nrow(); ++row) {
          const double candidate = pass2.f0_candidate(row, idx_start - 1);
          const double ratio = std::abs(candidate - f0[idx_start]) / (f0[idx_start] + kMatlabEps);
          if (ratio < delta_threshold) {
            any_match = true;
          }
          const double shift = std::abs(candidate - f0[idx_start]);
          if (shift < best_shift) {
            best_shift = shift;
            best_row = row;
          }
        }
        if (!any_match) {
          break;
        }
        f0[idx_start - 1] = pass2.f0_candidate(best_row, idx_start - 1);
        --idx_start;
      }
      idx_start = idx_start + 1;

      idx_end = idx_end + 1;
      while (idx_end <= f0_num - 1) {
        bool any_match = false;
        double best_shift = std::numeric_limits<double>::infinity();
        int best_row = 0;
        for (int row = 0; row < pass2.f0_candidate.nrow(); ++row) {
          const double candidate = pass2.f0_candidate(row, idx_end - 1);
          const double ratio = std::abs(candidate - f0[idx_end - 2]) / (f0[idx_end - 2] + kMatlabEps);
          if (ratio < delta_threshold) {
            any_match = true;
          }
          const double shift = std::abs(candidate - f0[idx_end - 2]);
          if (shift < best_shift) {
            best_shift = shift;
            best_row = row;
          }
        }
        if (!any_match) {
          break;
        }
        f0[idx_end - 1] = pass2.f0_candidate(best_row, idx_end - 1);
        ++idx_end;
      }
      idx_end = idx_end - 1;
      idx_stop = idx_end;
    }

    ii = idx_end + 1;
  }

  Rcpp::NumericVector f0_pre_smooth = Rcpp::clone(f0);
  f0 = movmedian3(f0);

  Rcpp::IntegerVector vad(f0_num);
  double vad_threshold = 0.07;
  if (mean_stddev(pass2.srh_val) > 0.05) {
    vad_threshold = 0.085;
  }
  for (int i = 0; i < f0_num; ++i) {
    if (pass2.srh_val[i] > vad_threshold) {
      vad[i] = 1;
    }
  }

  ii = 1;
  while (ii < f0_num - 1) {
    const int idx_start = ii;
    int idx_end = ii + 1;
    while (idx_end <= f0_num && vad[idx_end - 1]) {
      ++idx_end;
    }
    --idx_end;

    if (idx_end - idx_start <= length_threshold) {
      for (int idx = idx_start; idx <= idx_end; ++idx) {
        vad[idx - 1] = 0;
      }
    }

    ii = idx_end + 1;
  }

  return Rcpp::List::create(
    Rcpp::Named("f0") = f0,
    Rcpp::Named("vad") = vad,
    Rcpp::Named("residual_spectral") = residual,
    Rcpp::Named("f0_initial") = f0_initial,
    Rcpp::Named("f0_pre_smooth") = f0_pre_smooth,
    Rcpp::Named("f0_candidate_pass1") = pass1.f0_candidate,
    Rcpp::Named("srh_val_pass1") = pass1.srh_val,
    Rcpp::Named("f0_candidate_pass2") = pass2.f0_candidate,
    Rcpp::Named("srh_val_pass2") = pass2.srh_val,
    Rcpp::Named("edge_refined") = edge_refined
  );
}

}  // namespace

// [[Rcpp::export]]
Rcpp::NumericVector resample_polyphase_cpp(const Rcpp::NumericVector& x,
                                           const Rcpp::NumericVector& h,
                                           const int up,
                                           const int down,
                                           const int n_out,
                                           const int trim) {
  const int x_len = x.size();
  const int h_len = h.size();
  Rcpp::NumericVector y(n_out);

  for (int j = 0; j < n_out; ++j) {
    const long long m = static_cast<long long>(j + trim) * static_cast<long long>(down);
    const long long min_k_raw = m - static_cast<long long>(h_len - 1);
    long long k_min = 0;
    if (min_k_raw > 0) {
      k_min = (min_k_raw + up - 1) / up;
    }
    long long k_max = std::min<long long>(x_len - 1, m / up);

    double accum = 0.0;
    for (long long k = k_min; k <= k_max; ++k) {
      const int h_idx = static_cast<int>(m - k * up);
      if (h_idx >= 0 && h_idx < h_len) {
        accum += x[static_cast<int>(k)] * h[h_idx];
      }
    }
    y[j] = accum;
  }

  return y;
}

// [[Rcpp::export]]
Rcpp::List srh_core_cpp(const Rcpp::NumericVector& sig,
                        const int fs,
                        const Rcpp::IntegerVector& edge,
                        const int n_candidate) {
  CoreResult out = srh_core(sig, fs, edge[0], edge[1], n_candidate);
  return Rcpp::List::create(
    Rcpp::Named("f0_candidate") = out.f0_candidate,
    Rcpp::Named("srh_val") = out.srh_val
  );
}

// [[Rcpp::export]]
Rcpp::List srh_variant_cpp(const Rcpp::NumericVector& wave,
                           const int fs,
                           const Rcpp::IntegerVector& edge) {
  Rcpp::List out = srh_variant_impl(wave, fs, edge);
  return Rcpp::List::create(
    Rcpp::Named("f0") = out["f0"],
    Rcpp::Named("vad") = out["vad"]
  );
}

// [[Rcpp::export]]
Rcpp::List srh_variant_debug_cpp(const Rcpp::NumericVector& wave,
                                 const int fs,
                                 const Rcpp::IntegerVector& edge) {
  return srh_variant_impl(wave, fs, edge);
}
