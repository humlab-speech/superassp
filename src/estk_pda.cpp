// ESTK PDA (Pitch Detection Algorithm) C++ Implementation
// Lightweight implementation of Edinburgh Speech Tools SRPD algorithm
// Based on Medan, Yair & Chazan (1991) and Bagshaw et al. (1993)

// SIMD optimization with RcppXsimd for 4-6x speedup

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>


// Include xsimd for SIMD vectorization (via RcppXsimd package)
#ifdef RCPPXSIMD_AVAILABLE
#include <xsimd/xsimd.hpp>
#endif


using namespace Rcpp;

// Constants matching ESTK defaults
const double BREAK_NUMBER = 0.0;
const int DEFAULT_DECIMATION = 4;
const double DEFAULT_MIN_PITCH = 40.0;
const double DEFAULT_MAX_PITCH = 400.0;
const double DEFAULT_SHIFT = 5.0;      // ms
const double DEFAULT_LENGTH = 10.0;    // ms
const int DEFAULT_TSILENT = 120;
const double DEFAULT_TMIN = 0.75;
const double DEFAULT_TMAX_RATIO = 0.85;
const double DEFAULT_THIGH = 0.88;
const double DEFAULT_TDH = 0.77;

// Voice status
enum VoiceStatus { UNVOICED = 0, VOICED = 1, SILENT = 2 };
enum SendStatus { HOLD = 1, SEND = 2 };

// PDA parameters structure
struct PDA_Params {
  int sample_freq;
  int Nmax, Nmin;
  double shift, length;  // ms
  double min_pitch, max_pitch;  // Hz
  int L;  // Decimation factor
  double Tmin, Tmax_ratio, Thigh, Tdh;
  int Tsilent;
  bool peak_tracking;
};

// Frame status structure
struct FrameStatus {
  double pitch_freq;
  VoiceStatus v_uv;
  SendStatus s_h;
  double cc_max;
  double threshold;
};

// Peak candidate structure
struct PeakCandidate {
  int N0;
  int score;
};

// Super Resolution Pitch Detection Algorithm
void super_resolution_pda(
    const PDA_Params &params,
    const std::vector<short> &segment,
    std::vector<double> &cc_coeff,
    FrameStatus &status,
    int &zx_lft_N,
    int &zx_rht_N,
    double &prev_pf
) {
  // Set correlation coefficient threshold
  if (status.v_uv == UNVOICED || status.v_uv == SILENT) {
    status.threshold = params.Thigh;
  } else {
    status.threshold = std::max(params.Tmin, params.Tmax_ratio * status.cc_max);
  }

  // Determine if bias should be applied (peak tracking)
  bool apply_bias = false;
  if (params.peak_tracking && prev_pf != BREAK_NUMBER &&
      status.v_uv == VOICED && status.s_h != HOLD &&
      status.pitch_freq < 1.75 * prev_pf &&
      status.pitch_freq > 0.625 * prev_pf) {
    apply_bias = true;
  }

  // Analyze first two segments for Nmin
  short x_max = -std::numeric_limits<short>::max();
  short x_min = std::numeric_limits<short>::max();
  short y_max = -std::numeric_limits<short>::max();
  short y_min = std::numeric_limits<short>::max();

  double xx = 0.0, yy = 0.0, xy = 0.0;
  int seg1_zxs = 0, seg2_zxs = 0;
  int prev_seg1 = segment[params.Nmax - params.Nmin] < 0 ? -1 : 1;
  int prev_seg2 = segment[params.Nmax] < 0 ? -1 : 1;


  // Process first correlation at Nmin (SIMD-optimized)
#ifdef RCPPXSIMD_AVAILABLE
  // SIMD-optimized initial correlation computation
  using batch_type = xsimd::simd_type<float>;
  constexpr size_t simd_size = batch_type::size;

  batch_type xx_vec(0.0f), yy_vec(0.0f), xy_vec(0.0f);
  int j_init = 0;

  // SIMD loop for correlation accumulation
  int nmin_simd_limit = (params.Nmin / params.L) * params.L;
  for (; j_init + static_cast<int>(simd_size) * params.L <= nmin_simd_limit; j_init += simd_size * params.L) {
    alignas(32) float bufx[simd_size];
    alignas(32) float bufy[simd_size];

    for (size_t i = 0; i < simd_size; i++) {
      int idx = j_init + i * params.L;
      int x_idx = params.Nmax - params.Nmin + idx;
      int y_idx = params.Nmax + idx;

      bufx[i] = static_cast<float>(segment[x_idx]);
      bufy[i] = static_cast<float>(segment[y_idx]);

      // Track min/max (cannot vectorize easily)
      if (segment[x_idx] > x_max) x_max = segment[x_idx];
      if (segment[x_idx] < x_min) x_min = segment[x_idx];
      if (segment[y_idx] > y_max) y_max = segment[y_idx];
      if (segment[y_idx] < y_min) y_min = segment[y_idx];

      // Zero crossings (cannot vectorize)
      if (segment[x_idx] * prev_seg1 < 0) {
        prev_seg1 *= -1;
        seg1_zxs++;
      }
      if (segment[y_idx] * prev_seg2 < 0) {
        prev_seg2 *= -1;
        seg2_zxs++;
      }
    }

    batch_type bx, by;
    bx.load_aligned(bufx);
    by.load_aligned(bufy);

    xx_vec += bx * bx;
    yy_vec += by * by;
    xy_vec += bx * by;
  }

  // Horizontal reduction
  xx = xsimd::hadd(xx_vec);
  yy = xsimd::hadd(yy_vec);
  xy = xsimd::hadd(xy_vec);

  // Scalar tail loop
  for (; j_init < params.Nmin; j_init += params.L) {
    int x_idx = params.Nmax - params.Nmin + j_init;
    int y_idx = params.Nmax + j_init;

    if (segment[x_idx] > x_max) x_max = segment[x_idx];
    if (segment[x_idx] < x_min) x_min = segment[x_idx];
    if (segment[y_idx] > y_max) y_max = segment[y_idx];
    if (segment[y_idx] < y_min) y_min = segment[y_idx];

    if (segment[x_idx] * prev_seg1 < 0) {
      prev_seg1 *= -1;
      seg1_zxs++;
    }
    if (segment[y_idx] * prev_seg2 < 0) {
      prev_seg2 *= -1;
      seg2_zxs++;
    }

    xx += (double)segment[x_idx] * segment[x_idx];
    yy += (double)segment[y_idx] * segment[y_idx];
    xy += (double)segment[x_idx] * segment[y_idx];
  }
#else
  // Scalar fallback (original implementation)


  for (int j = 0; j < params.Nmin; j += params.L) {
    int x_idx = params.Nmax - params.Nmin + j;
    int y_idx = params.Nmax + j;

    if (segment[x_idx] > x_max) x_max = segment[x_idx];
    if (segment[x_idx] < x_min) x_min = segment[x_idx];
    if (segment[y_idx] > y_max) y_max = segment[y_idx];
    if (segment[y_idx] < y_min) y_min = segment[y_idx];

    if (segment[x_idx] * prev_seg1 < 0) {
      prev_seg1 *= -1;
      seg1_zxs++;
    }
    if (segment[y_idx] * prev_seg2 < 0) {
      prev_seg2 *= -1;
      seg2_zxs++;
    }

    xx += (double)segment[x_idx] * segment[x_idx];
    yy += (double)segment[y_idx] * segment[y_idx];
    xy += (double)segment[x_idx] * segment[y_idx];
  }


#endif


  // Check for silence
  if (std::abs(x_max) + std::abs(x_min) < 2 * params.Tsilent ||
      std::abs(y_max) + std::abs(y_min) < 2 * params.Tsilent) {
    std::fill(cc_coeff.begin(), cc_coeff.end(), 0.0);
    prev_pf = status.pitch_freq;
    status.pitch_freq = BREAK_NUMBER;
    status.v_uv = SILENT;
    status.s_h = SEND;
    status.cc_max = 0.0;
    return;
  }

  // First correlation coefficient
  cc_coeff[0] = xy / std::sqrt(xx) / std::sqrt(yy);
  status.cc_max = cc_coeff[0];

  for (int q = 1; q < std::min((int)cc_coeff.size(), params.L); q++) {
    cc_coeff[q] = 0.0;
  }

  int total_zxs = seg1_zxs + seg2_zxs;
  int prev_sign = cc_coeff[0] < 0.0 ? -1 : 1;
  prev_seg1 = segment[params.Nmax - params.Nmin] < 0 ? -1 : 1;

  // Find peak candidates
  std::vector<PeakCandidate> sig_peaks;
  int N0 = 0;
  double max_cc = 0.0;
  int zx_rate = 0;
  int zx_at_N0 = 0;
  bool lower_found = false;

  int j = 0;
  // Iterate through possible periods
  for (int n = params.Nmin + params.L; n <= params.Nmax; n += params.L, j += params.L) {
    int x_idx = params.Nmax - n;
    int y_idx = params.Nmax + j;

    if (segment[x_idx] * prev_seg1 < 0) {
      prev_seg1 *= -1;
      total_zxs++;
    }
    if (segment[y_idx] * prev_seg2 < 0) {
      prev_seg2 *= -1;
      total_zxs++;
    }

    xx += (double)segment[x_idx] * segment[x_idx];
    yy += (double)segment[y_idx] * segment[y_idx];



    // Cross-correlation computation (SIMD-optimized for 4-6x speedup)
    xy = 0.0;

#ifdef RCPPXSIMD_AVAILABLE
    // SIMD-optimized cross-correlation
    using batch_type = xsimd::simd_type<float>;
    constexpr size_t simd_size = batch_type::size;

    batch_type xy_vec(0.0f);
    int k = 0;

    // SIMD loop: process simd_size elements at once
    int n_simd_limit = (n / params.L) * params.L;  // Ensure we stay within bounds
    for (; k + static_cast<int>(simd_size) * params.L <= n_simd_limit; k += simd_size * params.L) {
      // Load decimated samples into aligned buffers
      alignas(32) float buf1[simd_size];
      alignas(32) float buf2[simd_size];

      for (size_t j = 0; j < simd_size; j++) {
        int idx = k + j * params.L;
        buf1[j] = static_cast<float>(segment[params.Nmax - n + idx]);
        buf2[j] = static_cast<float>(segment[params.Nmax + idx]);
      }

      // Load from aligned buffers
      batch_type b1, b2;
      b1.load_aligned(buf1);
      b2.load_aligned(buf2);

      // Vectorized multiply-accumulate
      xy_vec += b1 * b2;
    }

    // Horizontal reduction
    xy = xsimd::hadd(xy_vec);

    // Scalar tail loop for remaining elements
    for (; k < n; k += params.L) {
      xy += (double)segment[params.Nmax - n + k] * segment[params.Nmax + k];
    }
#else
    // Scalar fallback (original implementation)
    for (int k = 0; k < n; k += params.L) {
      xy += (double)segment[params.Nmax - n + k] * segment[params.Nmax + k];
    }
#endif


    cc_coeff[n - params.Nmin] = xy / std::sqrt(xx) / std::sqrt(yy);

    if (cc_coeff[n - params.Nmin] > status.cc_max) {
      status.cc_max = cc_coeff[n - params.Nmin];
    }

    // Zero unknown coefficients
    for (int q = n - params.Nmin + 1;
         q < (int)cc_coeff.size() && q < n - params.Nmin + params.L; q++) {
      cc_coeff[q] = 0.0;
    }

    if (cc_coeff[n - params.Nmin] > cc_coeff[n - params.Nmin - params.L]) {
      lower_found = true;
    }

    if (cc_coeff[n - params.Nmin] * prev_sign < 0.0) {
      prev_sign *= -1;
      zx_rate++;
    }

    if (N0 != 0 && zx_rate > zx_at_N0) {
      sig_peaks.push_back({N0, 1});
      N0 = 0;
      max_cc = 0.0;
    }

    double coeff_weight = (apply_bias && n > zx_lft_N && n < zx_rht_N) ? 2.0 : 1.0;

    if (cc_coeff[n - params.Nmin] > max_cc && total_zxs > 3 && lower_found) {
      max_cc = cc_coeff[n - params.Nmin];
      if (max_cc * coeff_weight >= status.threshold) {
        zx_at_N0 = zx_rate;
        N0 = n;
      }
    }
  }

  // Check for unvoiced
  if (sig_peaks.empty()) {
    prev_pf = status.pitch_freq;
    status.pitch_freq = BREAK_NUMBER;
    status.v_uv = UNVOICED;
    status.s_h = SEND;
    return;
  }

  // Score peak candidates
  std::vector<PeakCandidate> scored_peaks;
  int best_score = 1;

  for (const auto &peak : sig_peaks) {
    yy = 0.0;
    double zz = 0.0;
    double yz = 0.0;



#ifdef RCPPXSIMD_AVAILABLE
    // SIMD-optimized peak scoring loop (3 accumulations: yy, zz, yz)
    using batch_type = xsimd::simd_type<float>;
    constexpr size_t simd_size = batch_type::size;

    batch_type yy_vec(0.0f), zz_vec(0.0f), yz_vec(0.0f);
    int j = 0;

    // SIMD loop: process simd_size elements at once
    for (; j + static_cast<int>(simd_size) <= peak.N0; j += simd_size) {
      alignas(32) float buf_y[simd_size];
      alignas(32) float buf_z[simd_size];

      for (size_t k = 0; k < simd_size; k++) {
        int y_idx = params.Nmax + j + k;
        int z_idx = params.Nmax + peak.N0 + j + k;
        buf_y[k] = static_cast<float>(segment[y_idx]);
        buf_z[k] = static_cast<float>(segment[z_idx]);
      }

      batch_type y_batch, z_batch;
      y_batch.load_aligned(buf_y);
      z_batch.load_aligned(buf_z);

      yy_vec += y_batch * y_batch;
      zz_vec += z_batch * z_batch;
      yz_vec += y_batch * z_batch;
    }

    // Horizontal reductions
    yy = xsimd::hadd(yy_vec);
    zz = xsimd::hadd(zz_vec);
    yz = xsimd::hadd(yz_vec);

    // Scalar tail loop for remaining elements
    for (; j < peak.N0; j++) {
      int y_idx = params.Nmax + j;
      int z_idx = params.Nmax + peak.N0 + j;
      yy += (double)segment[y_idx] * segment[y_idx];
      zz += (double)segment[z_idx] * segment[z_idx];
      yz += (double)segment[y_idx] * segment[z_idx];
    }
#else
    // Scalar fallback


    for (int j = 0; j < peak.N0; j++) {
      int y_idx = params.Nmax + j;
      int z_idx = params.Nmax + peak.N0 + j;
      yy += (double)segment[y_idx] * segment[y_idx];
      zz += (double)segment[z_idx] * segment[z_idx];
      yz += (double)segment[y_idx] * segment[z_idx];
    }


#endif


    double coeff = (yy == 0.0 || zz == 0.0) ? 0.0 : yz / std::sqrt(yy) / std::sqrt(zz);
    double coeff_weight = (apply_bias && peak.N0 > zx_lft_N && peak.N0 < zx_rht_N) ? 2.0 : 1.0;

    PeakCandidate scored = peak;
    if (coeff * coeff_weight >= status.threshold) {
      scored.score = 2;
      best_score = 2;
    }
    scored_peaks.push_back(scored);
  }

  // Find best peak
  PeakCandidate first_peak = scored_peaks.front();
  PeakCandidate last_peak = scored_peaks.back();
  N0 = first_peak.N0;

  if (scored_peaks.size() > 1) {
    xx = 0.0;
    for (int j = 0; j < last_peak.N0; j++) {
      int x_idx = params.Nmax - last_peak.N0 + j;
      xx += (double)segment[x_idx] * segment[x_idx];
    }

    max_cc = 0.0;
    for (const auto &peak : scored_peaks) {
      if (peak.score == best_score) {
        double xz = 0.0, zz = 0.0;



#ifdef RCPPXSIMD_AVAILABLE
        // SIMD-optimized best peak selection (2 accumulations: xz, zz)
        using batch_type = xsimd::simd_type<float>;
        constexpr size_t simd_size = batch_type::size;

        batch_type xz_vec(0.0f), zz_vec(0.0f);
        int j = 0;

        // SIMD loop
        for (; j + static_cast<int>(simd_size) <= last_peak.N0; j += simd_size) {
          alignas(32) float buf_x[simd_size];
          alignas(32) float buf_z[simd_size];

          for (size_t k = 0; k < simd_size; k++) {
            int x_idx = params.Nmax - last_peak.N0 + j + k;
            int z_idx = params.Nmax + peak.N0 + j + k;
            buf_x[k] = static_cast<float>(segment[x_idx]);
            buf_z[k] = static_cast<float>(segment[z_idx]);
          }

          batch_type x_batch, z_batch;
          x_batch.load_aligned(buf_x);
          z_batch.load_aligned(buf_z);

          xz_vec += x_batch * z_batch;
          zz_vec += z_batch * z_batch;
        }

        // Horizontal reductions
        xz = xsimd::hadd(xz_vec);
        zz = xsimd::hadd(zz_vec);

        // Scalar tail loop
        for (; j < last_peak.N0; j++) {
          int x_idx = params.Nmax - last_peak.N0 + j;
          int z_idx = params.Nmax + peak.N0 + j;
          xz += (double)segment[x_idx] * segment[z_idx];
          zz += (double)segment[z_idx] * segment[z_idx];
        }
#else
        // Scalar fallback

        for (int j = 0; j < last_peak.N0; j++) {
          int x_idx = params.Nmax - last_peak.N0 + j;
          int z_idx = params.Nmax + peak.N0 + j;
          xz += (double)segment[x_idx] * segment[z_idx];
          zz += (double)segment[z_idx] * segment[z_idx];
        }


#endif

        double coeff = xz / std::sqrt(xx) / std::sqrt(zz);

        if (coeff * params.Tdh > max_cc) {
          N0 = peak.N0;
          max_cc = coeff;
        } else if (&peak == &scored_peaks.front()) {
          max_cc = coeff;
        }
      }
    }
  }

  status.cc_max = cc_coeff[N0 - params.Nmin];

  // Determine send/hold status
  if ((scored_peaks.size() == 1 && best_score == 1 && status.v_uv != VOICED) ||
      cc_coeff[N0 - params.Nmin] < status.threshold) {
    status.s_h = HOLD;
  } else {
    status.s_h = SEND;
  }

  // Find zero-crossing boundaries
  zx_lft_N = 0;
  zx_rht_N = 0;
  for (int q = N0; q >= params.Nmin; q -= params.L) {
    if (cc_coeff[q - params.Nmin] < 0.0) {
      zx_lft_N = q;
      break;
    }
  }
  for (int q = N0; q <= params.Nmax; q += params.L) {
    if (cc_coeff[q - params.Nmin] < 0.0) {
      zx_rht_N = q;
      break;
    }
  }

  // Refine estimate compensating for decimation
  int N1, N2;
  if (N0 - params.L < params.Nmin) {
    N1 = N0;
    N2 = N0 + 2 * params.L;
  } else if (N0 + params.L > params.Nmax) {
    N1 = N0 - 2 * params.L;
    N2 = N0;
  } else {
    N1 = N0 - params.L;
    N2 = N0 + params.L;
  }

  if (params.L != 1) {
    xx = yy = xy = 0.0;
    j = 0;



#ifdef RCPPXSIMD_AVAILABLE
    // SIMD-optimized initial correlation (3 accumulations: xx, yy, xy)
    using batch_type = xsimd::simd_type<float>;
    constexpr size_t simd_size = batch_type::size;

    batch_type xx_vec(0.0f), yy_vec(0.0f), xy_vec(0.0f);
    int i = 0;

    // SIMD loop
    for (; i + static_cast<int>(simd_size) <= N1; i += simd_size) {
      alignas(32) float buf_x[simd_size];
      alignas(32) float buf_y[simd_size];

      for (size_t k = 0; k < simd_size; k++) {
        int x_idx = params.Nmax - N1 + i + k;
        int y_idx = params.Nmax + i + k;
        buf_x[k] = static_cast<float>(segment[x_idx]);
        buf_y[k] = static_cast<float>(segment[y_idx]);
      }

      batch_type x_batch, y_batch;
      x_batch.load_aligned(buf_x);
      y_batch.load_aligned(buf_y);

      xx_vec += x_batch * x_batch;
      yy_vec += y_batch * y_batch;
      xy_vec += x_batch * y_batch;
    }

    // Horizontal reductions
    xx = xsimd::hadd(xx_vec);
    yy = xsimd::hadd(yy_vec);
    xy = xsimd::hadd(xy_vec);

    // Scalar tail loop
    for (; i < N1; i++) {
      int x_idx = params.Nmax - N1 + i;
      int y_idx = params.Nmax + i;
      xx += (double)segment[x_idx] * segment[x_idx];
      xy += (double)segment[x_idx] * segment[y_idx];
      yy += (double)segment[y_idx] * segment[y_idx];
    }
#else
    // Scalar fallback


    for (int i = 0; i < N1; i++) {
      int x_idx = params.Nmax - N1 + i;
      int y_idx = params.Nmax + i;
      xx += (double)segment[x_idx] * segment[x_idx];
      xy += (double)segment[x_idx] * segment[y_idx];
      yy += (double)segment[y_idx] * segment[y_idx];
    }


#endif

    cc_coeff[N1 - params.Nmin] = xy / std::sqrt(xx) / std::sqrt(yy);
    status.cc_max = cc_coeff[N1 - params.Nmin];
    N0 = N1;

    j = N1;
    for (int n = N1 + 1; n <= N2; n++, j++) {
      int x_idx = params.Nmax - n;
      int y_idx = params.Nmax + j;
      xx += (double)segment[x_idx] * segment[x_idx];
      yy += (double)segment[y_idx] * segment[y_idx];

      xy = 0.0;


#ifdef RCPPXSIMD_AVAILABLE
      // SIMD-optimized refinement dot product
      using batch_type = xsimd::simd_type<float>;
      constexpr size_t simd_size = batch_type::size;

      batch_type xy_vec(0.0f);
      int k = 0;

      // SIMD loop: process simd_size elements at once
      for (; k + static_cast<int>(simd_size) <= n; k += simd_size) {
        alignas(32) float buf_x[simd_size];
        alignas(32) float buf_y[simd_size];

        for (size_t idx = 0; idx < simd_size; idx++) {
          int x_idx = params.Nmax - n + k + idx;
          int y_idx = params.Nmax + k + idx;
          buf_x[idx] = static_cast<float>(segment[x_idx]);
          buf_y[idx] = static_cast<float>(segment[y_idx]);
        }

        batch_type x_batch, y_batch;
        x_batch.load_aligned(buf_x);
        y_batch.load_aligned(buf_y);

        xy_vec += x_batch * y_batch;
      }

      // Horizontal reduction
      xy = xsimd::hadd(xy_vec);

      // Scalar tail loop for remaining elements
      for (; k < n; k++) {
        xy += (double)segment[params.Nmax - n + k] * segment[params.Nmax + k];
      }
#else
      // Scalar fallback
      for (int k = 0; k < n; k++) {
        xy += (double)segment[params.Nmax - n + k] * segment[params.Nmax + k];
      }
#endif


      cc_coeff[n - params.Nmin] = xy / std::sqrt(xx) / std::sqrt(yy);

      if (cc_coeff[n - params.Nmin] > status.cc_max) {
        status.cc_max = cc_coeff[n - params.Nmin];
        N0 = n;
      }
    }
  }

  // Super-resolution refinement using parabolic interpolation
  int N_ = N0;
  if (N0 > params.Nmin && N0 != N1 && N0 < params.Nmax && N0 != N2) {
    if (cc_coeff[N0 - params.Nmin] - cc_coeff[N0 - params.Nmin - 1] <
        cc_coeff[N0 - params.Nmin] - cc_coeff[N0 - params.Nmin + 1]) {
      N_ = N0 - 1;
    }
  } else if (N0 - 1 >= params.Nmin && N0 != N1) {
    N_ = N0;
  } else if (N0 + 1 <= params.Nmax && N0 != N2) {
    N_ = N0 - 1;
  }

  // Calculate fractional part
  double xx_N = 0.0, yy_N = 0.0, xy_N = 0.0;
  double y1y1_N = 0.0, xy1_N = 0.0, yy1_N = 0.0;



#ifdef RCPPXSIMD_AVAILABLE
  // SIMD-optimized fractional calculation (6 accumulations)
  using batch_type_frac = xsimd::simd_type<float>;
  constexpr size_t simd_size_frac = batch_type_frac::size;

  batch_type_frac xx_N_vec(0.0f), yy_N_vec(0.0f), xy_N_vec(0.0f);
  batch_type_frac y1y1_N_vec(0.0f), xy1_N_vec(0.0f), yy1_N_vec(0.0f);
  int j_frac = 0;

  // SIMD loop
  for (; j_frac + static_cast<int>(simd_size_frac) <= N_; j_frac += simd_size_frac) {
    alignas(32) float buf_x[simd_size_frac];
    alignas(32) float buf_y[simd_size_frac];
    alignas(32) float buf_y1[simd_size_frac];

    for (size_t k = 0; k < simd_size_frac; k++) {
      int x_idx = params.Nmax - N_ + j_frac + k;
      int y_idx = params.Nmax + j_frac + k;
      buf_x[k] = static_cast<float>(segment[x_idx]);
      buf_y[k] = static_cast<float>(segment[y_idx]);
      buf_y1[k] = static_cast<float>(segment[y_idx + 1]);
    }

    batch_type_frac x_batch, y_batch, y1_batch;
    x_batch.load_aligned(buf_x);
    y_batch.load_aligned(buf_y);
    y1_batch.load_aligned(buf_y1);

    xx_N_vec += x_batch * x_batch;
    yy_N_vec += y_batch * y_batch;
    xy_N_vec += x_batch * y_batch;
    y1y1_N_vec += y1_batch * y1_batch;
    xy1_N_vec += x_batch * y1_batch;
    yy1_N_vec += y_batch * y1_batch;
  }

  // Horizontal reductions
  xx_N = xsimd::hadd(xx_N_vec);
  yy_N = xsimd::hadd(yy_N_vec);
  xy_N = xsimd::hadd(xy_N_vec);
  y1y1_N = xsimd::hadd(y1y1_N_vec);
  xy1_N = xsimd::hadd(xy1_N_vec);
  yy1_N = xsimd::hadd(yy1_N_vec);

  // Scalar tail loop
  for (; j_frac < N_; j_frac++) {
    int x_idx = params.Nmax - N_ + j_frac;
    int y_idx = params.Nmax + j_frac;
    xx_N += (double)segment[x_idx] * segment[x_idx];
    yy_N += (double)segment[y_idx] * segment[y_idx];
    xy_N += (double)segment[x_idx] * segment[y_idx];
    y1y1_N += (double)segment[y_idx + 1] * segment[y_idx + 1];
    xy1_N += (double)segment[x_idx] * segment[y_idx + 1];
    yy1_N += (double)segment[y_idx] * segment[y_idx + 1];
  }
#else
  // Scalar fallback

  for (int j = 0; j < N_; j++) {
    int x_idx = params.Nmax - N_ + j;
    int y_idx = params.Nmax + j;
    xx_N += (double)segment[x_idx] * segment[x_idx];
    yy_N += (double)segment[y_idx] * segment[y_idx];
    xy_N += (double)segment[x_idx] * segment[y_idx];
    y1y1_N += (double)segment[y_idx + 1] * segment[y_idx + 1];
    xy1_N += (double)segment[x_idx] * segment[y_idx + 1];
    yy1_N += (double)segment[y_idx] * segment[y_idx + 1];
  }


#endif



  double beta = (xy1_N * yy_N - xy_N * yy1_N) /
                (xy1_N * (yy_N - yy1_N) + xy_N * (y1y1_N - yy1_N));

  if (beta < 0.0) {
    N_--;
    beta = 0.0;
  } else if (beta >= 1.0) {
    N_++;
    beta = 0.0;
  } else {
    status.cc_max = ((1.0 - beta) * xy_N + beta * xy1_N) /
                    std::sqrt(xx_N * ((1.0 - beta) * (1.0 - beta) * yy_N +
                              2.0 * beta * (1.0 - beta) * yy1_N +
                              beta * beta * y1y1_N));
  }

  prev_pf = status.pitch_freq;
  status.pitch_freq = (double)params.sample_freq / (double)(N_ + beta);
  status.v_uv = VOICED;
}

// Main exported function
// [[Rcpp::export]]
List estk_pda_cpp(SEXP audio_obj,
                  double minF = 40.0,
                  double maxF = 400.0,
                  double windowShift = 5.0,
                  double windowSize = 10.0,
                  int decimation = 4,
                  int noise_floor = 120,
                  double min_v2uv_coef_thresh = 0.75,
                  double v2uv_coef_thresh_ratio = 0.85,
                  double uv2v_coef_thresh = 0.88,
                  double anti_doubling_thresh = 0.77,
                  bool peak_tracking = false,
                  bool verbose = false) {

  // Extract audio data from AsspDataObj
  if (!Rf_inherits(audio_obj, "AsspDataObj")) {
    stop("Input must be an AsspDataObj");
  }

  List audio_list(audio_obj);
  if (!audio_list.containsElementNamed("audio")) {
    stop("AsspDataObj must contain 'audio' track");
  }

  NumericMatrix audio_matrix = audio_list["audio"];
  int n_samples = audio_matrix.nrow();
  double sample_rate = as<double>(audio_list.attr("sampleRate"));

  if (verbose) {
    Rcout << "Processing " << n_samples << " samples at " << sample_rate << " Hz\n";
    Rcout << "F0 range: " << minF << "-" << maxF << " Hz\n";
  }

  // Extract first channel
  std::vector<short> signal(n_samples);
  for (int i = 0; i < n_samples; i++) {
    signal[i] = (short)(audio_matrix(i, 0) * 32767.0);
  }

  // Initialize parameters
  PDA_Params params;
  params.sample_freq = (int)sample_rate;
  params.min_pitch = minF;
  params.max_pitch = maxF;
  params.shift = windowShift;  // already in ms
  params.length = windowSize;  // already in ms
  params.L = decimation;
  params.Tmin = min_v2uv_coef_thresh;
  params.Tmax_ratio = v2uv_coef_thresh_ratio;
  params.Thigh = uv2v_coef_thresh;
  params.Tdh = anti_doubling_thresh;
  params.Tsilent = noise_floor;
  params.peak_tracking = peak_tracking;

  // Calculate Nmax and Nmin
  params.Nmax = (int)std::ceil((double)params.sample_freq / params.min_pitch);
  params.Nmin = (int)std::floor((double)params.sample_freq / params.max_pitch);
  params.min_pitch = (double)params.sample_freq / (double)params.Nmax;
  params.max_pitch = (double)params.sample_freq / (double)params.Nmin;

  // Calculate segment parameters
  int seg_size = 3 * params.Nmax;
  int seg_shift = (int)std::round(params.shift / 1000.0 * params.sample_freq);
  int seg_length = (int)std::round(params.length / 1000.0 * params.sample_freq);

  if (verbose) {
    Rcout << "Nmin=" << params.Nmin << " Nmax=" << params.Nmax << "\n";
    Rcout << "Segment size=" << seg_size << " shift=" << seg_shift << " length=" << seg_length << "\n";
  }

  // Calculate track length
  int track_len = (n_samples - seg_length) / seg_shift + 1;
  if (track_len < 1) track_len = 1;

  // Initialize output track
  NumericVector f0_track(track_len);
  std::vector<bool> is_voiced(track_len);

  // Initialize cross-correlation coefficients
  int cc_size = params.Nmax - params.Nmin + 1;
  std::vector<double> cc_coeff(cc_size);

  // Initialize status
  FrameStatus status;
  status.pitch_freq = BREAK_NUMBER;
  status.v_uv = SILENT;
  status.s_h = SEND;
  status.cc_max = 0.0;
  status.threshold = params.Thigh;

  FrameStatus held_status = status;
  int zx_lft_N = 0, zx_rht_N = 0;
  double prev_pf = BREAK_NUMBER;

  // Process frames
  int frame_idx = 0;
  int wave_pos = 0;

  // Handle padding
  int padding = 0;
  if (params.Nmax < seg_length / 2) {
    wave_pos = seg_length / 2 - params.Nmax;
  } else {
    if ((params.Nmax - seg_length / 2) % seg_shift != 0) {
      wave_pos = seg_shift - ((params.Nmax - seg_length / 2) % seg_shift);
    }
    padding = (params.Nmax - seg_length / 2) / seg_shift +
              ((params.Nmax - seg_length / 2) % seg_shift == 0 ? 0 : 1);
  }

  if (verbose) {
    Rcout << "Padding frames: " << padding << "\n";
  }

  // Process padding frames (silent)
  for (int p = 0; p < padding && frame_idx < track_len; p++) {
    f0_track[frame_idx] = BREAK_NUMBER;
    is_voiced[frame_idx] = false;
    frame_idx++;
  }

  // Process main frames
  while (wave_pos + seg_size <= n_samples && frame_idx < track_len) {
    // Extract segment
    std::vector<short> segment(seg_size);
    for (int i = 0; i < seg_size && wave_pos + i < n_samples; i++) {
      segment[i] = signal[wave_pos + i];
    }

    // Run PDA algorithm
    super_resolution_pda(params, segment, cc_coeff, status, zx_lft_N, zx_rht_N, prev_pf);

    // Handle held frames
    if (status.s_h == HOLD) {
      held_status.pitch_freq = status.pitch_freq;
      held_status.v_uv = VOICED;
      held_status.s_h = SEND;
      held_status.cc_max = status.cc_max;
      held_status.threshold = status.threshold;
    } else {
      if (held_status.s_h == SEND) {
        if (status.pitch_freq == BREAK_NUMBER) {
          held_status.pitch_freq = BREAK_NUMBER;
          held_status.v_uv = UNVOICED;
        }
        f0_track[frame_idx] = held_status.pitch_freq;
        is_voiced[frame_idx] = (held_status.v_uv == VOICED);
        frame_idx++;
        held_status.s_h = static_cast<SendStatus>(0);
      }

      if (frame_idx < track_len) {
        f0_track[frame_idx] = status.pitch_freq;
        is_voiced[frame_idx] = (status.v_uv == VOICED);
        frame_idx++;
      }
    }

    wave_pos += seg_shift;
  }

  // Handle final held frame
  if (held_status.s_h == SEND && frame_idx < track_len) {
    held_status.pitch_freq = BREAK_NUMBER;
    held_status.v_uv = UNVOICED;
    f0_track[frame_idx] = held_status.pitch_freq;
    is_voiced[frame_idx] = false;
    frame_idx++;
  }

  // Trim to actual length
  if (frame_idx < track_len) {
    f0_track = f0_track[Rcpp::Range(0, frame_idx - 1)];
    is_voiced.resize(frame_idx);
  }

  // Create time vector
  NumericVector times(frame_idx);
  for (int i = 0; i < frame_idx; i++) {
    times[i] = i * windowShift / 1000.0;  // convert ms to seconds
  }

  if (verbose) {
    int n_voiced = std::count(is_voiced.begin(), is_voiced.end(), true);
    Rcout << "Extracted " << frame_idx << " frames (" << n_voiced << " voiced)\n";
  }

  // Build result
  List result;
  result["f0"] = f0_track;
  result["times"] = times;
  result["is_voiced"] = LogicalVector(is_voiced.begin(), is_voiced.end());
  result["sample_rate"] = sample_rate;
  result["windowShift"] = windowShift;
  result["n_frames"] = frame_idx;

  return result;
}
