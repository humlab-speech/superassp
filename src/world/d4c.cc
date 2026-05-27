//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
// Last update: 2017/02/01
//
// Band-aperiodicity estimation on the basis of the idea of D4C.
//
// Adapted for superassp/SPTK integration:
// Wrapped in namespace sptk::world:: so that unqualified calls to
// fft_execute(), InitializeForwardRealFFT(), DCCorrection(), NuttallWindow()
// etc. resolve to implementations already compiled in sptk::world::.
//
// Performance changes (superassp):
//  - Scratch-buffer struct eliminates all per-frame heap allocs (~15/frame).
//  - OpenMP parallelism over the voiced-frame loop in D4C().
//-----------------------------------------------------------------------------
#include "d4c.h"

#include <math.h>
#include <algorithm>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "world/common.h"
#include "world/constantnumbers.h"
#include "world/matlabfunctions.h"

#if 1
namespace sptk {
namespace world {
#endif

namespace {

// ---- Scratch buffers: allocated once per thread, reused across frames ----
struct D4CScratch {
  // Fixed-size arrays (fft_size_d4c / 2 + 1)
  std::vector<double> tmp_real;
  std::vector<double> tmp_imag;
  std::vector<double> centroid1;
  std::vector<double> centroid2;
  std::vector<double> static_centroid;
  std::vector<double> smoothed_power_spectrum;
  std::vector<double> static_group_delay;
  std::vector<double> smoothed_group_delay;
  std::vector<double> power_spectrum_ap;

  // Variable-size arrays (max size set at construction from floor_f0)
  std::vector<int>    base_index;
  std::vector<int>    safe_index;
  std::vector<double> window_buf;

  D4CScratch() {}  // for std::vector resize

  void init(int fft_size, int fs, double floor_f0) {
    int sp = fft_size / 2 + 1;
    tmp_real.assign(sp, 0.0);
    tmp_imag.assign(sp, 0.0);
    centroid1.assign(sp, 0.0);
    centroid2.assign(sp, 0.0);
    static_centroid.assign(sp, 0.0);
    smoothed_power_spectrum.assign(sp, 0.0);
    static_group_delay.assign(sp, 0.0);
    smoothed_group_delay.assign(sp, 0.0);
    power_spectrum_ap.assign(sp, 0.0);

    int max_half = matlab_round(4.0 * fs / floor_f0 / 2.0) + 1;
    int max_win  = max_half * 2 + 2;
    base_index.assign(max_win, 0);
    safe_index.assign(max_win, 0);
    window_buf.assign(max_win, 0.0);
  }
};

static void SetParametersForGetWindowedWaveform(int half_window_length,
    int x_length, double current_position, int fs, double current_f0,
    int window_type, double window_length_ratio, int *base_index,
    int *safe_index, double *window) {
  for (int i = -half_window_length; i <= half_window_length; ++i)
    base_index[i + half_window_length] = i;
  int origin = matlab_round(current_position * fs + 0.001);
  for (int i = 0; i <= half_window_length * 2; ++i)
    safe_index[i] =
      MyMinInt(x_length - 1, MyMaxInt(0, origin + base_index[i]));

  double position;
  if (window_type == world::kHanning) {
    for (int i = 0; i <= half_window_length * 2; ++i) {
      position = (2.0 * base_index[i] / window_length_ratio) / fs;
      window[i] = 0.5 * cos(world::kPi * position * current_f0) + 0.5;
    }
  } else {
    for (int i = 0; i <= half_window_length * 2; ++i) {
      position = (2.0 * base_index[i] / window_length_ratio) / fs;
      window[i] = 0.42 + 0.5 * cos(world::kPi * position * current_f0) +
        0.08 * cos(world::kPi * position * current_f0 * 2);
    }
  }
}

static void GetWindowedWaveform(const double *x, int x_length, int fs,
    double current_f0, double current_position, int window_type,
    double window_length_ratio, double *waveform, D4CScratch *sc) {
  int half_window_length =
    matlab_round(window_length_ratio * fs / current_f0 / 2.0);
  int win_size = half_window_length * 2 + 1;

  // Use pre-allocated scratch when available; fall back to local alloc
  // (fallback used only by D4CLoveTrainSub in the serial love-train pass)
  std::vector<int>    local_base, local_safe;
  std::vector<double> local_win;
  int    *base_idx, *safe_idx;
  double *win_ptr;
  if (sc) {
    base_idx = sc->base_index.data();
    safe_idx = sc->safe_index.data();
    win_ptr  = sc->window_buf.data();
  } else {
    local_base.resize(win_size);
    local_safe.resize(win_size);
    local_win.resize(win_size);
    base_idx = local_base.data();
    safe_idx = local_safe.data();
    win_ptr  = local_win.data();
  }

  SetParametersForGetWindowedWaveform(half_window_length, x_length,
      current_position, fs, current_f0, window_type, window_length_ratio,
      base_idx, safe_idx, win_ptr);

  for (int i = 0; i <= half_window_length * 2; ++i)
    waveform[i] =
      x[safe_idx[i]] * win_ptr[i] + randn() * world::kMySafeGuardMinimum;

  double tmp_weight1 = 0;
  double tmp_weight2 = 0;
  for (int i = 0; i <= half_window_length * 2; ++i) {
    tmp_weight1 += waveform[i];
    tmp_weight2 += win_ptr[i];
  }
  double weighting_coefficient = tmp_weight1 / tmp_weight2;
  for (int i = 0; i <= half_window_length * 2; ++i)
    waveform[i] -= win_ptr[i] * weighting_coefficient;
}

static void GetCentroid(const double *x, int x_length, int fs,
    double current_f0, int fft_size, double current_position,
    const ForwardRealFFT *forward_real_fft, double *centroid,
    D4CScratch *sc) {
  for (int i = 0; i < fft_size; ++i) forward_real_fft->waveform[i] = 0.0;
  GetWindowedWaveform(x, x_length, fs, current_f0,
      current_position, world::kBlackman, 4.0, forward_real_fft->waveform, sc);
  double power = 0.0;
  for (int i = 0; i <= matlab_round(2.0 * fs / current_f0) * 2; ++i)
    power += forward_real_fft->waveform[i] * forward_real_fft->waveform[i];
  for (int i = 0; i <= matlab_round(2.0 * fs / current_f0) * 2; ++i)
    forward_real_fft->waveform[i] /= sqrt(power);

  fft_execute(forward_real_fft->forward_fft);
  for (int i = 0; i <= fft_size / 2; ++i) {
    sc->tmp_real[i] = forward_real_fft->spectrum[i][0];
    sc->tmp_imag[i] = forward_real_fft->spectrum[i][1];
  }

  for (int i = 0; i < fft_size; ++i)
    forward_real_fft->waveform[i] *= i + 1.0;
  fft_execute(forward_real_fft->forward_fft);
  for (int i = 0; i <= fft_size / 2; ++i)
    centroid[i] = forward_real_fft->spectrum[i][0] * sc->tmp_real[i] +
      sc->tmp_imag[i] * forward_real_fft->spectrum[i][1];
}

static void GetStaticCentroid(const double *x, int x_length, int fs,
    double current_f0, int fft_size, double current_position,
    const ForwardRealFFT *forward_real_fft, double *static_centroid,
    D4CScratch *sc) {
  GetCentroid(x, x_length, fs, current_f0, fft_size,
      current_position - 0.25 / current_f0, forward_real_fft,
      sc->centroid1.data(), sc);
  GetCentroid(x, x_length, fs, current_f0, fft_size,
      current_position + 0.25 / current_f0, forward_real_fft,
      sc->centroid2.data(), sc);

  for (int i = 0; i <= fft_size / 2; ++i)
    static_centroid[i] = sc->centroid1[i] + sc->centroid2[i];

  DCCorrection(static_centroid, current_f0, fs, fft_size, static_centroid);
}

static void GetSmoothedPowerSpectrum(const double *x, int x_length, int fs,
    double current_f0, int fft_size, double current_position,
    const ForwardRealFFT *forward_real_fft, double *smoothed_power_spectrum,
    D4CScratch *sc) {
  for (int i = 0; i < fft_size; ++i) forward_real_fft->waveform[i] = 0.0;
  GetWindowedWaveform(x, x_length, fs, current_f0,
      current_position, world::kHanning, 4.0, forward_real_fft->waveform, sc);

  fft_execute(forward_real_fft->forward_fft);
  for (int i = 0; i <= fft_size / 2; ++i)
    smoothed_power_spectrum[i] =
      forward_real_fft->spectrum[i][0] * forward_real_fft->spectrum[i][0] +
      forward_real_fft->spectrum[i][1] * forward_real_fft->spectrum[i][1];

  DCCorrection(smoothed_power_spectrum, current_f0, fs, fft_size,
      smoothed_power_spectrum);
  LinearSmoothing(smoothed_power_spectrum, current_f0, fs, fft_size,
      smoothed_power_spectrum);
}

static void GetStaticGroupDelay(const double *static_centroid,
    const double *smoothed_power_spectrum, int fs, double f0,
    int fft_size, double *static_group_delay, D4CScratch *sc) {
  for (int i = 0; i <= fft_size / 2; ++i)
    static_group_delay[i] = static_centroid[i] / smoothed_power_spectrum[i];
  LinearSmoothing(static_group_delay, f0 / 2.0, fs, fft_size,
      static_group_delay);

  LinearSmoothing(static_group_delay, f0, fs, fft_size,
      sc->smoothed_group_delay.data());

  for (int i = 0; i <= fft_size / 2; ++i)
    static_group_delay[i] -= sc->smoothed_group_delay[i];
}

static void GetCoarseAperiodicity(const double *static_group_delay, int fs,
    int fft_size, int number_of_aperiodicities, const double *window,
    int window_length, const ForwardRealFFT *forward_real_fft,
    double *coarse_aperiodicity, D4CScratch *sc) {
  int boundary =
    matlab_round(fft_size * 8.0 / window_length);
  int half_window_length = window_length / 2;

  for (int i = 0; i < fft_size; ++i) forward_real_fft->waveform[i] = 0.0;

  int center;
  for (int i = 0; i < number_of_aperiodicities; ++i) {
    center =
      static_cast<int>(world::kFrequencyInterval * (i + 1) * fft_size / fs);
    for (int j = 0; j <= half_window_length * 2; ++j)
      forward_real_fft->waveform[j] =
        static_group_delay[center - half_window_length + j] * window[j];
    fft_execute(forward_real_fft->forward_fft);
    for (int j = 0 ; j <= fft_size / 2; ++j)
      sc->power_spectrum_ap[j] =
        forward_real_fft->spectrum[j][0] * forward_real_fft->spectrum[j][0] +
        forward_real_fft->spectrum[j][1] * forward_real_fft->spectrum[j][1];
    std::sort(sc->power_spectrum_ap.begin(),
              sc->power_spectrum_ap.begin() + fft_size / 2 + 1);
    for (int j = 1 ; j <= fft_size / 2; ++j)
      sc->power_spectrum_ap[j] += sc->power_spectrum_ap[j - 1];
    coarse_aperiodicity[i] =
      10 * log10(sc->power_spectrum_ap[fft_size / 2 - boundary - 1] /
                 sc->power_spectrum_ap[fft_size / 2]);
  }
}

static double D4CLoveTrainSub(const double *x, int fs, int x_length,
    double current_f0, double current_position, int f0_length, int fft_size,
    int boundary0, int boundary1, int boundary2,
    ForwardRealFFT *forward_real_fft) {
  double *power_spectrum = new double[fft_size];

  int window_length = matlab_round(1.5 * fs / current_f0) * 2 + 1;
  GetWindowedWaveform(x, x_length, fs, current_f0, current_position,
    world::kBlackman, 3.0, forward_real_fft->waveform, nullptr);

  for (int i = window_length; i < fft_size; ++i)
    forward_real_fft->waveform[i] = 0.0;
  fft_execute(forward_real_fft->forward_fft);

  for (int i = 0; i <= boundary0; ++i) power_spectrum[i] = 0.0;
  for (int i = boundary0 + 1; i < fft_size / 2 + 1; ++i)
    power_spectrum[i] =
    forward_real_fft->spectrum[i][0] * forward_real_fft->spectrum[i][0] +
    forward_real_fft->spectrum[i][1] * forward_real_fft->spectrum[i][1];
  for (int i = boundary0; i <= boundary2; ++i)
    power_spectrum[i] += +power_spectrum[i - 1];

  double aperiodicity0 = power_spectrum[boundary1] / power_spectrum[boundary2];
  delete[] power_spectrum;
  return aperiodicity0;
}

static void D4CLoveTrain(const double *x, int fs, int x_length,
    const double *f0, int f0_length, const double *temporal_positions,
    double *aperiodicity0) {
  double lowest_f0 = 40.0;
  int fft_size = static_cast<int>(pow(2.0, 1.0 +
    static_cast<int>(log(3.0 * fs / lowest_f0 + 1) / world::kLog2)));
  ForwardRealFFT forward_real_fft = { 0 };
  InitializeForwardRealFFT(fft_size, &forward_real_fft);

  int boundary0 = static_cast<int>(ceil(100.0 * fft_size / fs));
  int boundary1 = static_cast<int>(ceil(4000.0 * fft_size / fs));
  int boundary2 = static_cast<int>(ceil(7900.0 * fft_size / fs));
  for (int i = 0; i < f0_length; ++i) {
    if (f0[i] == 0.0) {
      aperiodicity0[i] = 0.0;
      continue;
    }
    aperiodicity0[i] = D4CLoveTrainSub(x, fs, x_length,
      MyMaxDouble(f0[i], lowest_f0), temporal_positions[i], f0_length,
      fft_size, boundary0, boundary1, boundary2, &forward_real_fft);
  }

  DestroyForwardRealFFT(&forward_real_fft);
}

static void D4CGeneralBody(const double *x, int x_length, int fs,
    double current_f0, int fft_size, double current_position,
    int number_of_aperiodicities, const double *window, int window_length,
    const ForwardRealFFT *forward_real_fft, double *coarse_aperiodicity,
    D4CScratch *sc) {
  GetStaticCentroid(x, x_length, fs, current_f0, fft_size, current_position,
      forward_real_fft, sc->static_centroid.data(), sc);
  GetSmoothedPowerSpectrum(x, x_length, fs, current_f0, fft_size,
      current_position, forward_real_fft,
      sc->smoothed_power_spectrum.data(), sc);
  GetStaticGroupDelay(sc->static_centroid.data(),
      sc->smoothed_power_spectrum.data(),
      fs, current_f0, fft_size, sc->static_group_delay.data(), sc);
  GetCoarseAperiodicity(sc->static_group_delay.data(), fs, fft_size,
      number_of_aperiodicities, window, window_length, forward_real_fft,
      coarse_aperiodicity, sc);

  for (int i = 0; i < number_of_aperiodicities; ++i)
    coarse_aperiodicity[i] = MyMinDouble(0.0,
      coarse_aperiodicity[i] + (current_f0 - 100) / 50.0);
}

static void InitializeAperiodicity(int f0_length, int fft_size,
    double **aperiodicity) {
  for (int i = 0; i < f0_length; ++i)
    for (int j = 0; j < fft_size / 2 + 1; ++j)
      aperiodicity[i][j] = 1.0 - world::kMySafeGuardMinimum;
}

static void GetAperiodicity(const double *coarse_frequency_axis,
    const double *coarse_aperiodicity, int number_of_aperiodicities,
    const double *frequency_axis, int fft_size, double *aperiodicity) {
  interp1(coarse_frequency_axis, coarse_aperiodicity,
    number_of_aperiodicities + 2, frequency_axis, fft_size / 2 + 1,
    aperiodicity);
  for (int i = 0; i <= fft_size / 2; ++i)
    aperiodicity[i] = pow(10.0, aperiodicity[i] / 20.0);
}

}  // namespace

void D4C(const double *x, int x_length, int fs,
    const double *temporal_positions, const double *f0, int f0_length,
    int fft_size, const D4COption *option, double **aperiodicity) {
  InitializeAperiodicity(f0_length, fft_size, aperiodicity);

  int fft_size_d4c = static_cast<int>(pow(2.0, 1.0 +
    static_cast<int>(log(4.0 * fs / world::kFloorF0D4C + 1) /
      world::kLog2)));

  int number_of_aperiodicities =
    static_cast<int>(MyMinDouble(world::kUpperLimit, fs / 2.0 -
      world::kFrequencyInterval) / world::kFrequencyInterval);
  int window_length =
    static_cast<int>(world::kFrequencyInterval * fft_size_d4c / fs) * 2 + 1;
  double *window = new double[window_length];
  NuttallWindow(window_length, window);

  double *aperiodicity0 = new double[f0_length];
  D4CLoveTrain(x, fs, x_length, f0, f0_length, temporal_positions,
      aperiodicity0);

  double *coarse_aperiodicity = new double[number_of_aperiodicities + 2];
  coarse_aperiodicity[0] = -60.0;
  coarse_aperiodicity[number_of_aperiodicities + 1] =
    -world::kMySafeGuardMinimum;
  double *coarse_frequency_axis = new double[number_of_aperiodicities + 2];
  for (int i = 0; i <= number_of_aperiodicities; ++i)
    coarse_frequency_axis[i] = i * world::kFrequencyInterval;
  coarse_frequency_axis[number_of_aperiodicities + 1] = fs / 2.0;

  double *frequency_axis = new double[fft_size / 2 + 1];
  for (int i = 0; i <= fft_size / 2; ++i)
    frequency_axis[i] = static_cast<double>(i) * fs / fft_size;

#ifdef _OPENMP
  int n_threads = omp_get_max_threads();
#else
  int n_threads = 1;
#endif

  // Per-thread scratch buffers and FFT structs
  std::vector<D4CScratch> scratches(n_threads);
  std::vector<ForwardRealFFT> ffts(n_threads);
  for (int t = 0; t < n_threads; ++t) {
    scratches[t].init(fft_size_d4c, fs, world::kFloorF0D4C);
    ffts[t] = {0};
    InitializeForwardRealFFT(fft_size_d4c, &ffts[t]);
  }

  // Per-thread coarse_aperiodicity scratch (size: number_of_aperiodicities+2)
  std::vector<std::vector<double>> coarse_ap_threads(n_threads,
    std::vector<double>(number_of_aperiodicities + 2, 0.0));

#ifdef _OPENMP
  #pragma omp parallel for
#endif
  for (int i = 0; i < f0_length; ++i) {
    if (f0[i] == 0 || aperiodicity0[i] <= option->threshold) continue;

#ifdef _OPENMP
    int tid = omp_get_thread_num();
#else
    int tid = 0;
#endif
    D4CScratch *sc = &scratches[tid];
    ForwardRealFFT *fft = &ffts[tid];
    double *coarse_ap_local = coarse_ap_threads[tid].data();

    coarse_ap_local[0] = coarse_aperiodicity[0];
    coarse_ap_local[number_of_aperiodicities + 1] =
      coarse_aperiodicity[number_of_aperiodicities + 1];

    D4CGeneralBody(x, x_length, fs, MyMaxDouble(world::kFloorF0D4C, f0[i]),
        fft_size_d4c, temporal_positions[i], number_of_aperiodicities, window,
        window_length, fft, &coarse_ap_local[1], sc);

    GetAperiodicity(coarse_frequency_axis, coarse_ap_local,
        number_of_aperiodicities, frequency_axis, fft_size, aperiodicity[i]);
  }

  for (int t = 0; t < n_threads; ++t)
    DestroyForwardRealFFT(&ffts[t]);

  delete[] aperiodicity0;
  delete[] coarse_frequency_axis;
  delete[] coarse_aperiodicity;
  delete[] window;
  delete[] frequency_axis;
}

void InitializeD4COption(D4COption *option) {
  option->threshold = world::kThreshold;
}

#if 1
}  // namespace world
}  // namespace sptk
#endif
