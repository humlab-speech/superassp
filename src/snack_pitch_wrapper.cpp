// Snack Pitch Tracker — Full 4-track wrapper around dp_f0
//
// Calls sptk::snack::init_dp_f0 / dp_f0 / free_dp_f0 directly
// (already compiled from SPTK/third_party/Snack/jkGetF0.cc)
// to recover all 4 output tracks: f0, voicing, rms, acpeak.
//
// The existing cGet_f0() only keeps f0, discarding the other three.

#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include "SPTK/generation/normal_distributed_random_value_generation.h"
#include "Snack/jkGetF0.h"   // F0_params, ckalloc/ckfree, etc.

// These functions are defined in jkGetF0.cc (sptk::snack namespace)
// but not declared in the header.  Declare them here for linkage.
namespace sptk { namespace snack {
  int init_dp_f0(double freq, F0_params *par, long *buffsize, long *sdstep);
  int dp_f0(float *fdata, int buff_size, int sdstep, double freq,
            F0_params *par,
            float **f0p_pt, float **vuvp_pt,
            float **rms_speech_pt, float **acpkp_pt,
            int *vecsize, int last_time);
  void free_dp_f0();
}}

using namespace Rcpp;

//' Snack Pitch Extraction — full 4-track output (C++ implementation)
//'
//' Normalized cross-correlation + dynamic-programming pitch tracker from the
//' Snack Sound Toolkit (Talkin, 1995).  Returns F0, voicing probability,
//' RMS energy and autocorrelation-peak per frame.
//'
//' @param audio_obj  An AsspDataObj containing audio data
//' @param minF       Minimum F0 in Hz (default 50)
//' @param maxF       Maximum F0 in Hz (default 550)
//' @param windowShift Frame shift in milliseconds (default 10)
//' @param voiceBias  Bias toward voiced hypothesis (default 0.0)
//' @param verbose    Print processing info (default FALSE)
//' @return List with f0, voicing, rms, acpeak (matrices), times, sample_rate, n_frames
//' @keywords internal
// [[Rcpp::export]]
List snackp_cpp(SEXP audio_obj,
                double minF = 50.0,
                double maxF = 550.0,
                double windowShift = 10.0,
                double voiceBias = 0.0,
                bool verbose = false) {

  if (!Rf_inherits(audio_obj, "AsspDataObj"))
    stop("Input must be an AsspDataObj");

  List audio_list(audio_obj);
  if (!audio_list.containsElementNamed("audio"))
    stop("AsspDataObj must contain 'audio' track");

  NumericMatrix audio_matrix = audio_list["audio"];
  int sample_rate = as<int>(audio_list.attr("sampleRate"));
  int n_samples   = audio_matrix.nrow();

  if (verbose) {
    Rcout << "snackp_cpp: " << n_samples << " samples @ "
          << sample_rate << " Hz, F0 " << minF << "-" << maxF << " Hz\n";
  }

  // --- Convert to std::vector<double> (first channel) ---
  std::vector<double> waveform(n_samples);
  for (int i = 0; i < n_samples; i++)
    waveform[i] = audio_matrix(i, 0);

  int frame_shift = static_cast<int>(windowShift * sample_rate / 1000.0);

  // --- Pad input with low-level noise (same logic as cGet_f0) ---
  using sptk::snack::ckalloc;
  using sptk::snack::ckfree;

  long sound_length = static_cast<long>(n_samples);
  double fsp = sample_rate * 10.0 / frame_shift;
  int alpha = static_cast<int>(0.00275 * fsp + 0.5);
  int beta  = static_cast<int>((9600.0 / minF - 168.0) * fsp / 96000.0 + 0.5);
  if (beta < 0) beta = 0;
  long pad_length   = (alpha + beta + 3) * frame_shift;
  long total_length = sound_length + pad_length;

  float *buf = (float *)ckalloc(sizeof(float) * total_length);
  double noise_sdev = 50.0;
  sptk::NormalDistributedRandomValueGeneration generator(1);
  double noise;
  for (long i = 0; i < sound_length; i++) {
    if (!generator.Get(&noise)) { ckfree(buf); stop("RNG failed"); }
    buf[i] = static_cast<float>(waveform[i] + noise * noise_sdev);
  }
  for (long i = sound_length; i < total_length; i++) {
    if (!generator.Get(&noise)) { ckfree(buf); stop("RNG failed"); }
    buf[i] = static_cast<float>(noise * noise_sdev);
  }

  // --- F0_params (same defaults as cGet_f0) ---
  sptk::snack::F0_params *par =
      (sptk::snack::F0_params *)ckalloc(sizeof(sptk::snack::F0_params));
  par->cand_thresh    = 0.3f;
  par->lag_weight     = 0.3f;
  par->freq_weight    = 0.02f;
  par->trans_cost     = 0.005f;
  par->trans_amp      = 0.5f;
  par->trans_spec     = 0.5f;
  par->voice_bias     = static_cast<float>(voiceBias);
  par->double_cost    = 0.35f;
  par->min_f0         = static_cast<float>(minF);
  par->max_f0         = static_cast<float>(maxF);
  par->frame_step     = static_cast<float>((double)frame_shift / sample_rate);
  par->wind_dur       = 0.0075f;
  par->n_cands        = 20;
  par->mean_f0        = 200;
  par->mean_f0_weight = 0.0f;
  par->conditioning   = 0;

  double sf = static_cast<double>(sample_rate);
  long total_samps = sound_length;

  if (total_samps < ((par->frame_step * 2.0) + par->wind_dur) * sf) {
    ckfree(buf); ckfree(par);
    stop("Input too short for Snack pitch analysis");
  }

  long buff_size, sdstep;
  if (sptk::snack::init_dp_f0(sf, par, &buff_size, &sdstep)
      || buff_size > INT_MAX || sdstep > INT_MAX) {
    ckfree(buf); ckfree(par);
    stop("init_dp_f0 failed");
  }

  if (buff_size > total_samps)
    buff_size = total_samps;

  long actsize = std::min(buff_size, sound_length);
  float *fdata = (float *)ckalloc(sizeof(float) * std::max(buff_size, sdstep));

  // Pre-allocate output arrays (generous upper bound)
  long max_frames = 5 + sound_length / frame_shift;
  std::vector<float> out_f0, out_vuv, out_rms, out_acpk;
  out_f0.reserve(max_frames);
  out_vuv.reserve(max_frames);
  out_rms.reserve(max_frames);
  out_acpk.reserve(max_frames);

  // --- Main dp_f0 loop (mirrors cGet_f0 but keeps all 4 tracks) ---
  long ndone = 0;
  float *f0p, *vuvp, *rms_speech, *acpkp;
  int vecsize;

  while (true) {
    int done = (actsize < buff_size) || (total_samps == buff_size);

    for (long i = 0; i < actsize; i++)
      fdata[i] = buf[i + ndone];

    if (sptk::snack::dp_f0(fdata, (int)actsize, (int)sdstep, sf, par,
                            &f0p, &vuvp, &rms_speech, &acpkp,
                            &vecsize, done)) {
      break;
    }

    // dp_f0 returns frames in reverse order
    for (int i = vecsize - 1; i >= 0; i--) {
      out_f0.push_back(f0p[i]);
      out_vuv.push_back(vuvp[i]);
      out_rms.push_back(rms_speech[i]);
      out_acpk.push_back(acpkp[i]);
    }

    if (done) break;

    ndone += sdstep;
    actsize = std::min(buff_size, sound_length - ndone);
    total_samps -= sdstep;
    if (actsize > total_samps)
      actsize = total_samps;
  }

  ckfree(fdata);
  ckfree(par);
  ckfree(buf);
  sptk::snack::free_dp_f0();

  int n_frames = static_cast<int>(out_f0.size());
  if (verbose)
    Rcout << "snackp_cpp: extracted " << n_frames << " frames\n";

  // --- Pack into R matrices ---
  NumericMatrix f0_mat(n_frames, 1);
  NumericMatrix vuv_mat(n_frames, 1);
  NumericMatrix rms_mat(n_frames, 1);
  NumericMatrix acpk_mat(n_frames, 1);
  NumericVector times(n_frames);

  for (int i = 0; i < n_frames; i++) {
    f0_mat(i, 0)   = out_f0[i];
    vuv_mat(i, 0)  = out_vuv[i];
    rms_mat(i, 0)  = out_rms[i];
    acpk_mat(i, 0) = out_acpk[i];
    times[i]        = i * windowShift / 1000.0;
  }

  return List::create(
    Named("f0")          = f0_mat,
    Named("voicing")     = vuv_mat,
    Named("rms")         = rms_mat,
    Named("acpeak")      = acpk_mat,
    Named("times")       = times,
    Named("sample_rate") = sample_rate,
    Named("n_frames")    = n_frames
  );
}
