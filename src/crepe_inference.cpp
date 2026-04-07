// crepe_inference.cpp — C++ CREPE pitch inference via ONNX Runtime
//
// Implements: frame extraction, per-frame normalization, ORT inference,
// Viterbi decoding (360-state HMM), bin-to-Hz conversion, periodicity.
//
// Input: raw audio samples (float64 from R), sample_rate, model path, hop
// Output: list(f0, periodicity, times) — raw values before R post-processing

#include "crepe_inference.h"
#include "ort_session.h"

#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include <vector>
#include <string>

// ─── Constants ───────────────────────────────────────────────────────────────

static const int CREPE_FRAME_SIZE = 1024;
static const int CREPE_PITCH_BINS = 360;
static const double CREPE_CENTS_OFFSET = 1997.3794084376191;
static const double CREPE_CENTS_PER_BIN = 20.0;

// ─── Frame extraction ────────────────────────────────────────────────────────
// Zero-pad 512 each side, extract 1024-sample frames at hop_length intervals,
// per-frame normalize: subtract mean, divide by max(std, 1e-10).

static std::vector<float> extract_and_normalize_frames(
    const double* audio, int n_samples, int hop_length, int& n_frames_out
) {
    // Pad audio: 512 zeros on each side
    const int pad = CREPE_FRAME_SIZE / 2;  // 512
    const int padded_len = n_samples + 2 * pad;

    std::vector<float> padded(padded_len, 0.0f);
    for (int i = 0; i < n_samples; i++) {
        padded[pad + i] = static_cast<float>(audio[i]);
    }

    // Number of frames
    // torchcrepe: n_frames = 1 + (padded_len - frame_size) // hop_length
    // But actually, the original centers frames at each hop step starting at 0
    // in the original (unpadded) signal. With padding, center = pad + i*hop.
    // Frame start = center - pad = i * hop.
    n_frames_out = 1 + (n_samples - 1) / hop_length;
    if (n_frames_out < 1) n_frames_out = 1;

    // Extract and normalize frames
    std::vector<float> frames(n_frames_out * CREPE_FRAME_SIZE);

    for (int f = 0; f < n_frames_out; f++) {
        int start = f * hop_length;  // in padded signal
        float* frame = &frames[f * CREPE_FRAME_SIZE];

        // Copy frame
        for (int j = 0; j < CREPE_FRAME_SIZE; j++) {
            int idx = start + j;
            frame[j] = (idx < padded_len) ? padded[idx] : 0.0f;
        }

        // Compute mean
        double sum = 0.0;
        for (int j = 0; j < CREPE_FRAME_SIZE; j++) {
            sum += frame[j];
        }
        float mean = static_cast<float>(sum / CREPE_FRAME_SIZE);

        // Subtract mean and compute variance
        double var_sum = 0.0;
        for (int j = 0; j < CREPE_FRAME_SIZE; j++) {
            frame[j] -= mean;
            var_sum += (double)frame[j] * (double)frame[j];
        }
        float stddev = static_cast<float>(
            std::sqrt(var_sum / CREPE_FRAME_SIZE)
        );
        float scale = std::max(stddev, 1e-10f);

        // Divide by std
        for (int j = 0; j < CREPE_FRAME_SIZE; j++) {
            frame[j] /= scale;
        }
    }

    return frames;
}

// ─── Viterbi decoder ─────────────────────────────────────────────────────────
// 360-state HMM. Transition: T[i][j] = max(12 - |i-j|, 0), row-normalized.
// Uses log probabilities. Viterbi uses softmax (not sigmoid) probabilities.

static void softmax_inplace(float* probs, int n) {
    float max_val = *std::max_element(probs, probs + n);
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        probs[i] = std::exp(probs[i] - max_val);
        sum += probs[i];
    }
    float inv_sum = static_cast<float>(1.0 / sum);
    for (int i = 0; i < n; i++) {
        probs[i] *= inv_sum;
    }
}

static std::vector<int> viterbi_decode(
    const std::vector<float>& sigmoid_probs, int n_frames
) {
    const int N = CREPE_PITCH_BINS;  // 360
    const int BAND = 12;  // transition bandwidth

    // Precompute log transition weights (only within band)
    // T[i][j] = max(12 - |i-j|, 0) for |i-j| <= 11, normalized per row
    // Each row sums to 12+11+10+...+1+0+1+...+11 = 12 + 2*(11+10+...+1) = 12+110 = 132 for interior
    // But boundary rows have fewer nonzero entries.
    // We'll compute log(T[i][j]/row_sum) on the fly during transitions.

    // Actually, the band is 12, so T[i][j] = max(12-|i-j|, 0)
    // Nonzero when |i-j| < 12, i.e. |i-j| <= 11
    // For interior row: sum = 12 + 2*(11+10+...+1) = 12 + 110 = 122
    // Hmm, let me recalculate: values are 12,11,10,...,1 for d=0,1,...,11
    // Sum for interior = 12 + 2*(11+10+9+8+7+6+5+4+3+2+1) = 12 + 2*66 = 12+132 = 144
    // Wait: d=0 -> 12, d=1 -> 11, ..., d=11 -> 1, d=12 -> 0
    // So nonzero for d in [0,11], values: 12,11,10,9,8,7,6,5,4,3,2,1
    // Interior sum = 12+11+10+9+8+7+6+5+4+3+2+1 = 78, but both sides: 78 + 66 = 144? No.
    // d from -11 to +11: 1,2,...,11,12,11,...,2,1 -> sum = 2*(1+2+...+11)+12 = 2*66+12 = 144

    // Precompute log row sums for each state
    std::vector<float> log_row_sum(N);
    for (int i = 0; i < N; i++) {
        double s = 0.0;
        for (int j = std::max(0, i - BAND + 1); j <= std::min(N - 1, i + BAND - 1); j++) {
            s += (BAND - std::abs(i - j));
        }
        log_row_sum[i] = static_cast<float>(std::log(s));
    }

    // Convert sigmoid probs to log-softmax for Viterbi observation probabilities
    // torchcrepe uses softmax (not sigmoid) for Viterbi
    // We need to convert: sigmoid output → treat as logits → softmax → log
    // Since sigmoid(x) is monotonic with x, and we have sigmoid outputs,
    // we can recover logits: x = log(p/(1-p)), then apply log_softmax.
    // But simpler: just use log(softmax(logits)) where logits = log(p/(1-p))

    // Allocate Viterbi tables
    std::vector<float> V_prev(N), V_curr(N);
    std::vector<std::vector<int> > backptr(n_frames, std::vector<int>(N));

    // Working buffer for per-frame softmax
    std::vector<float> logits(N);

    // Initialize: first frame
    {
        const float* sp = &sigmoid_probs[0];
        for (int j = 0; j < N; j++) {
            // Convert sigmoid → logit for softmax
            float p = std::max(std::min(sp[j], 1.0f - 1e-7f), 1e-7f);
            logits[j] = std::log(p / (1.0f - p));
        }
        softmax_inplace(logits.data(), N);
        for (int j = 0; j < N; j++) {
            V_prev[j] = std::log(std::max(logits[j], 1e-30f));
        }
    }

    // Forward pass
    for (int t = 1; t < n_frames; t++) {
        // Compute observation log-probabilities for frame t
        const float* sp = &sigmoid_probs[t * N];
        for (int j = 0; j < N; j++) {
            float p = std::max(std::min(sp[j], 1.0f - 1e-7f), 1e-7f);
            logits[j] = std::log(p / (1.0f - p));
        }
        softmax_inplace(logits.data(), N);

        for (int j = 0; j < N; j++) {
            float log_obs = std::log(std::max(logits[j], 1e-30f));

            // Find best predecessor within band
            float best_val = -std::numeric_limits<float>::infinity();
            int best_prev = j;

            int lo = std::max(0, j - BAND + 1);
            int hi = std::min(N - 1, j + BAND - 1);

            for (int i = lo; i <= hi; i++) {
                float trans_weight = static_cast<float>(BAND - std::abs(i - j));
                float log_trans = std::log(trans_weight) - log_row_sum[i];
                float val = V_prev[i] + log_trans;
                if (val > best_val) {
                    best_val = val;
                    best_prev = i;
                }
            }

            V_curr[j] = best_val + log_obs;
            backptr[t][j] = best_prev;
        }

        std::swap(V_prev, V_curr);
    }

    // Backtrace
    std::vector<int> path(n_frames);
    path[n_frames - 1] = static_cast<int>(
        std::max_element(V_prev.begin(), V_prev.end()) - V_prev.begin()
    );
    for (int t = n_frames - 2; t >= 0; t--) {
        path[t] = backptr[t + 1][path[t + 1]];
    }

    return path;
}

// ─── Argmax decoder (simple alternative) ─────────────────────────────────────

static std::vector<int> argmax_decode(
    const std::vector<float>& probs, int n_frames
) {
    const int N = CREPE_PITCH_BINS;
    std::vector<int> bins(n_frames);
    for (int t = 0; t < n_frames; t++) {
        const float* p = &probs[t * N];
        bins[t] = static_cast<int>(
            std::max_element(p, p + N) - p
        );
    }
    return bins;
}

// ─── Bin to Hz conversion ────────────────────────────────────────────────────

static double bin_to_hz(int bin) {
    double cents = CREPE_CENTS_PER_BIN * bin + CREPE_CENTS_OFFSET;
    return 10.0 * std::pow(2.0, cents / 1200.0);
}

// ─── Main entry point ────────────────────────────────────────────────────────

// [[Rcpp::export]]
Rcpp::List crepe_inference_cpp(
    Rcpp::NumericVector audio,
    double sample_rate,
    std::string model_path,
    int hop_length,
    int batch_size,
    bool use_viterbi
) {
    int n_samples = audio.size();
    if (n_samples < CREPE_FRAME_SIZE) {
        Rcpp::stop("Audio too short for CREPE: need at least %d samples, got %d",
                    CREPE_FRAME_SIZE, n_samples);
    }

    // Step 1: Extract and normalize frames
    int n_frames = 0;
    std::vector<float> frames = extract_and_normalize_frames(
        audio.begin(), n_samples, hop_length, n_frames
    );

    // Step 2: Run ONNX inference in batches
    superassp::ort::OrtSessionWrapper session(model_path, 0);
    std::vector<float> all_probs(n_frames * CREPE_PITCH_BINS);

    for (int batch_start = 0; batch_start < n_frames; batch_start += batch_size) {
        int batch_end = std::min(batch_start + batch_size, n_frames);
        int this_batch = batch_end - batch_start;

        // Prepare input
        std::vector<float> batch_data(
            frames.begin() + batch_start * CREPE_FRAME_SIZE,
            frames.begin() + batch_end * CREPE_FRAME_SIZE
        );

        std::vector<std::string> input_names = {"frames"};
        std::vector<std::vector<float> > input_data = {batch_data};
        std::vector<std::vector<int64_t> > input_shapes = {
            {static_cast<int64_t>(this_batch),
             static_cast<int64_t>(CREPE_FRAME_SIZE)}
        };
        std::vector<std::string> output_names = {"probabilities"};

        Rcpp::List result = session.run(
            input_names, input_data, input_shapes, output_names
        );

        // Extract probabilities
        Rcpp::NumericVector probs_rv = result["probabilities"];
        float* dst = &all_probs[batch_start * CREPE_PITCH_BINS];
        for (int i = 0; i < this_batch * CREPE_PITCH_BINS; i++) {
            dst[i] = static_cast<float>(probs_rv[i]);
        }
    }

    // Step 3: Decode bins
    std::vector<int> bins;
    if (use_viterbi) {
        bins = viterbi_decode(all_probs, n_frames);
    } else {
        bins = argmax_decode(all_probs, n_frames);
    }

    // Step 4: Convert to Hz and extract periodicity
    Rcpp::NumericVector f0(n_frames);
    Rcpp::NumericVector periodicity(n_frames);
    Rcpp::NumericVector times(n_frames);

    double hop_seconds = static_cast<double>(hop_length) / sample_rate;

    for (int t = 0; t < n_frames; t++) {
        f0[t] = bin_to_hz(bins[t]);
        periodicity[t] = static_cast<double>(
            all_probs[t * CREPE_PITCH_BINS + bins[t]]
        );
        times[t] = t * hop_seconds;
    }

    return Rcpp::List::create(
        Rcpp::Named("f0") = f0,
        Rcpp::Named("periodicity") = periodicity,
        Rcpp::Named("times") = times,
        Rcpp::Named("n_frames") = n_frames,
        Rcpp::Named("hop_length") = hop_length,
        Rcpp::Named("sample_rate") = sample_rate
    );
}
