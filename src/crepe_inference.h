#ifndef SUPERASSP_CREPE_INFERENCE_H
#define SUPERASSP_CREPE_INFERENCE_H

#include <Rcpp.h>
#include <vector>
#include <string>

Rcpp::List crepe_inference_cpp(
    Rcpp::NumericVector audio,
    double sample_rate,
    std::string model_path,
    int hop_length,
    int batch_size,
    bool use_viterbi
);

// Release cached CREPE ONNX session. Called by ort_cleanup_cpp().
void crepe_cleanup_session();

#endif // SUPERASSP_CREPE_INFERENCE_H
