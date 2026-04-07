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

#endif // SUPERASSP_CREPE_INFERENCE_H
