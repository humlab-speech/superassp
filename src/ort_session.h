#ifndef SUPERASSP_ORT_SESSION_H
#define SUPERASSP_ORT_SESSION_H

#include "ort_loader.h"
#include <Rcpp.h>

#include <string>
#include <vector>

namespace superassp {
namespace ort {

/// RAII wrapper around an ORT inference session.
/// Owns the session, env (shared), and default allocator.
class OrtSessionWrapper {
public:
  OrtSessionWrapper(const std::string& model_path, int num_threads);
  ~OrtSessionWrapper();

  // Non-copyable
  OrtSessionWrapper(const OrtSessionWrapper&) = delete;
  OrtSessionWrapper& operator=(const OrtSessionWrapper&) = delete;

  /// Run inference.
  /// @param input_names  Names of input tensors
  /// @param input_data   Flat float data for each input (row-major)
  /// @param input_shapes Shape vector for each input (e.g. {1, 1024})
  /// @param output_names Names of output tensors to fetch (empty = all)
  /// @return List of output tensors, each as NumericVector with "shape" attribute
  Rcpp::List run(
    const std::vector<std::string>& input_names,
    const std::vector<std::vector<float> >& input_data,
    const std::vector<std::vector<int64_t> >& input_shapes,
    const std::vector<std::string>& output_names
  );

  /// Get input tensor metadata: list of (name, shape, type)
  Rcpp::List input_info();

  /// Get output tensor metadata: list of (name, shape, type)
  Rcpp::List output_info();

private:
  const OrtApi* api_;
  OrtSession* session_;
  OrtSessionOptions* session_options_;
  OrtAllocator* allocator_;

  // Shared ORT environment (process-wide singleton)
  static OrtEnv* shared_env();

  // Helper: check ORT status, throw on error
  void check_status(OrtStatus* status);

  // Helper: get node info (inputs or outputs)
  Rcpp::List get_node_info(bool is_input);
};

/// Release the shared ORT environment. Called during package unload.
void cleanup_env();

} // namespace ort
} // namespace superassp

#endif // SUPERASSP_ORT_SESSION_H
