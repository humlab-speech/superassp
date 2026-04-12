// ort_session.cpp — ONNX Runtime inference session wrapper
//
// Provides Rcpp-exported functions for creating sessions, running inference,
// and querying model I/O metadata. Sessions are stored as external pointers
// with custom release via Rcpp::XPtr.

#include "ort_session.h"
#include "crepe_inference.h"
#include <Rcpp.h>

#include <string>
#include <vector>
#include <stdexcept>
#include <mutex>

namespace superassp {
namespace ort {

// ---------- Shared ORT environment (singleton) ----------

static std::once_flag g_env_flag;
static OrtEnv* g_env = NULL;
static std::string g_env_error;  // deferred error from call_once

OrtEnv* OrtSessionWrapper::shared_env() {
  std::call_once(g_env_flag, [&]() {
    const OrtApi* api = get_api();
    if (!api) return;

    OrtStatus* status = api->CreateEnv(ORT_LOGGING_LEVEL_WARNING,
                                        "superassp", &g_env);
    if (status) {
      const char* msg = api->GetErrorMessage(status);
      g_env_error = std::string("ORT CreateEnv failed: ") + (msg ? msg : "unknown");
      api->ReleaseStatus(status);
      // Do NOT throw inside call_once — flag the error for caller
    }
  });
  // Report deferred error outside the call_once lambda
  if (!g_env && !g_env_error.empty()) {
    Rcpp::stop(g_env_error);
  }
  return g_env;
}

// ---------- Status checking ----------

void OrtSessionWrapper::check_status(OrtStatus* status) {
  if (status) {
    const char* msg = api_->GetErrorMessage(status);
    std::string err = msg ? std::string(msg) : "unknown ORT error";
    api_->ReleaseStatus(status);
    Rcpp::stop("onnxruntime error: " + err);
  }
}

// ---------- Constructor / Destructor ----------

OrtSessionWrapper::OrtSessionWrapper(const std::string& model_path, int num_threads)
  : api_(NULL), session_(NULL), session_options_(NULL), allocator_(NULL) {

  api_ = get_api();
  if (!api_) {
    Rcpp::stop("onnxruntime not available. Install with install_onnxruntime().");
  }

  OrtEnv* env = shared_env();
  if (!env) {
    Rcpp::stop("Failed to create ORT environment.");
  }

  // Session options
  check_status(api_->CreateSessionOptions(&session_options_));

  if (num_threads > 0) {
    check_status(api_->SetIntraOpNumThreads(session_options_, num_threads));
  }

  // Disable inter-op parallelism (single model, not needed)
  check_status(api_->SetInterOpNumThreads(session_options_, 1));

  // Use sequential execution (simpler, fine for our models)
  check_status(api_->SetSessionExecutionMode(session_options_,
                                              ORT_SEQUENTIAL));

  // Create session
#ifdef _WIN32
  // Windows needs wide string for model path
  int wlen = MultiByteToWideChar(CP_UTF8, 0, model_path.c_str(), -1, NULL, 0);
  std::vector<wchar_t> wpath(wlen);
  MultiByteToWideChar(CP_UTF8, 0, model_path.c_str(), -1, &wpath[0], wlen);
  check_status(api_->CreateSession(env, &wpath[0], session_options_, &session_));
#else
  check_status(api_->CreateSession(env, model_path.c_str(),
                                    session_options_, &session_));
#endif

  // Default allocator
  check_status(api_->GetAllocatorWithDefaultOptions(&allocator_));
}

OrtSessionWrapper::~OrtSessionWrapper() {
  // Release in reverse order of creation
  // Note: allocator from GetAllocatorWithDefaultOptions is not owned by us
  if (session_) {
    api_->ReleaseSession(session_);
    session_ = NULL;
  }
  if (session_options_) {
    api_->ReleaseSessionOptions(session_options_);
    session_options_ = NULL;
  }
}

// ---------- Inference ----------

Rcpp::List OrtSessionWrapper::run(
    const std::vector<std::string>& input_names,
    const std::vector<std::vector<float> >& input_data,
    const std::vector<std::vector<int64_t> >& input_shapes,
    const std::vector<std::string>& output_names) {

  size_t n_inputs = input_names.size();
  if (input_data.size() != n_inputs || input_shapes.size() != n_inputs) {
    Rcpp::stop("input_names, input_data, and input_shapes must have same length");
  }

  // --- Prepare input tensors ---
  OrtMemoryInfo* mem_info = NULL;
  check_status(api_->CreateCpuMemoryInfo(OrtArenaAllocator, OrtMemTypeDefault,
                                          &mem_info));

  std::vector<OrtValue*> input_tensors(n_inputs, NULL);
  std::vector<const char*> input_name_ptrs(n_inputs);

  for (size_t i = 0; i < n_inputs; i++) {
    input_name_ptrs[i] = input_names[i].c_str();

    // Compute expected element count from shape
    int64_t n_elem = 1;
    for (size_t d = 0; d < input_shapes[i].size(); d++) {
      n_elem *= input_shapes[i][d];
    }
    if ((int64_t)input_data[i].size() != n_elem) {
      api_->ReleaseMemoryInfo(mem_info);
      Rcpp::stop("Input '%s': data length %d does not match shape (expected %d)",
                  input_names[i].c_str(),
                  (int)input_data[i].size(), (int)n_elem);
    }

    check_status(api_->CreateTensorWithDataAsOrtValue(
      mem_info,
      (void*)input_data[i].data(),
      input_data[i].size() * sizeof(float),
      input_shapes[i].data(),
      input_shapes[i].size(),
      ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT,
      &input_tensors[i]
    ));
  }

  api_->ReleaseMemoryInfo(mem_info);

  // --- Determine output names ---
  size_t n_model_outputs = 0;
  check_status(api_->SessionGetOutputCount(session_, &n_model_outputs));

  std::vector<std::string> out_names_storage;
  std::vector<const char*> output_name_ptrs;
  size_t n_outputs;

  if (output_names.empty()) {
    // Fetch all output names from model
    n_outputs = n_model_outputs;
    out_names_storage.resize(n_outputs);
    output_name_ptrs.resize(n_outputs);
    for (size_t i = 0; i < n_outputs; i++) {
      char* name = NULL;
      check_status(api_->SessionGetOutputName(session_, i, allocator_, &name));
      out_names_storage[i] = name;
      output_name_ptrs[i] = out_names_storage[i].c_str();
      check_status(api_->AllocatorFree(allocator_, name));
    }
  } else {
    n_outputs = output_names.size();
    output_name_ptrs.resize(n_outputs);
    for (size_t i = 0; i < n_outputs; i++) {
      output_name_ptrs[i] = output_names[i].c_str();
    }
    out_names_storage = output_names;
  }

  // --- Run inference ---
  std::vector<OrtValue*> output_tensors(n_outputs, NULL);

  check_status(api_->Run(
    session_, NULL,
    input_name_ptrs.data(), (const OrtValue* const*)input_tensors.data(), n_inputs,
    output_name_ptrs.data(), n_outputs,
    output_tensors.data()
  ));

  // --- Extract outputs to R ---
  Rcpp::List result(n_outputs);
  Rcpp::CharacterVector result_names(n_outputs);

  for (size_t i = 0; i < n_outputs; i++) {
    result_names[i] = out_names_storage[i];

    // Get tensor type info
    OrtTensorTypeAndShapeInfo* type_info = NULL;
    check_status(api_->GetTensorTypeAndShape(output_tensors[i], &type_info));

    // Get shape
    size_t n_dims = 0;
    check_status(api_->GetDimensionsCount(type_info, &n_dims));
    std::vector<int64_t> shape(n_dims);
    check_status(api_->GetDimensions(type_info, shape.data(), n_dims));

    // Get element count
    size_t n_elem = 0;
    check_status(api_->GetTensorShapeElementCount(type_info, &n_elem));

    api_->ReleaseTensorTypeAndShapeInfo(type_info);

    // Get data pointer
    float* data = NULL;
    check_status(api_->GetTensorMutableData(output_tensors[i], (void**)&data));

    // Copy to R numeric vector
    Rcpp::NumericVector rv(n_elem);
    for (size_t j = 0; j < n_elem; j++) {
      rv[j] = (double)data[j];
    }

    // Attach shape as attribute
    Rcpp::IntegerVector shape_attr(n_dims);
    for (size_t d = 0; d < n_dims; d++) {
      shape_attr[d] = (int)shape[d];
    }
    rv.attr("shape") = shape_attr;

    result[i] = rv;
  }

  result.attr("names") = result_names;

  // --- Cleanup ---
  for (size_t i = 0; i < n_inputs; i++) {
    if (input_tensors[i]) api_->ReleaseValue(input_tensors[i]);
  }
  for (size_t i = 0; i < n_outputs; i++) {
    if (output_tensors[i]) api_->ReleaseValue(output_tensors[i]);
  }

  return result;
}

// ---------- Model metadata ----------

Rcpp::List OrtSessionWrapper::get_node_info(bool is_input) {
  size_t count = 0;
  if (is_input) {
    check_status(api_->SessionGetInputCount(session_, &count));
  } else {
    check_status(api_->SessionGetOutputCount(session_, &count));
  }

  Rcpp::List info(count);
  Rcpp::CharacterVector names(count);

  for (size_t i = 0; i < count; i++) {
    // Get name
    char* name = NULL;
    if (is_input) {
      check_status(api_->SessionGetInputName(session_, i, allocator_, &name));
    } else {
      check_status(api_->SessionGetOutputName(session_, i, allocator_, &name));
    }
    names[i] = std::string(name);
    check_status(api_->AllocatorFree(allocator_, name));

    // Get type info
    OrtTypeInfo* type_info = NULL;
    if (is_input) {
      check_status(api_->SessionGetInputTypeInfo(session_, i, &type_info));
    } else {
      check_status(api_->SessionGetOutputTypeInfo(session_, i, &type_info));
    }

    const OrtTensorTypeAndShapeInfo* tensor_info = NULL;
    check_status(api_->CastTypeInfoToTensorInfo(type_info, &tensor_info));

    // Element type
    ONNXTensorElementDataType elem_type;
    check_status(api_->GetTensorElementType(tensor_info, &elem_type));

    // Shape
    size_t n_dims = 0;
    check_status(api_->GetDimensionsCount(tensor_info, &n_dims));
    std::vector<int64_t> shape(n_dims);
    check_status(api_->GetDimensions(tensor_info, shape.data(), n_dims));

    Rcpp::IntegerVector shape_r(n_dims);
    for (size_t d = 0; d < n_dims; d++) {
      shape_r[d] = (int)shape[d];  // -1 for dynamic dims
    }

    // Type name mapping
    std::string type_str;
    switch (elem_type) {
      case ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT:   type_str = "float32"; break;
      case ONNX_TENSOR_ELEMENT_DATA_TYPE_DOUBLE:  type_str = "float64"; break;
      case ONNX_TENSOR_ELEMENT_DATA_TYPE_INT32:   type_str = "int32";   break;
      case ONNX_TENSOR_ELEMENT_DATA_TYPE_INT64:   type_str = "int64";   break;
      case ONNX_TENSOR_ELEMENT_DATA_TYPE_INT8:    type_str = "int8";    break;
      case ONNX_TENSOR_ELEMENT_DATA_TYPE_UINT8:   type_str = "uint8";   break;
      case ONNX_TENSOR_ELEMENT_DATA_TYPE_INT16:   type_str = "int16";   break;
      case ONNX_TENSOR_ELEMENT_DATA_TYPE_UINT16:  type_str = "uint16";  break;
      case ONNX_TENSOR_ELEMENT_DATA_TYPE_BOOL:    type_str = "bool";    break;
      case ONNX_TENSOR_ELEMENT_DATA_TYPE_STRING:  type_str = "string";  break;
      case ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT16: type_str = "float16"; break;
      default: type_str = "unknown"; break;
    }

    api_->ReleaseTypeInfo(type_info);

    info[i] = Rcpp::List::create(
      Rcpp::Named("name") = names[i],
      Rcpp::Named("shape") = shape_r,
      Rcpp::Named("type") = type_str
    );
  }

  info.attr("names") = names;
  return info;
}

Rcpp::List OrtSessionWrapper::input_info() {
  return get_node_info(true);
}

Rcpp::List OrtSessionWrapper::output_info() {
  return get_node_info(false);
}

// Release the shared OrtEnv and reset state so it can be re-created.
// Must be called BEFORE unload_library() since it uses g_api.
void cleanup_env() {
  if (g_env) {
    const OrtApi* api = get_api();
    if (api) {
      api->ReleaseEnv(g_env);
    }
    g_env = NULL;
  }
  g_env_error.clear();
  new (&g_env_flag) std::once_flag();
}

} // namespace ort
} // namespace superassp

// ---------- Rcpp exports ----------

// [[Rcpp::export]]
void ort_cleanup_cpp() {
  // 1. Release cached sessions (they hold OrtSession pointers)
  crepe_cleanup_session();
  // 2. Release the shared OrtEnv
  superassp::ort::cleanup_env();
  // 3. Unload the onnxruntime shared library
  superassp::ort::unload_library();
}

// [[Rcpp::export]]
SEXP ort_create_session_cpp(std::string model_path, int num_threads) {
  superassp::ort::OrtSessionWrapper* wrapper =
    new superassp::ort::OrtSessionWrapper(model_path, num_threads);

  Rcpp::XPtr<superassp::ort::OrtSessionWrapper> xptr(wrapper, true);
  xptr.attr("class") = "ort_session";
  return xptr;
}

// [[Rcpp::export]]
Rcpp::List ort_run_cpp(SEXP session_xptr,
                        std::vector<std::string> input_names,
                        Rcpp::List input_data_list,
                        Rcpp::List input_shapes_list,
                        Rcpp::Nullable<std::vector<std::string> > output_names) {

  Rcpp::XPtr<superassp::ort::OrtSessionWrapper> sess(session_xptr);

  // Convert R lists to C++ vectors
  size_t n = input_names.size();
  std::vector<std::vector<float> > input_data(n);
  std::vector<std::vector<int64_t> > input_shapes(n);

  for (size_t i = 0; i < n; i++) {
    Rcpp::NumericVector rv = Rcpp::as<Rcpp::NumericVector>(input_data_list[i]);
    input_data[i].resize(rv.size());
    for (int j = 0; j < rv.size(); j++) {
      input_data[i][j] = (float)rv[j];
    }

    Rcpp::IntegerVector sv = Rcpp::as<Rcpp::IntegerVector>(input_shapes_list[i]);
    input_shapes[i].resize(sv.size());
    for (int j = 0; j < sv.size(); j++) {
      input_shapes[i][j] = (int64_t)sv[j];
    }
  }

  std::vector<std::string> out_names;
  if (output_names.isNotNull()) {
    out_names = Rcpp::as<std::vector<std::string> >(output_names.get());
  }

  return sess->run(input_names, input_data, input_shapes, out_names);
}

// [[Rcpp::export]]
Rcpp::List ort_session_input_info_cpp(SEXP session_xptr) {
  Rcpp::XPtr<superassp::ort::OrtSessionWrapper> sess(session_xptr);
  return sess->input_info();
}

// [[Rcpp::export]]
Rcpp::List ort_session_output_info_cpp(SEXP session_xptr) {
  Rcpp::XPtr<superassp::ort::OrtSessionWrapper> sess(session_xptr);
  return sess->output_info();
}
