#ifndef SUPERASSP_ORT_LOADER_H
#define SUPERASSP_ORT_LOADER_H

// onnxruntime C API — loaded at runtime via dlopen, not linked at build time
#include "onnxruntime_c_api.h"

#include <string>

namespace superassp {
namespace ort {

/// Get the OrtApi function-pointer table. Loads the library on first call.
/// Returns NULL if the library cannot be loaded.
const OrtApi* get_api();

/// Check whether onnxruntime is available (loadable) on this system.
bool is_available();

/// Get the onnxruntime version string. Empty if not available.
std::string version();

/// Override the library search path. Must be called before first get_api().
/// @param dir Directory containing libonnxruntime.{dylib,so,dll}
void set_lib_dir(const std::string& dir);

/// Get the currently resolved library path (after loading). Empty if not loaded.
std::string lib_path();

} // namespace ort
} // namespace superassp

#endif // SUPERASSP_ORT_LOADER_H
