// ort_loader.cpp — Runtime loading of onnxruntime shared library via dlopen
//
// This file implements dynamic loading of the onnxruntime C API.
// The package compiles without onnxruntime present. The shared library
// is loaded at runtime when first needed. Users install it via
// install_onnxruntime() which downloads prebuilt binaries.

#include "ort_loader.h"
#include <Rcpp.h>

#include <string>
#include <mutex>

#ifdef _WIN32
#  define WIN32_LEAN_AND_MEAN
#  include <windows.h>
#  define ORT_LIB_NAME "onnxruntime.dll"
#else
#  include <dlfcn.h>
#  ifdef __APPLE__
#    define ORT_LIB_NAME "libonnxruntime.dylib"
#  else
#    define ORT_LIB_NAME "libonnxruntime.so"
#  endif
#endif

namespace superassp {
namespace ort {

// ---------- Module-level state ----------

static std::once_flag g_init_flag;
static const OrtApi* g_api = NULL;
static std::string g_lib_path;      // resolved path after successful load
static std::string g_override_dir;  // user-specified dir (before load)

// Platform-specific handle
#ifdef _WIN32
static HMODULE g_lib_handle = NULL;
#else
static void* g_lib_handle = NULL;
#endif

// ---------- Internal helpers ----------

// Resolve library directory using R-level options/env/defaults.
// Called from C++ but queries R for the cache path.
static std::string resolve_lib_dir() {
  // 1. User override via set_lib_dir()
  if (!g_override_dir.empty()) {
    return g_override_dir;
  }

  // 2. R option: superassp.onnxruntime.path
  try {
    Rcpp::Function getOption("getOption");
    SEXP val = getOption("superassp.onnxruntime.path");
    if (!Rf_isNull(val) && Rf_isString(val) && Rf_length(val) > 0) {
      std::string opt = Rcpp::as<std::string>(val);
      if (!opt.empty()) return opt;
    }
  } catch (...) {}

  // 3. Environment variable ONNXRUNTIME_LIB_PATH
  const char* env = std::getenv("ONNXRUNTIME_LIB_PATH");
  if (env && env[0] != '\0') {
    return std::string(env);
  }

  // 4. Default: tools::R_user_dir("superassp", "cache")/onnxruntime/lib
  try {
    Rcpp::Function r_user_dir("R_user_dir", Rcpp::Environment::namespace_env("tools"));
    std::string cache_dir = Rcpp::as<std::string>(r_user_dir("superassp", "cache"));
    return cache_dir + "/onnxruntime/lib";
  } catch (...) {}

  return "";
}

// Try to load the library from a specific directory. Returns true on success.
static bool try_load_from_dir(const std::string& dir) {
  if (dir.empty()) return false;

  std::string full_path = dir + "/" + ORT_LIB_NAME;

#ifdef _WIN32
  // Convert to wide string for LoadLibraryW
  int wlen = MultiByteToWideChar(CP_UTF8, 0, full_path.c_str(), -1, NULL, 0);
  if (wlen <= 0) return false;
  std::vector<wchar_t> wpath(wlen);
  MultiByteToWideChar(CP_UTF8, 0, full_path.c_str(), -1, &wpath[0], wlen);
  g_lib_handle = LoadLibraryW(&wpath[0]);
#else
  g_lib_handle = dlopen(full_path.c_str(), RTLD_NOW | RTLD_LOCAL);
#endif

  if (g_lib_handle) {
    g_lib_path = full_path;
    return true;
  }
  return false;
}

// Try system-level load (no path — relies on LD_LIBRARY_PATH, DYLD_LIBRARY_PATH, etc.)
static bool try_load_system() {
#ifdef _WIN32
  g_lib_handle = LoadLibraryA(ORT_LIB_NAME);
#else
  g_lib_handle = dlopen(ORT_LIB_NAME, RTLD_NOW | RTLD_LOCAL);
#endif
  if (g_lib_handle) {
    g_lib_path = ORT_LIB_NAME;  // system path, no full resolution
    return true;
  }
  return false;
}

// Resolve OrtGetApiBase from loaded library, get OrtApi*.
static bool resolve_api() {
  if (!g_lib_handle) return false;

  typedef const OrtApiBase* (*OrtGetApiBaseFn)(void);
  OrtGetApiBaseFn get_api_base = NULL;

#ifdef _WIN32
  get_api_base = (OrtGetApiBaseFn)GetProcAddress(g_lib_handle, "OrtGetApiBase");
#else
  get_api_base = (OrtGetApiBaseFn)dlsym(g_lib_handle, "OrtGetApiBase");
#endif

  if (!get_api_base) return false;

  const OrtApiBase* api_base = get_api_base();
  if (!api_base) return false;

  g_api = api_base->GetApi(ORT_API_VERSION);
  return (g_api != NULL);
}

// The one-shot initializer
static void do_init() {
  // Try resolved dir first, then system paths
  std::string dir = resolve_lib_dir();
  if (!try_load_from_dir(dir)) {
    try_load_system();
  }

  if (g_lib_handle) {
    if (!resolve_api()) {
      // Library loaded but API resolution failed — clean up
#ifdef _WIN32
      FreeLibrary(g_lib_handle);
#else
      dlclose(g_lib_handle);
#endif
      g_lib_handle = NULL;
      g_lib_path.clear();
      g_api = NULL;
    }
  }
}

// ---------- Public API ----------

const OrtApi* get_api() {
  std::call_once(g_init_flag, do_init);
  return g_api;
}

bool is_available() {
  return get_api() != NULL;
}

std::string version() {
  const OrtApi* api = get_api();
  if (!api) return "";

  // OrtApiBase has GetVersionString
  typedef const OrtApiBase* (*OrtGetApiBaseFn)(void);
  OrtGetApiBaseFn get_api_base = NULL;

#ifdef _WIN32
  get_api_base = (OrtGetApiBaseFn)GetProcAddress(g_lib_handle, "OrtGetApiBase");
#else
  get_api_base = (OrtGetApiBaseFn)dlsym(g_lib_handle, "OrtGetApiBase");
#endif

  if (get_api_base) {
    const OrtApiBase* base = get_api_base();
    if (base) {
      const char* ver = base->GetVersionString();
      if (ver) return std::string(ver);
    }
  }
  return "";
}

void set_lib_dir(const std::string& dir) {
  g_override_dir = dir;
}

std::string lib_path() {
  get_api();  // ensure initialized
  return g_lib_path;
}

void unload_library() {
  g_api = NULL;
  if (g_lib_handle) {
#ifdef _WIN32
    FreeLibrary(g_lib_handle);
#else
    dlclose(g_lib_handle);
#endif
    g_lib_handle = NULL;
  }
  g_lib_path.clear();
  // Reset once_flag so the library can be re-loaded if the package is re-loaded
  new (&g_init_flag) std::once_flag();
}

} // namespace ort
} // namespace superassp

// ---------- Rcpp exports ----------

// [[Rcpp::export]]
bool ort_available_cpp() {
  return superassp::ort::is_available();
}

// [[Rcpp::export]]
std::string ort_version_cpp() {
  return superassp::ort::version();
}

// [[Rcpp::export]]
std::string ort_lib_path_cpp() {
  return superassp::ort::lib_path();
}

// [[Rcpp::export]]
void ort_set_lib_dir_cpp(std::string dir) {
  superassp::ort::set_lib_dir(dir);
}
