// OpenSMILE C++ Wrapper for R
// Provides direct C++ integration for OpenSMILE feature extraction
// [[Rcpp::plugins(cpp11)]]

#include <Rcpp.h>
#include <vector>
#include <string>
#include <memory>
#include <cstring>
#include <fstream>

// OpenSMILE C API
extern "C" {
#include <smileapi/SMILEapi.h>
}

using namespace Rcpp;

// Feature collector class to store results from callback
class FeatureCollector {
public:
  std::vector<float> features;
  std::vector<std::string> feature_names;
  bool collected;
  bool verbose;
  
  FeatureCollector(bool v = false) : collected(false), verbose(v) {}
  
  void collect(const float *data, long vectorSize) {
    if (verbose) {
      Rcpp::Rcout << "Callback triggered! vectorSize = " << vectorSize << "\n";
    }
    features.clear();
    features.reserve(vectorSize);
    for (long i = 0; i < vectorSize; i++) {
      features.push_back(data[i]);
    }
    collected = true;
  }
};

// C callback wrapper (OpenSMILE expects C linkage)
extern "C" {
  static bool feature_callback(const float *data, long vectorSize, void *param) {
    FeatureCollector *collector = static_cast<FeatureCollector*>(param);
    collector->collect(data, vectorSize);
    return true;  // Continue processing
  }
}

//' Generic OpenSMILE Feature Extraction (C++ Implementation)
//'
//' Extracts acoustic features using any OpenSMILE configuration file
//' via the OpenSMILE C++ library directly.
//'
//' @param audio_obj An AsspDataObj containing audio data
//' @param config_file Path to OpenSMILE configuration file
//' @param feature_set_name Name of feature set (for verbose output)
//' @param verbose Print processing information (default: FALSE)
//' @return Named list with acoustic features
//' @export
// [[Rcpp::export]]
List opensmile_extract_cpp(SEXP audio_obj,
                            std::string config_file,
                            std::string feature_set_name = "features",
                            bool verbose = false) {
  
  // Validate input
  if (!Rf_inherits(audio_obj, "AsspDataObj")) {
    stop("Input must be an AsspDataObj");
  }
  
  List audio_list(audio_obj);
  if (!audio_list.containsElementNamed("audio")) {
    stop("AsspDataObj must contain 'audio' track");
  }
  
  // Extract audio data
  IntegerMatrix audio_matrix = audio_list["audio"];
  int sample_rate = as<int>(audio_list.attr("sampleRate"));
  int n_samples = audio_matrix.nrow();
  int n_channels = audio_matrix.ncol();
  
  if (verbose) {
    Rcout << "OpenSMILE " << feature_set_name << " extraction\n";
    Rcout << "Audio: " << n_samples << " samples at " << sample_rate << " Hz\n";
    Rcout << "Channels: " << n_channels << "\n";
    Rcout << "Config: " << config_file << "\n";
  }
  
  // Convert INT16 audio to 16-bit PCM bytes (little-endian)
  // OpenSMILE ExternalAudioSource expects raw PCM data
  std::vector<int16_t> pcm_data(n_samples);
  for (int i = 0; i < n_samples; i++) {
    // Mix channels if stereo (take first channel for now)
    pcm_data[i] = static_cast<int16_t>(audio_matrix(i, 0));
  }
  
  // Initialize OpenSMILE
  smileobj_t *smile = smile_new();
  if (smile == NULL) {
    stop("Failed to create OpenSMILE instance");
  }
  
  // Initialize with config file
  smileres_t res = smile_initialize(smile, config_file.c_str(), 
                                     0, NULL,  // No command-line options
                                     1,        // Log level 1 = warnings only
                                     0,        // Debug off
                                     0,        // Console output off
                                     NULL);    // No log file
  
  if (res != SMILE_SUCCESS) {
    const char *error = smile_error_msg(smile);
    std::string error_str = error ? error : "Unknown error";
    smile_free(smile);
    stop("Failed to initialize OpenSMILE: " + error_str);
  }
  
  // Set up feature collector and callback
  FeatureCollector collector(verbose);
  res = smile_extsink_set_data_callback(smile, "functionals", 
                                         feature_callback, &collector);
  
  if (res != SMILE_SUCCESS) {
    const char *error = smile_error_msg(smile);
    std::string error_str = error ? error : "Component not found";
    smile_free(smile);
    stop("Failed to set data callback: " + error_str);
  }
  
  if (verbose) {
    Rcout << "Writing audio data...\n";
  }
  
  // Write ALL audio data before running
  // The external audio source will buffer it and feed it to the processing chain
  int data_size = n_samples * sizeof(int16_t);
  
  // Write data in one go
  res = smile_extaudiosource_write_data(smile, "externalAudio", 
                                         pcm_data.data(), data_size);
  
  if (verbose) {
    if (res == SMILE_SUCCESS) {
      Rcout << "Successfully wrote " << data_size << " bytes\n";
    } else if (res == SMILE_NOT_WRITTEN) {
      Rcout << "Data buffered (NOT_WRITTEN - this is normal)\n";
    } else {
      Rcout << "Warning: write returned code " << res << "\n";
    }
  }
  
  // Signal end of input
  res = smile_extaudiosource_set_external_eoi(smile, "externalAudio");
  if (res != SMILE_SUCCESS) {
    const char *error = smile_error_msg(smile);
    std::string error_str = error ? error : "Failed to set EOI";
    smile_free(smile);
    stop("Failed to set end of input: " + error_str);
  }
  
  if (verbose) {
    Rcout << "Running OpenSMILE...\n";
  }
  
  // Run OpenSMILE processing
  res = smile_run(smile);
  if (res != SMILE_SUCCESS) {
    const char *error = smile_error_msg(smile);
    std::string error_str = error ? error : "Processing failed";
    smile_free(smile);
    stop("OpenSMILE processing failed: " + error_str);
  }
  
  if (verbose) {
    Rcout << "OpenSMILE processing completed, requesting abort to trigger functionals...\n";
  }
  
  // Request abort to trigger final functional processing
  // This is needed for frameMode=full configs that compute functionals at EOI
  res = smile_abort(smile);
  if (res != SMILE_SUCCESS && verbose) {
    Rcout << "Warning: smile_abort returned " << res << "\n";
  }
  
  // Give it a moment to finish processing
  if (verbose) {
    Rcout << "Checking if callback was triggered...\n";
    Rcout << "Collector.collected = " << collector.collected << "\n";
  }
  
  // Check if features were collected
  if (!collector.collected) {
    smile_free(smile);
    stop("No features were extracted (callback not triggered)");
  }
  
  // Get feature names from sink
  long numElements = 0;
  res = smile_extsink_get_num_elements(smile, "functionals", &numElements);
  if (res != SMILE_SUCCESS || numElements <= 0) {
    smile_free(smile);
    stop("Failed to get number of features");
  }
  
  collector.feature_names.reserve(numElements);
  for (long i = 0; i < numElements; i++) {
    const char *name = NULL;
    res = smile_extsink_get_element_name(smile, "functionals", i, &name);
    if (res == SMILE_SUCCESS && name != NULL) {
      collector.feature_names.push_back(name);
    } else {
      // Fallback to generic name
      collector.feature_names.push_back("feature_" + std::to_string(i));
    }
  }
  
  if (verbose) {
    Rcout << "Extracted " << collector.features.size() << " features\n";
    Rcout << "Feature names: " << collector.feature_names.size() << "\n";
  }
  
  // Clean up OpenSMILE
  smile_free(smile);
  
  // Build result list with named features
  List result = List::create();
  
  if (collector.features.size() != collector.feature_names.size()) {
    warning("Feature count mismatch: " + 
            std::to_string(collector.features.size()) + " values, " +
            std::to_string(collector.feature_names.size()) + " names");
  }
  
  size_t n_features = std::min(collector.features.size(), 
                                collector.feature_names.size());
  
  for (size_t i = 0; i < n_features; i++) {
    result[collector.feature_names[i]] = collector.features[i];
  }
  
  if (verbose) {
    Rcout << feature_set_name << " extraction complete\n";
  }
  
  return result;
}

//' OpenSMILE GeMAPS Feature Extraction (C++ Implementation)
//'
//' Wrapper function for GeMAPS feature extraction
//'
//' @param audio_obj An AsspDataObj containing audio data
//' @param config_file Path to OpenSMILE configuration file
//' @param verbose Print processing information (default: FALSE)
//' @return Named list with 62 GeMAPS features
//' @export
// [[Rcpp::export]]
List opensmile_gemaps_cpp(SEXP audio_obj,
                          std::string config_file,
                          bool verbose = false) {
  return opensmile_extract_cpp(audio_obj, config_file, "GeMAPS", verbose);
}
