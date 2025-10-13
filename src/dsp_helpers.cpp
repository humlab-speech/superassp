// File: src/dsp_helpers.cpp
// Complete Rcpp implementation for dsp helper functions
// [[Rcpp::plugins(cpp11)]]

#include <Rcpp.h>
#include <string>
#include <vector>
#include <unordered_set>
#include <algorithm>

using namespace Rcpp;

//' Fast file extension extraction
 //' 
 //' @param paths Character vector of file paths
 //' @return Character vector of file extensions (without dot)
 //' @export
 // [[Rcpp::export]]
 CharacterVector fast_file_ext(CharacterVector paths) {
   int n = paths.size();
   CharacterVector result(n);
   
   for(int i = 0; i < n; i++) {
     std::string path = as<std::string>(paths[i]);
     size_t dot_pos = path.find_last_of('.');
     size_t slash_pos = path.find_last_of("/\\");
     
     // Only extract extension if dot comes after last slash
     if(dot_pos != std::string::npos && 
        (slash_pos == std::string::npos || dot_pos > slash_pos)) {
       result[i] = path.substr(dot_pos + 1);
     } else {
       result[i] = "";
     }
   }
   
   return result;
 }
 
 //' Fast check if extensions are in native formats
 //' 
 //' @param extensions Character vector of file extensions
 //' @param nativeFormats Character vector of native format extensions
 //' @return Logical vector indicating if each extension is native
 //' @export
 // [[Rcpp::export]]
 LogicalVector fast_is_native(CharacterVector extensions, 
                              CharacterVector nativeFormats) {
   int n = extensions.size();
   LogicalVector result(n);
   
   // Create hash set for O(1) lookup
   std::unordered_set<std::string> nativeSet;
   for(int i = 0; i < nativeFormats.size(); i++) {
     nativeSet.insert(as<std::string>(nativeFormats[i]));
   }
   
   // Check each extension
   for(int i = 0; i < n; i++) {
     std::string ext = as<std::string>(extensions[i]);
     result[i] = nativeSet.find(ext) != nativeSet.end();
   }
   
   return result;
 }
 
 //' Fast check if extensions are lossless
 //' 
 //' @param extensions Character vector of file extensions
 //' @param losslessFormats Character vector of known lossless format extensions
 //' @return Logical vector indicating if each extension is lossless
 //' @export
 // [[Rcpp::export]]
 LogicalVector fast_is_lossless(CharacterVector extensions,
                                CharacterVector losslessFormats) {
   int n = extensions.size();
   LogicalVector result(n);
   
   // Create hash set for O(1) lookup
   std::unordered_set<std::string> losslessSet;
   for(int i = 0; i < losslessFormats.size(); i++) {
     losslessSet.insert(as<std::string>(losslessFormats[i]));
   }
   
   // Check each extension
   for(int i = 0; i < n; i++) {
     std::string ext = as<std::string>(extensions[i]);
     result[i] = losslessSet.find(ext) != losslessSet.end();
   }
   
   return result;
 }
 
 //' Strip file:// protocol prefix efficiently
 //' 
 //' @param paths Character vector of paths
 //' @return Character vector with file:// removed
 //' @export
 // [[Rcpp::export]]
 CharacterVector fast_strip_file_protocol(CharacterVector paths) {
   int n = paths.size();
   CharacterVector result(n);
   
   const std::string prefix = "file://";
   const size_t prefix_len = prefix.length();
   
   for(int i = 0; i < n; i++) {
     std::string path = as<std::string>(paths[i]);
     
     if(path.compare(0, prefix_len, prefix) == 0) {
       result[i] = path.substr(prefix_len);
     } else {
       result[i] = path;
     }
   }
   
   return result;
 }
 
 //' Recycle time parameters efficiently
 //' 
 //' @param times Numeric vector of times (can be length 1 or n)
 //' @param n Target length
 //' @return Numeric vector of length n
 //' @export
 // [[Rcpp::export]]
 NumericVector fast_recycle_times(NumericVector times, int n) {
   int len = times.size();
   
   if(len == n) {
     return times;
   } else if(len == 1) {
     return NumericVector(n, times[0]);
   } else {
     stop("Time vector length must be 1 or match number of files");
   }
 }
 
 //' Validate time windows efficiently
 //' 
 //' @param beginTimes Numeric vector of begin times
 //' @param endTimes Numeric vector of end times
 //' @param filenames Character vector of filenames (for error messages)
 //' @return Logical TRUE if all valid, otherwise throws error
 //' @export
 // [[Rcpp::export]]
 bool fast_validate_times(NumericVector beginTimes,
                          NumericVector endTimes,
                          CharacterVector filenames) {
   int n = beginTimes.size();
   
   std::vector<int> invalid_indices;
   for(int i = 0; i < n; i++) {
     if(endTimes[i] != 0.0 && endTimes[i] <= beginTimes[i]) {
       invalid_indices.push_back(i);
     }
   }
   
   if(!invalid_indices.empty()) {
     std::string msg = "Invalid time windows for files: ";
     for(size_t i = 0; i < invalid_indices.size() && i < 5; i++) {
       int idx = invalid_indices[i];
       msg += as<std::string>(filenames[idx]) + " (" + 
         std::to_string(beginTimes[idx]) + "->" + 
         std::to_string(endTimes[idx]) + "), ";
     }
     if(invalid_indices.size() > 5) {
       msg += "and " + std::to_string(invalid_indices.size() - 5) + " more...";
     }
     stop(msg);
   }
   
   return true;
 }
 
 //' Batch rename list elements
 //' 
 //' @param objList List of objects with names to rename
 //' @param newNames Character vector of new names
 //' @return List with renamed elements
 //' @export
 // [[Rcpp::export]]
 List fast_rename_tracks(List objList, CharacterVector newNames) {
   int n = objList.size();
   List result = clone(objList);
   
   for(int i = 0; i < n; i++) {
     if(Rf_isNull(objList[i])) continue;
     
     SEXP obj = objList[i];
     if(TYPEOF(obj) == VECSXP) {
       List subObj = as<List>(obj);
       subObj.names() = newNames;
       result[i] = subObj;
     }
   }
   
   return result;
 }
 
 //' Generate output file paths efficiently
 //' 
 //' @param inputPaths Character vector of input file paths
 //' @param outputExt String for output extension
 //' @param isNative Logical vector indicating which files are already native
 //' @param convertTimewindow Logical vector indicating which need time window conversion
 //' @return Character vector of output paths
 //' @export
 // [[Rcpp::export]]
 CharacterVector fast_generate_output_paths(CharacterVector inputPaths,
                                            std::string outputExt,
                                            LogicalVector isNative,
                                            LogicalVector convertTimewindow) {
   int n = inputPaths.size();
   CharacterVector result(n);
   Function tempfile("tempfile");
   Function basename_r("basename");
   
   for(int i = 0; i < n; i++) {
     std::string path = as<std::string>(inputPaths[i]);
     
     if(convertTimewindow[i]) {
       // Generate temp file for time window conversion
       CharacterVector temp = tempfile(Named("pattern") = basename_r(inputPaths[i]),
                                       Named("fileext") = "." + outputExt);
       result[i] = as<std::string>(temp[0]);
     } else {
       // Replace extension with output extension
       size_t dot_pos = path.find_last_of('.');
       size_t slash_pos = path.find_last_of("/\\");
       
       if(dot_pos != std::string::npos && 
          (slash_pos == std::string::npos || dot_pos > slash_pos)) {
         result[i] = path.substr(0, dot_pos) + "." + outputExt;
       } else {
         result[i] = path + "." + outputExt;
       }
     }
   }
   
   return result;
 }
 
 //' Calculate adjusted times for file conversion
 //' 
 //' @param beginTimes Numeric vector of begin times
 //' @param endTimes Numeric vector of end times  
 //' @param durations Numeric vector of file durations
 //' @param windowShift Numeric window shift in milliseconds
 //' @param indices Integer vector of indices to process (0-based)
 //' @return DataFrame with start_time and total_time columns
 //' @export
 // [[Rcpp::export]]
 DataFrame fast_calculate_conversion_times(NumericVector beginTimes,
                                           NumericVector endTimes,
                                           NumericVector durations,
                                           double windowShift,
                                           IntegerVector indices) {
   int n = indices.size();
   NumericVector start_times(n);
   NumericVector total_times(n);
   
   double window_shift_sec = windowShift / 1000.0;
   
   for(int i = 0; i < n; i++) {
     int idx = indices[i];
     
     // Start time with buffer, but not negative
     start_times[i] = std::max(0.0, beginTimes[idx] - window_shift_sec);
     
     // End time (use duration if endTime is 0)
     double end = (endTimes[idx] == 0.0) ? durations[i] : endTimes[idx];
     
     // Total time needed
     total_times[i] = (end - beginTimes[idx]) + window_shift_sec;
   }
   
   return DataFrame::create(
     Named("start_time") = start_times,
     Named("total_time") = total_times,
     Named("stringsAsFactors") = false
   );
 }
 
 //' Extract basenames efficiently in batch
 //' 
 //' @param paths Character vector of file paths
 //' @return Character vector of base names without paths
 //' @export
 // [[Rcpp::export]]
 CharacterVector fast_basename(CharacterVector paths) {
   int n = paths.size();
   CharacterVector result(n);
   
   for(int i = 0; i < n; i++) {
     std::string path = as<std::string>(paths[i]);
     size_t last_slash = path.find_last_of("/\\");
     
     if(last_slash != std::string::npos) {
       result[i] = path.substr(last_slash + 1);
     } else {
       result[i] = path;
     }
   }
   
   return result;
 }
 
 //' Extract file names without extension in batch
 //' 
 //' @param paths Character vector of file paths
 //' @return Character vector of filenames without extension
 //' @export
 // [[Rcpp::export]]
 CharacterVector fast_file_path_sans_ext(CharacterVector paths) {
   int n = paths.size();
   CharacterVector result(n);
   
   for(int i = 0; i < n; i++) {
     std::string path = as<std::string>(paths[i]);
     size_t last_dot = path.find_last_of('.');
     size_t last_slash = path.find_last_of("/\\");
     
     // Only remove extension if dot comes after last slash
     if(last_dot != std::string::npos && 
        (last_slash == std::string::npos || last_dot > last_slash)) {
       result[i] = path.substr(0, last_dot);
     } else {
       result[i] = path;
     }
   }
   
   return result;
 }
 
 //' Build conversion dataframe efficiently
 //'
 //' @param audioPaths Character vector of input audio paths
 //' @param outputPaths Character vector of output paths
 //' @param startTimes Numeric vector of start times
 //' @param totalTimes Numeric vector of total/duration times
 //' @param indices Integer vector of indices in original file list (0-based)
 //' @return DataFrame ready for av::av_audio_convert
 //' @export
 // [[Rcpp::export]]
 DataFrame fast_build_conversion_df(CharacterVector audioPaths,
                                    CharacterVector outputPaths,
                                    NumericVector startTimes,
                                    NumericVector totalTimes,
                                    IntegerVector indices) {
   int n = indices.size();

   if(n == 0) {
     return DataFrame::create(
       Named("audio") = CharacterVector(),
       Named("output") = CharacterVector(),
       Named("start_time") = NumericVector(),
       Named("total_time") = NumericVector(),
       Named("stringsAsFactors") = false
     );
   }

   CharacterVector audio(n);
   CharacterVector output(n);

   for(int i = 0; i < n; i++) {
     int idx = indices[i];
     audio[i] = audioPaths[idx];
     output[i] = outputPaths[idx];
   }

   return DataFrame::create(
     Named("audio") = audio,
     Named("output") = output,
     Named("start_time") = startTimes,
     Named("total_time") = totalTimes,
     Named("stringsAsFactors") = false
   );
 }

 //' Fast check for native formats and lossless in one pass
 //'
 //' @param paths Character vector of file paths
 //' @param nativeFormats Character vector of native format extensions
 //' @param losslessFormats Character vector of known lossless format extensions
 //' @return List with extensions, is_native, and is_lossless vectors
 //' @export
 // [[Rcpp::export]]
 List fast_check_file_formats(CharacterVector paths,
                               CharacterVector nativeFormats,
                               CharacterVector losslessFormats) {
   int n = paths.size();
   CharacterVector extensions(n);
   LogicalVector is_native(n);
   LogicalVector is_lossless(n);

   // Create hash sets for O(1) lookup
   std::unordered_set<std::string> nativeSet;
   for(int i = 0; i < nativeFormats.size(); i++) {
     nativeSet.insert(as<std::string>(nativeFormats[i]));
   }

   std::unordered_set<std::string> losslessSet;
   for(int i = 0; i < losslessFormats.size(); i++) {
     losslessSet.insert(as<std::string>(losslessFormats[i]));
   }

   // Single pass through all files
   for(int i = 0; i < n; i++) {
     std::string path = as<std::string>(paths[i]);
     size_t dot_pos = path.find_last_of('.');
     size_t slash_pos = path.find_last_of("/\\");

     std::string ext = "";
     if(dot_pos != std::string::npos &&
        (slash_pos == std::string::npos || dot_pos > slash_pos)) {
       ext = path.substr(dot_pos + 1);
     }

     extensions[i] = ext;
     is_native[i] = nativeSet.find(ext) != nativeSet.end();
     is_lossless[i] = losslessSet.find(ext) != losslessSet.end();
   }

   return List::create(
     Named("extensions") = extensions,
     Named("is_native") = is_native,
     Named("is_lossless") = is_lossless
   );
 }