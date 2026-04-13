#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/*
The following name(s) appear with different usages
e.g., with different numbers of arguments:

performAssp

This needs to be resolved in the tables and any declarations.
*/

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP AsspLpTypes_(void);
extern SEXP AsspSpectTypes_(void);
extern SEXP AsspWindowTypes_(void);
extern SEXP writeDObj_(SEXP, SEXP);

/* Rcpp exports */
extern SEXP _superassp_fast_file_ext(SEXP);
extern SEXP _superassp_fast_is_native(SEXP, SEXP);
extern SEXP _superassp_fast_is_lossless(SEXP, SEXP);
extern SEXP _superassp_fast_strip_file_protocol(SEXP);
extern SEXP _superassp_fast_recycle_times(SEXP, SEXP);
extern SEXP _superassp_fast_validate_times(SEXP, SEXP, SEXP);
extern SEXP _superassp_fast_rename_tracks(SEXP, SEXP);
extern SEXP _superassp_fast_generate_output_paths(SEXP, SEXP, SEXP, SEXP);
extern SEXP _superassp_fast_calculate_conversion_times(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _superassp_fast_basename(SEXP);
extern SEXP _superassp_fast_file_path_sans_ext(SEXP);
extern SEXP _superassp_fast_build_conversion_df(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _superassp_fast_check_file_formats(SEXP, SEXP, SEXP);
extern SEXP _superassp_estk_pda_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _superassp_estk_pitchmark_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _superassp_opensmile_gemaps_cpp(SEXP, SEXP, SEXP);
extern SEXP _superassp_opensmile_extract_cpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP _superassp_rapt_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _superassp_swipe_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _superassp_reaper_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _superassp_dio_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _superassp_harvest_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _superassp_d4c_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _superassp_sptk_mfcc_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _superassp_yin_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _superassp_pyin_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _superassp_momel_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP _superassp_tandem_pitch_cpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _superassp_iaif_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _superassp_extract_vq_params_cpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _superassp_gfmiaif_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* onnxruntime exports */
extern SEXP _superassp_ort_available_cpp(void);
extern SEXP _superassp_ort_version_cpp(void);
extern SEXP _superassp_ort_lib_path_cpp(void);
extern SEXP _superassp_ort_set_lib_dir_cpp(SEXP);
extern SEXP _superassp_ort_create_session_cpp(SEXP, SEXP);
extern SEXP _superassp_ort_run_cpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _superassp_ort_session_input_info_cpp(SEXP);
extern SEXP _superassp_ort_session_output_info_cpp(SEXP);

/* Snack pitch + formant exports */
extern SEXP _superassp_snackp_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _superassp_snackf_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* CREPE inference export */
extern SEXP _superassp_crepe_inference_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* SRH Variant pitch exports */
extern SEXP _superassp_resample_polyphase_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _superassp_srh_core_cpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP _superassp_srh_variant_cpp(SEXP, SEXP, SEXP);
extern SEXP _superassp_srh_variant_debug_cpp(SEXP, SEXP, SEXP);

/* .External calls */
extern SEXP getDObj2(SEXP);
extern SEXP performAssp(SEXP);
extern SEXP performAsspMemory(SEXP);
extern SEXP debugDObjConversion(SEXP);
extern SEXP debugDObjFromFile(SEXP);
//extern SEXP performAssp2(SEXP);


static const R_CallMethodDef CallEntries[] = {
  {"AsspLpTypes_",     (DL_FUNC) &AsspLpTypes_,     0},
  {"AsspSpectTypes_",  (DL_FUNC) &AsspSpectTypes_,  0},
  {"AsspWindowTypes_", (DL_FUNC) &AsspWindowTypes_, 0},
  {"writeDObj_",       (DL_FUNC) &writeDObj_,       2},
  /* Rcpp exports */
  {"_superassp_fast_file_ext",                  (DL_FUNC) &_superassp_fast_file_ext,                  1},
  {"_superassp_fast_is_native",                 (DL_FUNC) &_superassp_fast_is_native,                 2},
  {"_superassp_fast_is_lossless",               (DL_FUNC) &_superassp_fast_is_lossless,               2},
  {"_superassp_fast_strip_file_protocol",       (DL_FUNC) &_superassp_fast_strip_file_protocol,       1},
  {"_superassp_fast_recycle_times",             (DL_FUNC) &_superassp_fast_recycle_times,             2},
  {"_superassp_fast_validate_times",            (DL_FUNC) &_superassp_fast_validate_times,            3},
  {"_superassp_fast_rename_tracks",             (DL_FUNC) &_superassp_fast_rename_tracks,             2},
  {"_superassp_fast_generate_output_paths",     (DL_FUNC) &_superassp_fast_generate_output_paths,     4},
  {"_superassp_fast_calculate_conversion_times",(DL_FUNC) &_superassp_fast_calculate_conversion_times,5},
  {"_superassp_fast_basename",                  (DL_FUNC) &_superassp_fast_basename,                  1},
  {"_superassp_fast_file_path_sans_ext",        (DL_FUNC) &_superassp_fast_file_path_sans_ext,        1},
  {"_superassp_fast_build_conversion_df",       (DL_FUNC) &_superassp_fast_build_conversion_df,       5},
  {"_superassp_fast_check_file_formats",        (DL_FUNC) &_superassp_fast_check_file_formats,        3},
  {"_superassp_estk_pda_cpp",                   (DL_FUNC) &_superassp_estk_pda_cpp,                  13},
  {"_superassp_estk_pitchmark_cpp",             (DL_FUNC) &_superassp_estk_pitchmark_cpp,            15},
  {"_superassp_opensmile_gemaps_cpp",           (DL_FUNC) &_superassp_opensmile_gemaps_cpp,           3},
  {"_superassp_opensmile_extract_cpp",          (DL_FUNC) &_superassp_opensmile_extract_cpp,          4},
  {"_superassp_rapt_cpp",                       (DL_FUNC) &_superassp_rapt_cpp,                       6},
  {"_superassp_swipe_cpp",                      (DL_FUNC) &_superassp_swipe_cpp,                      6},
  {"_superassp_reaper_cpp",                     (DL_FUNC) &_superassp_reaper_cpp,                     6},
  {"_superassp_dio_cpp",                        (DL_FUNC) &_superassp_dio_cpp,                        6},
  {"_superassp_harvest_cpp",                    (DL_FUNC) &_superassp_harvest_cpp,                    6},
  {"_superassp_d4c_cpp",                        (DL_FUNC) &_superassp_d4c_cpp,                        7},
  {"_superassp_sptk_mfcc_cpp",                  (DL_FUNC) &_superassp_sptk_mfcc_cpp,                 10},
  {"_superassp_yin_cpp",                        (DL_FUNC) &_superassp_yin_cpp,                        7},
  {"_superassp_pyin_cpp",                       (DL_FUNC) &_superassp_pyin_cpp,                       7},
  {"_superassp_momel_c",                        (DL_FUNC) &_superassp_momel_c,                        8},
  {"_superassp_tandem_pitch_cpp",               (DL_FUNC) &_superassp_tandem_pitch_cpp,               5},
  {"_superassp_iaif_cpp",                       (DL_FUNC) &_superassp_iaif_cpp,                       6},
  {"_superassp_extract_vq_params_cpp",          (DL_FUNC) &_superassp_extract_vq_params_cpp,          5},
  {"_superassp_gfmiaif_cpp",                    (DL_FUNC) &_superassp_gfmiaif_cpp,                    8},
  /* onnxruntime exports */
  {"_superassp_ort_available_cpp",              (DL_FUNC) &_superassp_ort_available_cpp,              0},
  {"_superassp_ort_version_cpp",                (DL_FUNC) &_superassp_ort_version_cpp,                0},
  {"_superassp_ort_lib_path_cpp",               (DL_FUNC) &_superassp_ort_lib_path_cpp,               0},
  {"_superassp_ort_set_lib_dir_cpp",            (DL_FUNC) &_superassp_ort_set_lib_dir_cpp,            1},
  {"_superassp_ort_create_session_cpp",         (DL_FUNC) &_superassp_ort_create_session_cpp,         2},
  {"_superassp_ort_run_cpp",                    (DL_FUNC) &_superassp_ort_run_cpp,                    5},
  {"_superassp_ort_session_input_info_cpp",     (DL_FUNC) &_superassp_ort_session_input_info_cpp,     1},
  {"_superassp_ort_session_output_info_cpp",    (DL_FUNC) &_superassp_ort_session_output_info_cpp,    1},
  /* Snack pitch + formant */
  {"_superassp_snackp_cpp",                     (DL_FUNC) &_superassp_snackp_cpp,                     6},
  {"_superassp_snackf_cpp",                     (DL_FUNC) &_superassp_snackf_cpp,                    11},
  /* CREPE inference */
  {"_superassp_crepe_inference_cpp",            (DL_FUNC) &_superassp_crepe_inference_cpp,            6},
  /* SRH Variant pitch */
  {"_superassp_resample_polyphase_cpp",         (DL_FUNC) &_superassp_resample_polyphase_cpp,         6},
  {"_superassp_srh_core_cpp",                   (DL_FUNC) &_superassp_srh_core_cpp,                   4},
  {"_superassp_srh_variant_cpp",                (DL_FUNC) &_superassp_srh_variant_cpp,                3},
  {"_superassp_srh_variant_debug_cpp",          (DL_FUNC) &_superassp_srh_variant_debug_cpp,          3},

  {NULL, NULL, 0}
};

static const R_ExternalMethodDef ExternalEntries[] = {
  {"getDObj2",          (DL_FUNC) &getDObj2,          4},
  {"performAssp",       (DL_FUNC) &performAssp,      -1}, // -1 specifies a variable number of arguments
  {"performAsspMemory", (DL_FUNC) &performAsspMemory,-1}, // -1 specifies a variable number of arguments
  {"debugDObjConversion", (DL_FUNC) &debugDObjConversion, -1},
  {"debugDObjFromFile", (DL_FUNC) &debugDObjFromFile, -1},
  //{"performAssp2", (DL_FUNC) &performAssp2,  -1}, // -1 specifies a variable number of arguments
  //{"performAssp", (DL_FUNC) &performAssp,  8},
  //{"performAssp", (DL_FUNC) &performAssp, 11},
  //{"performAssp", (DL_FUNC) &performAssp, 12},
  //{"performAssp", (DL_FUNC) &performAssp, 14},
  //{"performAssp", (DL_FUNC) &performAssp, 15},
  //{"performAssp", (DL_FUNC) &performAssp, 16},
  //{"performAssp", (DL_FUNC) &performAssp, 19},
  {NULL, NULL, 0}
};

void R_init_superassp(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, ExternalEntries);
  R_useDynamicSymbols(dll, FALSE);
}
