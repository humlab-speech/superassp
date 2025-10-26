#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _superassp_d4c_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _superassp_dio_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _superassp_estk_pda_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _superassp_estk_pitchmark_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _superassp_fast_basename(SEXP);
extern SEXP _superassp_fast_build_conversion_df(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _superassp_fast_calculate_conversion_times(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _superassp_fast_check_file_formats(SEXP, SEXP, SEXP);
extern SEXP _superassp_fast_file_ext(SEXP);
extern SEXP _superassp_fast_file_path_sans_ext(SEXP);
extern SEXP _superassp_fast_generate_output_paths(SEXP, SEXP, SEXP, SEXP);
extern SEXP _superassp_fast_is_lossless(SEXP, SEXP);
extern SEXP _superassp_fast_is_native(SEXP, SEXP);
extern SEXP _superassp_fast_recycle_times(SEXP, SEXP);
extern SEXP _superassp_fast_rename_tracks(SEXP, SEXP);
extern SEXP _superassp_fast_strip_file_protocol(SEXP);
extern SEXP _superassp_fast_validate_times(SEXP, SEXP, SEXP);
extern SEXP _superassp_harvest_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _superassp_opensmile_gemaps_cpp(SEXP, SEXP, SEXP);
extern SEXP _superassp_rapt_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _superassp_reaper_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _superassp_sptk_mfcc_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _superassp_swipe_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP AsspLpTypes_(void);
extern SEXP AsspSpectTypes_(void);
extern SEXP AsspWindowTypes_(void);
extern SEXP writeDObj_(SEXP, SEXP);

/* .External calls */
extern SEXP getDObj2(SEXP);
extern SEXP performAsspMemory(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_superassp_d4c_cpp",                         (DL_FUNC) &_superassp_d4c_cpp,                          7},
    {"_superassp_dio_cpp",                         (DL_FUNC) &_superassp_dio_cpp,                          6},
    {"_superassp_estk_pda_cpp",                    (DL_FUNC) &_superassp_estk_pda_cpp,                    13},
    {"_superassp_estk_pitchmark_cpp",              (DL_FUNC) &_superassp_estk_pitchmark_cpp,              15},
    {"_superassp_fast_basename",                   (DL_FUNC) &_superassp_fast_basename,                    1},
    {"_superassp_fast_build_conversion_df",        (DL_FUNC) &_superassp_fast_build_conversion_df,         5},
    {"_superassp_fast_calculate_conversion_times", (DL_FUNC) &_superassp_fast_calculate_conversion_times,  5},
    {"_superassp_fast_check_file_formats",         (DL_FUNC) &_superassp_fast_check_file_formats,          3},
    {"_superassp_fast_file_ext",                   (DL_FUNC) &_superassp_fast_file_ext,                    1},
    {"_superassp_fast_file_path_sans_ext",         (DL_FUNC) &_superassp_fast_file_path_sans_ext,          1},
    {"_superassp_fast_generate_output_paths",      (DL_FUNC) &_superassp_fast_generate_output_paths,       4},
    {"_superassp_fast_is_lossless",                (DL_FUNC) &_superassp_fast_is_lossless,                 2},
    {"_superassp_fast_is_native",                  (DL_FUNC) &_superassp_fast_is_native,                   2},
    {"_superassp_fast_recycle_times",              (DL_FUNC) &_superassp_fast_recycle_times,               2},
    {"_superassp_fast_rename_tracks",              (DL_FUNC) &_superassp_fast_rename_tracks,               2},
    {"_superassp_fast_strip_file_protocol",        (DL_FUNC) &_superassp_fast_strip_file_protocol,         1},
    {"_superassp_fast_validate_times",             (DL_FUNC) &_superassp_fast_validate_times,              3},
    {"_superassp_harvest_cpp",                     (DL_FUNC) &_superassp_harvest_cpp,                      6},
    {"_superassp_opensmile_gemaps_cpp",            (DL_FUNC) &_superassp_opensmile_gemaps_cpp,             3},
    {"_superassp_rapt_cpp",                        (DL_FUNC) &_superassp_rapt_cpp,                         6},
    {"_superassp_reaper_cpp",                      (DL_FUNC) &_superassp_reaper_cpp,                       6},
    {"_superassp_sptk_mfcc_cpp",                   (DL_FUNC) &_superassp_sptk_mfcc_cpp,                   10},
    {"_superassp_swipe_cpp",                       (DL_FUNC) &_superassp_swipe_cpp,                        6},
    {"AsspLpTypes_",                               (DL_FUNC) &AsspLpTypes_,                                0},
    {"AsspSpectTypes_",                            (DL_FUNC) &AsspSpectTypes_,                             0},
    {"AsspWindowTypes_",                           (DL_FUNC) &AsspWindowTypes_,                            0},
    {"writeDObj_",                                 (DL_FUNC) &writeDObj_,                                  2},
    {NULL, NULL, 0}
};

static const R_ExternalMethodDef ExternalEntries[] = {
    {"getDObj2",          (DL_FUNC) &getDObj2,           4},
    {"performAsspMemory", (DL_FUNC) &performAsspMemory, -1},
    {NULL, NULL, 0}
};

void R_init_superassp(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, ExternalEntries);
    R_useDynamicSymbols(dll, FALSE);
}
