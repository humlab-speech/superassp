# Pladdrr Integration Project - Session 7 Continuation Prompt

## Project Context

We are performing a **complete migration** of the superassp R package from Python's `parselmouth` to R's `pladdrr` for all Praat-based DSP functionality. This eliminates Python dependencies and achieves pure R/C implementation.

**Branch**: `pladdrr-integration` (created from `cpp_optimization`)  
**Package Version**: 0.11.1  
**pladdrr Version Required**: >= 4.8.16  
**Working Directory**: `/Users/frkkan96/Documents/src/superassp/`

## What We've Accomplished So Far

### ✅ BATCH 1 COMPLETE (3/18 functions = 17%)
1. **trk_intensityp** - ✅ DONE (Session 3)
2. **trk_pitchp** - ✅ DONE (Session 3)
3. **trk_formantp** - ✅ DONE (Session 3)

### ✅ BATCH 2 COMPLETE (4/18 functions = 22%)
4. **lst_voice_reportp** - ✅ DONE (Session 4)
5. **lst_dsip** - ✅ DONE (Session 4)
6. **lst_voice_tremorp** - ✅ DONE (Session 5)
7. **lst_avqip** - ✅ DONE (Session 5)

### ✅ BATCH 3 COMPLETE (2/18 functions = 11%)
8. **trk_spectral_momentsp** - ✅ DONE (Session 5)
9. **trk_praatsaucep** - ✅ DONE (Session 6) - **MAJOR MILESTONE**
   - **Most complex function**: 36 output tracks
   - VoiceSauce-compatible voice quality analysis
   - Hawks-Miller bandwidth estimation + Iseli-Alwan formant correction
   - F0, F1-F3, B1-B3, uncorrected/corrected harmonics, CPP, HNR (4 bands)
   - ~680 lines including helper functions

**OVERALL PROGRESS: 50% COMPLETE (9/18 functions)** ✅

## Current State - WHERE WE ARE NOW

**Last Session (6)**: Completed `trk_praatsaucep` - the most complex track function with 36 measures

**Last Commits**:
```
025a092 - feat: trk_praatsaucep (pladdrr) - 36 voice quality tracks
b3f16e5 - docs: Session 5 summary - Batch 2 complete + spectral moments
a94f2cc - docs: Update status - 8/18 complete (44%)
```

**Git status**: All changes committed, working directory clean

**Testing Status**: ⚠️ **NO FUNCTIONS TESTED** (pladdrr v4.8.16+ not installed - user deferred installation)

## What We're Going To Do Next

### START PHASE 4: New Functions from plabench

We've completed all the direct migrations (Batches 1-3). Now we implement **new functions** from plabench that don't exist in superassp yet.

**Target for Session 7**: Implement 2-3 functions from Phase 4

### Priority Order for Phase 4 (9 remaining functions)

**HIGH PRIORITY** (implement first):

1. **trk_cpps** - Cepstral Peak Prominence Smoothed
   - Source: `/Users/frkkan96/Documents/src/plabench/R_implementations/cpp.R`
   - Complexity: MEDIUM (~200 lines expected)
   - Output: SSFF track with CPP values over time
   - Algorithm: Smooth CPP calculation from power cepstrum
   - Extension: `.cpp` or `.cps`

2. **lst_vq** - Voice Quality Measures (summary)
   - Source: `/Users/frkkan96/Documents/src/plabench/R_implementations/vq.R`
   - Complexity: MEDIUM (~300 lines expected)
   - Output: Data frame with voice quality scalars
   - Measures: Jitter, shimmer, HNR, CPP (aggregate statistics)
   - Extension: `.vq` (JSTF format)

3. **trk_vuv** - Voice/Unvoiced Detection
   - Source: `/Users/frkkan96/Documents/src/plabench/R_implementations/vuv.R`
   - Complexity: LOW (~150 lines expected)
   - Output: SSFF track with binary voiced/unvoiced decisions
   - Algorithm: Based on pitch + intensity thresholds
   - Extension: `.vuv`

**MEDIUM PRIORITY** (implement next):

4. **lst_pharyngeal** - Pharyngeal Analysis
   - Source: `/Users/frkkan96/Documents/src/plabench/R_implementations/pharyngeal.R`
   - Complexity: MEDIUM (~250 lines expected)
   - Output: Data frame with pharyngeal measures
   - Extension: `.pha` (JSTF format)

**LOW PRIORITY** (implement if time):

5. **utils_momel** - MOMEL pitch target extraction
   - Source: `/Users/frkkan96/Documents/src/plabench/R_implementations/momel_pure_r.R`
   - Note: Pure R implementation, check if already in dysprosody
   - May not need separate implementation

6. **utils_intsint** - INTSINT pitch coding
   - Source: `/Users/frkkan96/Documents/src/plabench/R_implementations/intsint_pure_r.R`
   - Note: Pure R implementation, check if already in dysprosody
   - May not need separate implementation

**SPECIAL CASES**:

7. **lst_dysprosody** - UPDATE existing function
   - Already exists in superassp (uses Parselmouth)
   - Check plabench version for any improvements
   - May just need Python → pladdrr port (already done?)

8. **Additional plabench utilities** - Review shared code
   - Check for any helper functions we need

9. **trk_formantpathp** - MERGED into trk_formantp
   - Already handled by HMM tracking in trk_formantp
   - No separate implementation needed

## File Locations

### Source Material (plabench)
- **CPP**: `/Users/frkkan96/Documents/src/plabench/R_implementations/cpp.R`
- **VQ**: `/Users/frkkan96/Documents/src/plabench/R_implementations/vq.R`
- **VUV**: `/Users/frkkan96/Documents/src/plabench/R_implementations/vuv.R`
- **Pharyngeal**: `/Users/frkkan96/Documents/src/plabench/R_implementations/pharyngeal.R`
- **MOMEL**: `/Users/frkkan96/Documents/src/plabench/R_implementations/momel_pure_r.R`
- **INTSINT**: `/Users/frkkan96/Documents/src/plabench/R_implementations/intsint_pure_r.R`
- **Dysprosody**: `/Users/frkkan96/Documents/src/plabench/R_implementations/dysprosody.R`

### superassp Target Files (to be created)
- **CPP**: `R/ssff_pladdrr_cpps.R` (new file)
- **VQ**: `R/list_pladdrr_vq.R` (new file)
- **VUV**: `R/ssff_pladdrr_vuv.R` (new file)
- **Pharyngeal**: `R/list_pladdrr_pharyngeal.R` (new file)

### Helper Files (already complete)
- `R/pladdrr_helpers.R` - Core helpers
- `R/install_pladdrr.R` - Installation functions
- `R/jstf_helpers.R` - JSTF writing (FIXED in Session 4)
- `R/json_track_core.R` - JSTF core (FIXED in Session 4)
- `R/json_track_methods.R` - JSTF methods (FIXED in Session 4)

### Documentation
- `PLADDRR_MIGRATION_STATUS.md` - Progress tracker (currently 9/18 = 50%)
- `PLADDRR_IMPLEMENTATION_PLAN.md` - Detailed implementation guide
- `SESSION_6_SUMMARY.md` - Last session summary
- `NEWS.md` - Version 0.11.1 notes

## Key Technical Patterns Established

### For Track Functions (trk_*)

**Function signature**:
```r
trk_function <- function(listOfFiles,
                         beginTime = 0.0,
                         endTime = 0.0,
                         # ... DSP parameters ...
                         toFile = TRUE,
                         explicitExt = "ext",
                         outputDirectory = NULL,
                         verbose = TRUE)
```

**Standard structure** (from trk_praatsaucep):
1. Check pladdrr availability
2. Validate inputs
3. Progress bar for multiple files
4. Loop through files:
   - Load with `av_load_for_pladdrr()`
   - Create pladdrr objects (Pitch, Formant, etc.)
   - Loop through frames/timepoints
   - Extract measures
   - Build AsspDataObj
   - Write to SSFF if `toFile=TRUE`
5. Return results

**AsspDataObj creation**:
```r
assp_obj <- list(track1 = values1, track2 = values2, ...)
attr(assp_obj, "sampleRate") <- sample_rate
attr(assp_obj, "startTime") <- start_time
attr(assp_obj, "startRecord") <- 1L
attr(assp_obj, "endRecord") <- as.integer(num_frames)
attr(assp_obj, "trackFormats") <- rep("REAL32", num_tracks)
attr(assp_obj, "origFreq") <- original_freq
class(assp_obj) <- "AsspDataObj"
```

**Function attributes** (REQUIRED):
```r
attr(trk_function, "ext") <- "ext"
attr(trk_function, "tracks") <- c("track1", "track2", ...)
attr(trk_function, "outputType") <- "SSFF"
attr(trk_function, "nativeFiletypes") <- c("wav")
```

### For Summary Functions (lst_*)

**Function signature**:
```r
lst_function <- function(listOfFiles,
                         beginTime = 0.0,
                         endTime = 0.0,
                         # ... DSP parameters ...
                         toFile = FALSE,
                         explicitExt = "ext",
                         outputDirectory = NULL,
                         verbose = TRUE)
```

**JSTF output** (for lst_* functions):
```r
if (toFile) {
  output_paths <- write_lst_results_to_jstf(
    results = results_list,
    listOfFiles = listOfFiles,
    function_name = "lst_function",
    explicitExt = explicitExt,
    outputDirectory = outputDirectory,
    verbose = verbose
  )
  return(invisible(output_paths))
}
return(results)
```

**Extension registration** (inst/extdata/json_extensions.csv):
```csv
function,extension,description,fields,format
lst_function,ext,Description,N,JSTF
```

### pladdrr API Patterns

**Direct API** (fastest):
- `to_pitch_cc_direct()`, `to_pitch_ac_direct()`
- `to_formant_direct()`
- `to_intensity_direct()`
- `to_harmonicity_direct()`

**Ultra API** (batch/optimized, when available):
- `get_durations_batch()`
- `get_jitter_shimmer_batch()`
- `calculate_f0_stats_ultra()`
- `get_voice_quality_ultra()`
- `calculate_cpps_ultra()`

**R6 methods**:
- `sound$extract_part()`, `sound$filter_pass_hann_band()`
- `sound$to_spectrum()`, `spectrum$to_ltas_1to1()`, `spectrum$to_power_cepstrum()`
- `pitch$get_value_at_time()`, `formant$get_value_at_time()`
- `cepstrum$get_peak_prominence()`

## Timeline & Goals

- **Started**: 2026-02-06 (Session 3)
- **Current**: Session 6 complete (9/18 = 50%)
- **Target**: 2026-02-27 (3 weeks)
- **Current pace**: ~2-3 functions per session
- **Remaining**: 9 functions × ~1 hour = ~9 hours ≈ 3-4 sessions
- **Expected completion**: 2026-02-09 (18 days ahead of schedule!)

**Session 7 goals**: Implement 2-3 new functions (trk_cpps, lst_vq, trk_vuv)

## Quick Start Commands for Session 7

```bash
# Navigate and check status
cd /Users/frkkan96/Documents/src/superassp
git status
git log --oneline -5
git branch

# Check current progress
cat PLADDRR_MIGRATION_STATUS.md | grep "Progress:"
cat SESSION_6_SUMMARY.md | head -30

# Start implementing new functions
# 1. Read plabench source: /Users/frkkan96/Documents/src/plabench/R_implementations/cpp.R
# 2. Create superassp file: R/ssff_pladdrr_cpps.R
# 3. Implement with established patterns
# 4. Add function attributes
# 5. Test structure (don't run - pladdrr not installed)
# 6. Repeat for lst_vq and trk_vuv
```

## Important Notes

1. **pladdrr not installed**: Cannot run tests, focus on implementation
2. **Phase 4 = new functions**: No existing code to replace, create from scratch
3. **Follow established patterns**: Use trk_praatsaucep and lst_avqip as templates
4. **JSTF for lst_* functions**: Use json track format for file output
5. **Extension registration**: Add to inst/extdata/json_extensions.csv for new lst_* functions
6. **Helper functions**: Mark with `.` prefix and `@keywords internal`
7. **Function attributes**: REQUIRED for all trk_* functions

## Success Criteria for Session 7

**trk_cpps** (Cepstral Peak Prominence Smoothed):
- ✅ Main function accepts single/multiple files
- ✅ Returns AsspDataObj with CPP track
- ✅ Time windowing support
- ✅ Batch processing with progress bar
- ✅ SSFF file output support
- ✅ Documented with examples
- ✅ Function attributes set

**lst_vq** (Voice Quality summary):
- ✅ Returns data frame with scalar measures
- ✅ Jitter, shimmer, HNR, CPP included
- ✅ JSTF file output support
- ✅ Batch processing
- ✅ Extension registered in json_extensions.csv

**trk_vuv** (Voice/Unvoiced Detection):
- ✅ Returns AsspDataObj with binary track
- ✅ Threshold-based algorithm
- ✅ SSFF output support
- ✅ Simple and fast

**ESTIMATED TIME**: 2-3 hours for all 3 functions

## Key Differences from Previous Batches

1. **No existing code**: Creating new functions, not replacing Parselmouth versions
2. **Simpler algorithms**: trk_cpps, trk_vuv are less complex than trk_praatsaucep
3. **More autonomy**: Less reference code, more adaptation from plabench
4. **JSTF expertise**: We've fixed all bugs, now apply to new functions

## Motivation

We're **HALFWAY DONE** and **18 days ahead of schedule**! 🎉

The hardest function (trk_praatsaucep with 36 tracks) is complete. The remaining functions are simpler and should go faster.

**LET'S FINISH STRONG!**

---

**START BY**: Reading `/Users/frkkan96/Documents/src/plabench/R_implementations/cpp.R` to understand the CPP algorithm, then create `R/ssff_pladdrr_cpps.R`.
