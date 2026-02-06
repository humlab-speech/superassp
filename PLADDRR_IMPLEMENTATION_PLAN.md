# Pladdrr Migration Implementation Plan

## Summary

**Complete**: Phase 1-2 + trk_intensityp ✅
**Remaining**: 17 functions across phases 3-4

## Implementation Strategy

Given scope, implement in batches:
- **Batch 1 (HIGH)**: Simple track functions (2-3 functions)
- **Batch 2 (HIGH)**: Summary functions (4 functions)  
- **Batch 3 (MEDIUM)**: Complex track functions (3 functions)
- **Batch 4 (NEW)**: New functions from plabench (8 functions)

---

## Batch 1: Simple Track Functions (Next Priority)

### 1. trk_pitchp - `R/ssff_python_pm_ppitch.R`

**Source**: `plabench/R_implementations/pitch.R`  
**Complexity**: MEDIUM  
**Priority**: HIGH

**Key Implementation Points**:
```r
# Use pladdrr direct API
pitch_cc_ptr <- pladdrr::to_pitch_cc_direct(sound, ...)
pitch_ac_ptr <- pladdrr::to_pitch_ac_direct(sound, ...)
pitch_cc <- pladdrr::Pitch(.xptr = pitch_cc_ptr)
pitch_ac <- pladdrr::Pitch(.xptr = pitch_ac_ptr)

# Get data frames
cc_df <- pitch_cc$as_data_frame()
ac_df <- pitch_ac$as_data_frame()

# Create multiple tracks: cc, ac (SPINET/SHS not available)
# Convert to AsspDataObj with tracks: pitch_cc, pitch_ac
```

**Output**: `.pit` file with 2 tracks (cc, ac)

---

### 2. trk_formantp - `R/ssff_python_pm_pformantb.R` 

**Source**: `plabench/R_implementations/formant.R`  
**Complexity**: MEDIUM-HIGH  
**Priority**: HIGH

**Key Implementation Points**:
```r
# Burg's method
formant_ptr <- pladdrr::to_formant_direct(
  sound,
  time_step = time_step,
  max_formants = num_formants,
  max_formant = max_formant_hz,
  window_length = window_length,
  pre_emphasis = pre_emphasis
)
formant <- pladdrr::Formant(.xptr = formant_ptr)

# Optional HMM tracking
if (track_formants) {
  formant <- formant$track(
    num_tracks, ref_f1, ref_f2, ref_f3, ref_f4, ref_f5,
    frequency_cost, bandwidth_cost, transition_cost
  )
}

# Get long format data
formant_df_long <- formant$as_data_frame()

# Convert to wide format: F1, F2, F3, F4, F5, B1, B2, B3, B4, B5
# Use pladdrr_df_to_superassp(formant_df_long, "formant")

# Optional: Add spectral intensity (L1-L5)
if (include_intensity) {
  spectrogram <- sound$to_spectrogram(...)
  for each formant: L_vals[i] <- spectrogram$get_power_at(time, freq)
}
```

**Critical**: Test formant bug fix in v4.8.16!  
**Output**: `.pfm` file with F1-F5, B1-B5, optionally L1-L5

---

### 3. Merge formantpath into formantp

**Source**: `R/ssff_python_pm_pformantpathb.R`  
**Strategy**: Add `track_formants` parameter to trk_formantp  
**Remove**: Separate formantpath function

---

## Batch 2: Summary Functions

### 4. lst_avqip - `R/list_python_pm_pavqi.R`

**Source**: `plabench/R_implementations/avqi.R`  
**Complexity**: HIGH  
**Priority**: CRITICAL (3x faster than Python!)

**Key Implementation Points**:
```r
# Load/concatenate CS and SV files
cs_sound <- load_and_concatenate_sounds(cs_files)
sv_sound <- load_and_concatenate_sounds(sv_files)

# Filter
cs_filtered <- cs_sound$filter_stop_hann_band(0, 34, 0.1)
sv_filtered <- sv_sound$filter_stop_hann_band(0, 34, 0.1)

# Extract voiced (version-dependent)
voiced_cs <- extract_voiced_segments(cs_filtered, version)

# Use pladdrr Ultra API (v4.8.15+)
cpps <- calculate_cpps_ultra(...)
vq <- get_voice_quality_ultra(...)  # HNR + shimmer combined

# LTAS slope/tilt
ltas <- sound$to_ltas(...)
slope <- ltas$get_slope(...)
tilt <- ltas$get_slope(...)  # Different freq ranges

# Calculate AVQI
avqi <- calculate_avqi_formula(cpps, hnr, shimmer, slope, tilt, version)
```

**Output**: List with AVQI score + 6 components  
**Extensions**: `.avq` (JSTF format)

---

### 5. lst_dsip - `R/list_python_pm_pdsi.R`

**Source**: `plabench/R_implementations/dsi.R`  
**Complexity**: MEDIUM  
**Priority**: HIGH

**Components**:
- MPT: Maximum phonation time (from file duration)
- F0-high: Maximum F0 from /a/ vowel
- I-low: Minimum intensity from high pitch
- Jitter ppq5: From sustained vowel

**Output**: DSI = 0.13 × MPT + 0.0053 × F0-high - 0.26 × I-low - 1.18 × (jitter% × 100) + 12.4

---

### 6. lst_voice_tremorp - `R/list_python_pm_pvoice_tremor.R`

**Source**: `plabench/R_implementations/tremor.R`  
**Complexity**: HIGH  
**Priority**: MEDIUM

**18 measures**: Frequency tremor (9) + Amplitude tremor (9)  
**Output**: Named list with all measures

---

### 7. lst_voice_reportp - `R/list_python_pm_pvoice_report.R`

**Source**: `plabench/R_implementations/voice_report.R`  
**Complexity**: MEDIUM  
**Priority**: MEDIUM

**Measures**: Jitter, shimmer, HNR, mean/SD pitch, intensity  
**Output**: Named list

---

## Batch 3: Complex Track Functions

### 8. trk_praatsaucep - `R/ssff_python_pm_psauce.R`

**Source**: `plabench/R_implementations/praatsauce.R`  
**Complexity**: VERY HIGH  
**Priority**: HIGH

**Outputs** (time series):
- F0, formants (F1-F4, B1-B4)
- Harmonics H1, H2, H4 (uncorrected)
- Harmonics H1c, H2c, H4c (Iseli-Alwan corrected)
- HNR at 5 bands (0.5, 1.5, 2.5, 3.5 kHz)
- CPP

**Functions needed**:
- hawks_miller_bandwidth()
- iseli_alwan_correction()

**Output**: `.psa` file with 20+ tracks

---

### 9. trk_spectral_momentsp - `R/ssff_python_pm_pspectral_moments.R`

**Source**: `plabench/R_implementations/spectral_moments.R`  
**Complexity**: MEDIUM  
**Priority**: LOW

**Output**: Spectral moments (center of gravity, SD, skewness, kurtosis)

---

## Batch 4: New Functions

### 10. trk_cpps - NEW from `plabench/R_implementations/cpp.R`

**Function**: Cepstral Peak Prominence (Smoothed)  
**Use**: Base for AVQI and voice quality

```r
trk_cpps <- function(listOfFiles, ..., per_frame = FALSE) {
  # PowerCepstrogram
  cepstrogram <- sound$to_power_cepstrogram(...)
  
  if (per_frame) {
    # Return time series
  } else {
    # Return single CPPS value
    cpps <- calculate_cpps_fast(sound, ...)
  }
}
```

---

### 11. trk_vuv - NEW from `plabench/R_implementations/vuv.R`

**Function**: Voice/Unvoiced detection  
**Output**: Binary classification per frame

---

### 12. lst_vq - NEW from `plabench/R_implementations/vq.R`

**Function**: Voice quality measures  
**Expanded set beyond voice_report**

---

### 13. lst_pharyngeal - NEW from `plabench/R_implementations/pharyngeal.R`

**Function**: Pharyngeal analysis  
**Specialized acoustic measures**

---

### 14-16. Prosody Utils

Check if already present in superassp:
- momel_pure_r.R
- intsint_pure_r.R  
- Update lst_dysprosody if needed

---

### 17-18. Additional helpers

Port shared utilities from plabench

---

## Testing Strategy

For each function:

```r
test_that("function works with pladdrr", {
  test_file <- system.file("samples/sustained/a1.wav", package = "superassp")
  skip_if(test_file == "")
  skip_if(!pladdrr_available())
  
  # Single file, toFile=FALSE
  result <- trk_function(test_file, toFile = FALSE)
  expect_s3_class(result, "AsspDataObj")
  
  # Multiple files, toFile=TRUE
  files <- rep(test_file, 2)
  n <- trk_function(files, toFile = TRUE, outputDirectory = tempdir())
  expect_equal(n, 2)
  
  # Time windowing
  result <- trk_function(test_file, beginTime = 0.1, endTime = 0.5, toFile = FALSE)
  expect_true(nrow(result$track) > 0)
  
  # Verify SSFF format
  expect_true("sampleRate" %in% names(attributes(result)))
})
```

---

## Performance Monitoring

Track execution time for each function:

```r
# Benchmark template
library(microbenchmark)

mb <- microbenchmark(
  pladdrr = trk_function(test_file, toFile = FALSE),
  times = 10
)

print(summary(mb))
```

Document in `PLADDRR_PERFORMANCE_ANALYSIS.md`

---

## Next Steps

1. **Immediate**: Implement Batch 1 (pitch, formant)
2. **Day 2**: Implement Batch 2 (AVQI, DSI, tremor, voice_report)
3. **Day 3**: Implement Batch 3 (praatsauce, spectral_moments)
4. **Day 4-5**: Implement Batch 4 (new functions)
5. **Day 6**: Testing
6. **Day 7**: Performance analysis
7. **Day 8**: Documentation
8. **Day 9**: Final polish

---

## Reference Patterns

### Standard trk_* function structure:

```r
trk_function <- function(listOfFiles, beginTime = 0.0, endTime = 0.0,
                         ..., toFile = TRUE, explicitExt = "ext",
                         outputDirectory = NULL, verbose = TRUE) {
  
  # 1. Validation
  if (!pladdrr_available()) stop("pladdrr required")
  if (length(listOfFiles) > 1 && !toFile) stop("toFile=FALSE for single files only")
  
  # 2. Setup file processing
  fileBeginEnd <- data.frame(listOfFiles, beginTime, endTime)
  
  # 3. Loop over files
  for (i in 1:nrow(fileBeginEnd)) {
    # Load via av
    sound <- av_load_for_pladdrr(file, start_time, end_time)
    
    # Process with pladdrr
    result_ptr <- pladdrr::to_xxx_direct(sound, ...)
    result_obj <- pladdrr::XXX(.xptr = result_ptr)
    
    # Extract data
    df <- result_obj$as_data_frame()
    
    # Convert to AsspDataObj
    outDataObj <- create_assp_data_obj(df, ...)
    
    # Write if toFile=TRUE
    if (toFile) {
      wrassp::write.AsspDataObj(outDataObj, file)
    }
  }
  
  # 4. Return
  if (toFile) return(n_files) else return(outDataObj)
}

# Set attributes
attr(trk_function, "ext") <- "ext"
attr(trk_function, "tracks") <- c("track1", "track2")
attr(trk_function, "outputType") <- "SSFF"
attr(trk_function, "nativeFiletypes") <- c("wav")
```

### Standard lst_* function structure:

```r
lst_function <- function(listOfFiles, beginTime = 0.0, endTime = 0.0,
                         ..., toFile = FALSE, explicitExt = "ext",
                         outputDirectory = NULL, verbose = TRUE) {
  
  # Similar validation and loop
  
  # Compute summary measures (not time series)
  measures <- list(
    measure1 = value1,
    measure2 = value2,
    ...
  )
  
  # If toFile, write as JSTF
  if (toFile) {
    json_obj <- create_json_track_obj(measures, ...)
    write_json_track(json_obj, output_path)
    return(invisible(output_path))
  }
  
  return(measures)
}
```

---

## Critical Success Factors

1. **Formant bug verification**: Test v4.8.16 fix immediately
2. **AVQI performance**: Confirm 3x speedup maintained
3. **SSFF compatibility**: Ensure emuR can read output
4. **Track naming**: Consistent with existing conventions
5. **Error handling**: Graceful fallbacks, clear messages
6. **Documentation**: Every function fully documented

---

## Risk Mitigation

- **Test early, test often**: Don't batch-implement without testing
- **Commit frequently**: After each function or small batch
- **Keep plabench reference**: For validation and debugging
- **Document issues**: For pladdrr developer report
- **Maintain TODO**: Update PLADDRR_MIGRATION_STATUS.md regularly
