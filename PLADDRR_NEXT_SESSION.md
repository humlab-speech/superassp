# Pladdrr Integration - Session 4 Quick Start

**Current State**: Branch `pladdrr-integration`, **Batch 1 COMPLETE** (3/18 = 17%)  
**Last Commit**: `e3d3d47` - Batch 1 complete, formant bug verified fixed  
**Achievement**: ✅ All simple track functions migrated (intensity, pitch, formants)

## Quick Status Check

```bash
cd /Users/frkkan96/Documents/src/superassp
git log --oneline -5
git status

# Test completed functions
R
devtools::load_all()
test_file <- system.file("samples/sustained/a1.wav", package = "superassp")

# All should work:
trk_intensityp(test_file, toFile = FALSE)
trk_pitchp(test_file, toFile = FALSE)
trk_formantp(test_file, toFile = FALSE)
```

## Batch 2: Summary Functions (START HERE)

### Function 1: lst_voice_reportp (RECOMMENDED START)

**Why start here**: Simplest summary function, good template for others

**Source**: `../plabench/R_implementations/voice_report.R` (~320 lines)  
**File**: `R/list_python_pm_pvoice_report.R`  
**Priority**: MEDIUM  
**Estimated time**: 1-2 hours

**What it does**: Basic voice quality measures
- Jitter (local, absolute, rap, ppq5)
- Shimmer (local, absolute, apq3, apq5, apq11)
- HNR (harmonics-to-noise ratio)
- Voicing statistics (fraction voiced, breaks, etc.)

**Key pladdrr calls**:
```r
sound <- pladdrr::Sound(file)
pitch <- sound$to_pitch(...)
point_process <- pitch$to_point_process()

# Jitter/shimmer
jitter_local <- point_process$get_jitter_local(...)
shimmer_local <- point_process$get_shimmer_local(...)

# HNR
hnr <- sound$to_harmonicity_cc(...)
mean_hnr <- hnr$get_mean(0, 0)
```

**Output**: Data frame with columns:
- jitter_local, jitter_abs, jitter_rap, jitter_ppq5
- shimmer_local, shimmer_abs, shimmer_apq3, shimmer_apq5, shimmer_apq11
- hnr_mean, hnr_std
- fraction_voiced, num_breaks, degree_breaks

### Function 2: lst_dsip

**Source**: `../plabench/R_implementations/dsi.R` (~270 lines)  
**File**: `R/list_python_pm_pdsi.R`  
**Priority**: HIGH  
**Estimated time**: 2-3 hours

**What it does**: Dysphonia Severity Index (Wuyts et al., 2000)

**DSI Formula**:
```
DSI = 1.127 + 0.164*MPT - 0.038*I-low + 0.0053*F0-high - 5.30*Jitter
```

**Components**:
1. MPT (Maximum Phonation Time) - from file durations
2. I-low (Minimum Intensity) - from intensity analysis
3. F0-high (Maximum F0) - from pitch analysis
4. Jitter PPQ5 - from voice quality

**Note**: plabench version has "ultra" optimizations - may need simplification

**Key functions needed**:
- `calculate_mpt()` - max duration across files
- `calculate_minimum_intensity()` - voiced intensity minimum
- `calculate_maximum_f0()` - highest F0 value
- `calculate_jitter_ppq5()` - jitter from sustained vowel

### Function 3: lst_voice_tremorp

**Source**: `../plabench/R_implementations/tremor.R` (~700 lines)  
**File**: `R/list_python_pm_pvoice_tremor.R`  
**Priority**: MEDIUM  
**Estimated time**: 2-3 hours

**What it does**: Voice tremor analysis (18 measures)
- Pitch tremor (frequency, amplitude)
- Amplitude tremor (rate, intensity)
- Power spectral density analysis

**Complex** - many FFT operations on extracted tracks

### Function 4: lst_avqip

**Source**: `../plabench/R_implementations/avqi.R` (~400 lines)  
**File**: `R/list_python_pm_pavqi.R`  
**Priority**: HIGH  
**Estimated time**: 3-4 hours

**What it does**: Acoustic Voice Quality Index (AVQI v2.03 & v3.01)

**Most complex** - multiple features combined:
- CPPS (Cepstral Peak Prominence Smoothed)
- HNR
- Shimmer
- Slope of LTAS
- Tilt of regression line
- Dysphonia measures

**Note**: 3x faster than Python version (benchmark from plabench)

## Implementation Pattern for Summary Functions

```r
lst_function <- function(listOfFiles,
                         beginTime = 0.0,
                         endTime = 0.0,
                         ...,
                         verbose = TRUE) {
  
  # Check pladdrr
  if (!pladdrr_available()) {
    stop("pladdrr not available")
  }
  
  # Validate files
  if (!all(file.exists(listOfFiles))) {
    stop("Some files don't exist")
  }
  
  # Progress bar
  if (verbose && length(listOfFiles) > 1) {
    pb <- txtProgressBar(max = length(listOfFiles), style = 3)
  }
  
  results <- list()
  
  for (i in seq_along(listOfFiles)) {
    file <- listOfFiles[i]
    
    # Load with pladdrr
    sound <- pladdrr::Sound(file)
    
    # Time windowing if needed
    if (beginTime > 0 || endTime > 0) {
      duration <- sound$.cpp$duration
      actual_end <- if (endTime > 0) endTime else duration
      sound <- sound$extract_part(beginTime, actual_end, preserve_times = TRUE)
    }
    
    # Perform analyses
    pitch <- sound$to_pitch(...)
    intensity <- sound$to_intensity(...)
    point_process <- pitch$to_point_process()
    
    # Extract measures
    measure1 <- pitch$get_mean(0, 0, "Hertz")
    measure2 <- point_process$get_jitter_local(...)
    
    # Store results
    results[[i]] <- data.frame(
      file = basename(file),
      measure1 = measure1,
      measure2 = measure2
    )
    
    if (verbose && length(listOfFiles) > 1) {
      setTxtProgressBar(pb, i)
    }
  }
  
  if (verbose && length(listOfFiles) > 1) {
    close(pb)
  }
  
  # Combine to data frame
  do.call(rbind, results)
}
```

## Key Differences: Track vs Summary Functions

| Aspect | Track Functions | Summary Functions |
|--------|----------------|-------------------|
| Output | AsspDataObj (SSFF) | data.frame or list |
| Processing | Frame-by-frame | Entire file |
| File output | SSFF files | CSV or R objects |
| Complexity | Extract + format | Multiple analyses + combine |
| Testing | Check tracks | Check statistics |

## Testing Strategy for Summary Functions

```r
test_that("lst_function works", {
  skip_if(!pladdrr_available())
  
  test_file <- system.file("samples/sustained/a1.wav", package = "superassp")
  
  # Single file
  result <- lst_function(test_file, verbose = FALSE)
  
  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) == 1)
  expect_true("measure1" %in% names(result))
  
  # Check values in reasonable range
  expect_gt(result$measure1, 0)
  expect_lt(result$measure1, 1000)  # Adjust for measure type
})
```

## Timeline for Batch 2

- **Function 1** (voice_report): 1-2 hours (simple)
- **Function 2** (dsi): 2-3 hours (moderate, needs simplification)
- **Function 3** (tremor): 2-3 hours (complex)
- **Function 4** (avqi): 3-4 hours (most complex)

**Total**: 8-12 hours (1-2 days if focused)

## After Batch 2

Move to **Batch 3** (complex tracks):
1. `trk_praatsaucep` - Comprehensive voice quality (F0, formants, harmonics, HNR, CPP)
2. `trk_spectral_momentsp` - Spectral moments analysis

Then **Batch 4** (new functions from plabench).

## Current Project Structure

```
R/
├── install_pladdrr.R              # ✅ Setup complete
├── pladdrr_helpers.R              # ✅ Helpers complete
├── ssff_python_pm_pintensity.R   # ✅ Track function 1
├── ssff_python_pm_ppitch.R        # ✅ Track function 2
├── ssff_python_pm_pformantb.R    # ✅ Track function 3
├── list_python_pm_pvoice_report.R # ⏳ Next
├── list_python_pm_pdsi.R          # 📋 TODO
├── list_python_pm_pvoice_tremor.R # 📋 TODO
└── list_python_pm_pavqi.R         # 📋 TODO
```

## Important Notes

1. **Formant bug is FIXED** ✅ - v4.8.16 values are correct
2. **Spectrogram crashes** ⚠️ - Avoid `to_spectrogram()` + `get_power_at()`
3. **Ultra API functions** - plabench uses advanced optimizations, may need simplification
4. **Progress tracking** - Always add progress bars for batch operations
5. **Error handling** - Wrap pladdrr calls in tryCatch where appropriate

## References

- **Completed template**: `R/ssff_python_pm_ppitch.R` (good pattern for track functions)
- **Summary function source**: `../plabench/R_implementations/voice_report.R` (simplest example)
- **pladdrr docs**: `?pladdrr::PointProcess`, `?pladdrr::Harmonicity`
- **Session summary**: `PLADDRR_SESSION_3_SUMMARY.md`

---

**Start with**: Read `../plabench/R_implementations/voice_report.R` and create `R/list_python_pm_pvoice_report.R`.

**Goal**: Complete Batch 2 (4 summary functions) in 1-2 days.
