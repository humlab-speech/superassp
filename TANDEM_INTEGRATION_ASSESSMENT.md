# TANDEM Integration Assessment for superassp

**Date**: 2025-11-06  
**Library**: TANDEM (Tandem Algorithm for pitch estimation and voiced speech segregation)  
**Authors**: G. Hu & D. L. Wang (Ohio State University)  
**Status**: ✅ **ACOUSTIC-ONLY ANALYSIS - APPROVED FOR INTEGRATION**

---

## Executive Summary

TANDEM is a **pitch estimation and voiced speech segregation** algorithm that works with **monaural acoustic audio signals**. It uses computational auditory scene analysis (CASA) principles to separate voiced speech from background noise and estimate pitch simultaneously.

**Integration Recommendation**: ⭐⭐⭐⭐ **RECOMMENDED**

This is appropriate for superassp as it:
1. ✅ Works with standard microphone recordings (acoustic-only)
2. ✅ Provides unique voiced speech separation capability
3. ✅ Includes robust pitch tracking in noisy conditions
4. ✅ Pure C++ implementation (no external dependencies)
5. ✅ Well-documented with published papers

---

## What is TANDEM?

### Algorithm Purpose

**TANDEM = Pitch estimation + Voiced speech segregation**

- **Input**: Monaural audio signal (single-channel acoustic recording)
- **Output**: 
  - Pitch contours (F0 tracking)
  - Binary voiced speech masks (time-frequency segregation)
  - Separated voiced speech streams

### Key Features

1. **Gammatone Filterbank** (64 or 128 channels)
   - Auditory peripheral model
   - Frequency range: 50-8000 Hz
   - Sample rate: 20 kHz

2. **Neural Network-based Pitch Tracking**
   - 3 trained MLPs for pitch estimation
   - Handles multiple simultaneous pitch sources
   - Robust to noise and interference

3. **Voiced Speech Segregation**
   - Time-frequency masking
   - Segment-based grouping
   - Cross-channel correlation

### Published Papers

**128-channel version**:
- Hu, G., & Wang, D. L. (2010). "A tandem algorithm for pitch estimation and voiced speech segregation." *IEEE Trans. Audio, Speech, Lang. Process.*, 18(8), 2067-2079. DOI: [10.1109/TASL.2010.2041110](https://doi.org/10.1109/TASL.2010.2041110)

**64-channel version**:
- Hu, K., & Wang, D. L. (2011). "Unvoiced speech separation from nonspeech interference via CASA and spectral subtraction." *IEEE Trans. Audio, Speech, Lang. Process.*, 19(6), 1600-1609. DOI: [10.1109/TASL.2010.2094211](https://doi.org/10.1109/TASL.2010.2094211)

**Note**: Both references are available in `inst/REFERENCES.bib` as `hu2010tandem` and `hu2011unvoiced`.

---

## Technical Details

### Code Structure

**Location**: `src/tandem/`

**Two versions available**:
1. `tandem_128/` - 128-channel version (original, more accurate)
2. `tandem_64/` - 64-channel version (faster, still robust)

**Recommendation**: Integrate `tandem_64` (better speed/accuracy trade-off for superassp users)

### Source Files (tandem_64)

| File | Lines | Purpose |
|------|-------|---------|
| `tandem.cpp` | 397 | Main algorithm, file I/O |
| `pitch.cpp` | 519 | Pitch tracking with neural networks |
| `segmentation.cpp` | 1482 | Segment-based grouping |
| `voicedMask.cpp` | 723 | Voiced mask estimation |
| `feature.cpp` | 235 | Feature extraction |
| `gammaTone.cpp` | 124 | Gammatone filterbank |
| `mScaleInten.cpp` | 136 | Multi-scale intensity analysis |
| `filter.cpp` | 209 | Signal filtering |
| `tool.cpp` | 216 | Utility functions (FFT, etc.) |
| **Total** | **4041** | **~4K lines of C++** |

### Neural Network Weights

**Location**: `net/`
- `MLP1.64.dat` - Single-unit pitch network (10 KB)
- `MLP2.64.dat` - Multi-unit pitch network (24 KB)
- `MLP3.64.dat` - Mask pitch network (124 bytes)

**Total**: ~35 KB of trained weights

---

## Integration Strategy

### Approach: Rcpp Wrapper

**Goal**: Create `trk_tandem()` function that:
1. Accepts audio files (any format via av package)
2. Converts to required format (20 kHz mono)
3. Calls TANDEM C++ code
4. Returns pitch track + voiced mask as AsspDataObj

### Implementation Plan

#### Step 1: Core Library Compilation (1-2 hours)

**Add to `src/Makevars`**:

```makefile
# TANDEM 64-channel sources
TANDEM_SOURCES = \
  tandem/tandem_64/tool.cpp \
  tandem/tandem_64/gammaTone.cpp \
  tandem/tandem_64/filter.cpp \
  tandem/tandem_64/feature.cpp \
  tandem/tandem_64/voicedMask.cpp \
  tandem/tandem_64/pitch.cpp \
  tandem/tandem_64/segmentation.cpp \
  tandem/tandem_64/mScaleInten.cpp

# Add to PKG_CPPFLAGS
PKG_CPPFLAGS += -Itandem/tandem_64

# Add to sources
CXX_SOURCES += tandem_wrapper.cpp
```

**Modify TANDEM code**:
- Remove `main()` function from `tandem.cpp`
- Extract core logic into callable functions
- Replace file I/O with in-memory processing

---

#### Step 2: Rcpp Wrapper (2-3 hours)

**File**: `src/tandem_wrapper.cpp`

```cpp
#include <Rcpp.h>
#include "tandem/tandem_64/pitch.h"
#include "tandem/tandem_64/voicedMask.h"
#include "tandem/tandem_64/gammaTone.h"

// [[Rcpp::export]]
Rcpp::List tandem_pitch_cpp(
    Rcpp::NumericVector audio_signal,
    int sample_rate,
    double min_pitch = 50.0,
    double max_pitch = 500.0
) {
    // Convert audio to 20 kHz if needed
    // ... resampling code ...
    
    // Initialize TANDEM
    gammaToneFilterBank *audioPery;
    voicedMask *tGroup;
    initVoicedMask(audioPery, tGroup);
    
    // Process audio (in-memory, no file I/O)
    processAudioSamples(audio_signal.begin(), audio_signal.size(), 
                       audioPery, tGroup);
    
    // Extract pitch contours
    std::vector<double> pitch_values;
    std::vector<double> pitch_times;
    std::vector<double> voicing_prob;
    
    for (int i = 0; i < tGroup->numContour; i++) {
        // Extract pitch track from contour i
        // ...
    }
    
    // Extract voiced mask
    Rcpp::NumericMatrix voiced_mask(64, tGroup->corrLgm.size());
    // ... populate mask ...
    
    // Cleanup
    delete audioPery;
    delete tGroup;
    
    return Rcpp::List::create(
        Rcpp::Named("pitch") = pitch_values,
        Rcpp::Named("time") = pitch_times,
        Rcpp::Named("voicing_prob") = voicing_prob,
        Rcpp::Named("voiced_mask") = voiced_mask,
        Rcpp::Named("sample_rate") = sample_rate
    );
}
```

---

#### Step 3: R Wrapper Function (1-2 hours)

**File**: `R/ssff_cpp_tandem.R`

```r
#' TANDEM Pitch Tracking and Voiced Speech Segregation
#'
#' Estimate pitch and separate voiced speech from background using the TANDEM
#' algorithm (Computational Auditory Scene Analysis).
#'
#' @details
#' TANDEM uses a gammatone filterbank and neural networks to simultaneously:
#' - Track pitch in noisy conditions
#' - Separate voiced speech from background interference
#' - Handle multiple simultaneous speakers/sources
#'
#' The algorithm works best with:
#' - Monaural recordings (single-channel)
#' - Sample rates around 16-20 kHz
#' - Speech with background noise or interference
#'
#' @param listOfFiles Character vector of audio file paths
#' @param minF Numeric, minimum F0 in Hz (default: 50)
#' @param maxF Numeric, maximum F0 in Hz (default: 500)
#' @param target_sample_rate Numeric, internal processing rate (default: 20000 Hz)
#' @param return_mask Logical, return time-frequency voiced mask (default: FALSE)
#' @param toFile Logical, write output to SSFF file (default: FALSE)
#' @param explicitExt Character, output file extension (default: "tnd")
#' @param outputDirectory Character, output directory (default: NULL)
#' @param verbose Logical, print progress (default: TRUE)
#' @param ... Additional arguments
#'
#' @return If toFile=FALSE: AsspDataObj or list with:
#'   - `pitch`: F0 track in Hz
#'   - `voicing_prob`: Voicing probability (0-1)
#'   - `voiced_mask`: Time-frequency mask (if return_mask=TRUE)
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic pitch tracking in noise
#' result <- trk_tandem("noisy_speech.wav")
#' plot(result$pitch, type = "l")
#'
#' # With voiced mask for segregation
#' result <- trk_tandem("noisy_speech.wav", return_mask = TRUE)
#' image(result$voiced_mask)  # Visualize segregation
#'
#' # Batch processing
#' files <- c("speaker1.wav", "speaker2.wav")
#' results <- trk_tandem(files, minF = 80, maxF = 400)
#' }
#'
#' @references
#' Hu, G., & Wang, D. L. (2010). A tandem algorithm for pitch estimation 
#' and voiced speech segregation. IEEE Trans. Audio, Speech, Lang. Process., 
#' 18(8), 2067-2079.
trk_tandem <- function(
  listOfFiles,
  minF = 50,
  maxF = 500,
  target_sample_rate = 20000,
  return_mask = FALSE,
  toFile = FALSE,
  explicitExt = "tnd",
  outputDirectory = NULL,
  verbose = TRUE,
  ...
) {
  # Validate inputs
  if (!is.character(listOfFiles)) {
    stop("listOfFiles must be a character vector")
  }
  
  n_files <- length(listOfFiles)
  results <- vector("list", n_files)
  
  if (verbose) {
    message("Processing ", n_files, " file(s) with TANDEM...")
  }
  
  for (i in seq_along(listOfFiles)) {
    if (verbose && n_files > 1) {
      message("  [", i, "/", n_files, "] ", basename(listOfFiles[i]))
    }
    
    # Load audio via av package
    audio_data <- av::read_audio_bin(
      listOfFiles[i],
      channels = 1  # TANDEM requires mono
    )
    
    orig_sr <- attr(audio_data, "sample_rate")
    
    # Resample to 20 kHz if needed
    if (orig_sr != target_sample_rate) {
      # Simple resampling (TODO: use better method)
      ratio <- target_sample_rate / orig_sr
      audio_data <- signal::resample(as.numeric(audio_data), ratio, 1)
    } else {
      audio_data <- as.numeric(audio_data)
    }
    
    # Call TANDEM C++ wrapper
    tandem_result <- tandem_pitch_cpp(
      audio_signal = audio_data,
      sample_rate = target_sample_rate,
      min_pitch = minF,
      max_pitch = maxF
    )
    
    # Convert to AsspDataObj
    tracks <- list(
      pitch = tandem_result$pitch,
      voicing_prob = tandem_result$voicing_prob
    )
    
    if (return_mask) {
      # Include voiced mask as additional tracks
      # (or return separately)
      tracks$voiced_mask <- tandem_result$voiced_mask
    }
    
    assp_obj <- list(
      pitch = tandem_result$pitch,
      voicing_prob = tandem_result$voicing_prob
    )
    
    # Set attributes
    attr(assp_obj, "sampleRate") <- 100  # Analysis rate
    attr(assp_obj, "startTime") <- 0
    attr(assp_obj, "startRecord") <- 1L
    attr(assp_obj, "endRecord") <- length(tandem_result$pitch)
    attr(assp_obj, "trackFormats") <- c("REAL64", "REAL64")
    class(assp_obj) <- "AsspDataObj"
    
    # Write to file if requested
    if (toFile) {
      base_name <- tools::file_path_sans_ext(basename(listOfFiles[i]))
      out_dir <- if (is.null(outputDirectory)) dirname(listOfFiles[i]) else outputDirectory
      output_path <- file.path(out_dir, paste0(base_name, ".", explicitExt))
      
      wrassp::write.AsspDataObj(assp_obj, output_path)
      results[[i]] <- output_path
    } else {
      results[[i]] <- assp_obj
    }
  }
  
  # Simplify output for single file
  if (n_files == 1) {
    return(results[[1]])
  } else {
    return(results)
  }
}

# Set function attributes
attr(trk_tandem, "ext") <- "tnd"
attr(trk_tandem, "tracks") <- c("pitch", "voicing_prob")
attr(trk_tandem, "outputType") <- "SSFF"
attr(trk_tandem, "nativeFiletypes") <- c("wav")
```

---

#### Step 4: Copy Neural Network Weights (5 minutes)

```r
# Copy network files to inst/
dir.create("inst/tandem_net", recursive = TRUE)
file.copy(
  "src/tandem/tandem_64/net/*.dat",
  "inst/tandem_net/"
)
```

---

#### Step 5: Testing (1-2 hours)

**File**: `tests/testthat/test-tandem.R`

```r
test_that("trk_tandem tracks pitch in clean speech", {
  skip_if_not_installed("superassp")
  
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")
  
  result <- trk_tandem(test_wav, toFile = FALSE, verbose = FALSE)
  
  expect_s3_class(result, "AsspDataObj")
  expect_true("pitch" %in% names(result))
  expect_true("voicing_prob" %in% names(result))
  
  # Pitch should be in reasonable range
  valid_pitch <- result$pitch[!is.na(result$pitch)]
  expect_true(all(valid_pitch >= 50 & valid_pitch <= 500))
})

test_that("trk_tandem handles noisy speech", {
  # Test with speech + noise
  # Verify segregation works
})
```

---

## Key Implementation Challenges

### Challenge 1: File I/O → In-Memory Processing

**Current**: TANDEM reads audio from text file
**Solution**: Modify `readInput()` to accept memory buffer

```cpp
// Replace this:
double* readInput(char *filename, int *numSample)
{
  ifstream is(filename);
  // ...
}

// With this:
double* processAudioBuffer(double *samples, int numSamples)
{
  // Normalize and scale like original
  double sumE = 0;
  for (int i = 0; i < numSamples; i++) {
    sumE += samples[i] * samples[i];
  }
  sumE = sqrt(sumE / numSamples);
  double scale = 1000.0 / sumE;  // 60-dB loudness level
  
  double *scaled = new double[numSamples];
  for (int i = 0; i < numSamples; i++) {
    scaled[i] = samples[i] * scale;
  }
  
  return scaled;
}
```

### Challenge 2: Neural Network File Loading

**Current**: Loads from `net/*.dat` files
**Solution**: Load from `inst/tandem_net/` at runtime

```cpp
// Update path resolution
std::string net_path = R_PACKAGE_DIR + "/tandem_net/";
std::string mlp1_path = net_path + "MLP1.64.dat";
// ...
```

### Challenge 3: Sample Rate Requirement

**Current**: Fixed at 20 kHz
**Solution**: Resample input audio

**Options**:
1. Use R's `signal::resample()` before calling C++
2. Implement resampling in C++ (more complex)
3. Document 20 kHz requirement, warn if different

**Recommendation**: Option 1 (R-based resampling)

---

## Expected Benefits

### 1. Unique Capability

**Voiced Speech Segregation** - No other superassp function does this:
- Separate speech from background noise
- Handle multiple simultaneous speakers
- Time-frequency masking for downstream processing

### 2. Robust Pitch Tracking

**Neural network-based** pitch estimation:
- More robust to noise than autocorrelation methods
- Handles reverberant conditions
- Multiple pitch contour tracking

### 3. Research Applications

- **Speech enhancement**: Use voiced mask to filter noise
- **Multi-speaker analysis**: Track multiple F0 contours
- **Noisy speech**: Analyze degraded recordings

---

## Comparison with Existing Functions

| Feature | TANDEM | Existing |
|---------|--------|----------|
| **Pitch tracking** | Yes | Yes (17 methods) |
| **Noise robustness** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐ |
| **Voiced mask** | ✅ Yes | ❌ No |
| **Multi-pitch** | ✅ Yes | ❌ No |
| **Segregation** | ✅ Yes | ❌ No |
| **Speed** | ⭐⭐⭐ Medium | ⭐⭐⭐⭐ Fast |

**Conclusion**: TANDEM fills gaps (segregation, multi-pitch) while providing competitive pitch tracking.

---

## Timeline Estimate

| Phase | Time | Tasks |
|-------|------|-------|
| **Phase 1** | 1-2 hours | Modify TANDEM code for library use |
| **Phase 2** | 2-3 hours | Create Rcpp wrapper |
| **Phase 3** | 1-2 hours | Create R function |
| **Phase 4** | 5 min | Copy network weights |
| **Phase 5** | 1-2 hours | Testing & validation |
| **Total** | **5-9 hours** | **Complete integration** |

---

## Recommendation

✅ **PROCEED with TANDEM integration**

**Rationale**:
1. ✅ Acoustic-only (meets requirement)
2. ✅ Unique capabilities (voiced segregation, multi-pitch)
3. ✅ Well-documented (2 IEEE papers)
4. ✅ Reasonable code size (~4K lines)
5. ✅ No external dependencies
6. ✅ Fills gaps in superassp

**Priority**: Medium-High (valuable for noisy speech analysis)

---

**Document Version**: 1.0  
**Date**: 2025-11-06  
**Status**: Ready for implementation
