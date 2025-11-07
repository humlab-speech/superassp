# TANDEM Full Integration - Implementation Notes

**Date**: 2025-11-07  
**Status**: 📋 **DEFERRED - Framework Complete, Awaiting Full Implementation**  
**Estimated Time**: 6-8 hours (complex refactoring required)

---

## Executive Summary

After detailed analysis of the TANDEM source code, **full integration is more complex than initially estimated** and should be done carefully to avoid introducing bugs. The **framework is complete and ready** - the function `trk_tandem()` exists and can be used now with placeholder results.

**Recommendation**: Defer full TANDEM core integration until:
1. Time is available for careful refactoring (6-8 hours)
2. Test suite is expanded to validate accuracy
3. Performance benchmarks are in place

---

## Why Deferral is Prudent

### Complexity Discovered

1. **File I/O deeply embedded**: TANDEM uses `readInput()` and file-based processing throughout
2. **Global state management**: Multiple initialization functions with complex interdependencies
3. **Neural network loading**: MLP weights loaded from files at initialization
4. **Memory management**: Complex object lifecycle (new/delete in C++)
5. **No clear API boundary**: Code designed as standalone executable, not library

### Risk Assessment

| Risk | Severity | Impact |
|------|----------|--------|
| Memory leaks | High | Package crashes, R session instability |
| Incorrect results | High | Invalid pitch tracking data |
| Build failures | Medium | Package won't compile |
| Platform issues | Medium | May work on Mac, fail on Linux/Windows |

### Current Status

✅ **What's Working**:
- Function interface complete (`trk_tandem()`)
- Documentation complete
- Build system integrated
- References added to bibliography
- Placeholder implementation functional

⏸️ **What's Pending**:
- Replace `readInput()` with memory-based processing
- Extract pitch contours from TANDEM data structures
- Load neural networks from `inst/tandem_net/`
- Validate accuracy against reference implementation

---

## Integration Plan (When Time Allows)

### Phase 1: Refactor TANDEM Source (3-4 hours)

**Goal**: Make TANDEM code usable as a library

#### Task 1.1: Create Memory-Based Input Function

**File**: `src/tandem/tandem_64/tandem.cpp`

```cpp
// Replace readInput() with:
double* processAudioBuffer(double *samples, int numSamples, double *scale_out) {
    // Calculate RMS for scaling
    double sumE = 0;
    for (int i = 0; i < numSamples; i++) {
        sumE += samples[i] * samples[i];
    }
    sumE = sqrt(sumE / numSamples);
    double scale = 1000.0 / sumE;  // 60-dB loudness level
    *scale_out = scale;
    
    // Scale audio
    double *signal = new double[numSamples];
    for (int i = 0; i < numSamples; i++) {
        signal[i] = scale * samples[i];
    }
    
    return signal;
}
```

#### Task 1.2: Refactor voicedMaskEst for Memory Input

```cpp
// New signature
void voicedMaskEstMemory(double *audio, int sigLength, 
                        gammaToneFilterBank *AudiPery, 
                        voicedMask *TGroup) {
    double scale;
    double *Input = processAudioBuffer(audio, sigLength, &scale);
    
    TGroup->numFrame = sigLength / (TGroup->fs / 100);
    TGroup->newPitchPara();
    
    AudiPery->sigLength = sigLength;
    for(int chan = 0; chan < TGroup->numberChannel; chan++) {
        AudiPery->filtering(Input, sigLength, chan);
    }
    
    TGroup->computeFeature(AudiPery, sigLength, 1);
    TGroup->dtmPitchMask();
    
    delete [] Input;
}
```

#### Task 1.3: Extract Pitch Contours

```cpp
struct TandemResults {
    std::vector<double> pitch;       // F0 values (Hz)
    std::vector<double> voicing;     // Voicing probability (0-1)
    std::vector<int> pitch_idx;      // Pitch period indices
    int n_frames;
    int n_contours;
};

TandemResults extractPitchContours(voicedMask *TGroup) {
    TandemResults results;
    results.n_frames = TGroup->numFrame;
    results.n_contours = TGroup->numContour;
    
    // Initialize with NaN
    results.pitch.resize(TGroup->numFrame, R_NaReal);
    results.voicing.resize(TGroup->numFrame, 0.0);
    
    // Extract from TANDEM pitch contours
    for (int c = 0; c < TGroup->numContour; c++) {
        pitchContour *pc = &(TGroup->Pitch[c]);
        for (int f = pc->sFrame; f <= pc->eFrame; f++) {
            if (pc->indicate[f - pc->sFrame] > 0) {
                // Convert delay index to Hz
                int delay = pc->value[f - pc->sFrame];
                double freq = (double)TGroup->fs / (double)delay;
                
                results.pitch[f] = freq;
                results.voicing[f] = pc->mProb[f - pc->sFrame].value[0];  // Probability
            }
        }
    }
    
    return results;
}
```

---

### Phase 2: Update C++ Wrapper (2 hours)

**File**: `src/tandem_wrapper.cpp`

```cpp
#include <Rcpp.h>
#include "tandem/tandem_64/pitch.h"
#include "tandem/tandem_64/voicedMask.h"
#include "tandem/tandem_64/gammaTone.h"

// Include new refactored functions
extern "C" {
    void initVoicedMask(gammaToneFilterBank *&AudiPery, voicedMask *&TGroup);
    void voicedMaskEstMemory(double *audio, int sigLength, 
                            gammaToneFilterBank *AudiPery, 
                            voicedMask *TGroup);
    TandemResults extractPitchContours(voicedMask *TGroup);
}

// [[Rcpp::export]]
Rcpp::List tandem_pitch_cpp(
    Rcpp::NumericVector audio_signal,
    int sample_rate = 20000,
    double min_pitch = 50.0,
    double max_pitch = 500.0,
    std::string net_path = ""
) {
    // Validate sample rate
    if (sample_rate != 20000) {
        Rcpp::warning("TANDEM requires 20 kHz. Results may be suboptimal.");
    }
    
    int n_samples = audio_signal.size();
    if (n_samples < 1000) {
        Rcpp::stop("Audio too short (minimum 1000 samples)");
    }
    
    // Initialize TANDEM
    gammaToneFilterBank *audioPery;
    voicedMask *tGroup;
    initVoicedMask(audioPery, tGroup);
    
    // Process audio
    double *samples = &audio_signal[0];
    voicedMaskEstMemory(samples, n_samples, audioPery, tGroup);
    
    // Extract results
    TandemResults results = extractPitchContours(tGroup);
    
    // Convert to R vectors
    Rcpp::NumericVector pitch(results.pitch.begin(), results.pitch.end());
    Rcpp::NumericVector voicing(results.voicing.begin(), results.voicing.end());
    
    // Cleanup
    delete audioPery;
    delete tGroup;
    
    return Rcpp::List::create(
        Rcpp::Named("pitch") = pitch,
        Rcpp::Named("voicing_prob") = voicing,
        Rcpp::Named("sample_rate") = sample_rate,
        Rcpp::Named("n_frames") = results.n_frames,
        Rcpp::Named("n_contours") = results.n_contours,
        Rcpp::Named("status") = "full_tandem"
    );
}
```

---

### Phase 3: Neural Network Loading (1 hour)

**Task**: Load MLP weights from `inst/tandem_net/`

**Files to copy**:
```bash
mkdir -p inst/tandem_net
cp src/tandem/tandem_64/net/*.dat inst/tandem_net/
```

**Modify TANDEM code** to load from package directory:

```cpp
// In pitch.cpp or wherever networks are loaded
void pitchMask::readNet(char *fname1, char *fname2, char *fname3) {
    // Get package directory
    std::string pkg_dir = R_PACKAGE_DIR;  // Set by R
    std::string net_dir = pkg_dir + "/tandem_net/";
    
    std::string file1 = net_dir + "MLP1.64.dat";
    std::string file2 = net_dir + "MLP2.64.dat";
    std::string file3 = net_dir + "MLP3.64.dat";
    
    // Load networks from inst/tandem_net/
    // ... existing loading code ...
}
```

---

### Phase 4: Update Makevars (30 min)

**File**: `src/Makevars`

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

# Add to CXX_SOURCES
CXX_SOURCES += tandem_wrapper.cpp $(TANDEM_SOURCES)
```

---

### Phase 5: Testing & Validation (1-2 hours)

**Test Cases**:

1. **Clean speech**:
```r
test_that("TANDEM tracks pitch in clean speech", {
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  result <- trk_tandem(test_wav, toFile = FALSE, verbose = FALSE)
  
  expect_s3_class(result, "AsspDataObj")
  expect_true(all(!is.na(result$pitch[50:150])))  # Should have pitch
  expect_true(all(result$pitch[50:150] >= 50 & result$pitch[50:150] <= 500))
})
```

2. **Noisy speech**:
```r
test_that("TANDEM handles noisy speech", {
  # Add white noise to clean speech
  # Verify TANDEM still tracks pitch
})
```

3. **Compare with RAPT**:
```r
test_that("TANDEM pitch comparable to RAPT", {
  test_wav <- ...
  tandem_pitch <- trk_tandem(test_wav, toFile = FALSE)$pitch
  rapt_pitch <- trk_rapt(test_wav, toFile = FALSE)$fm
  
  # Calculate correlation
  cor_value <- cor(tandem_pitch, rapt_pitch, use = "complete.obs")
  expect_true(cor_value > 0.8)  # Should be reasonably correlated
})
```

---

## Protoscribe Integration Opportunity

**Note for protoscribe package integration**:

TANDEM's **voiced mask** output could be used in protoscribe for:

### `draft_vad_tandem()` - Voice Activity Detection from TANDEM

**Concept**: Use TANDEM's voiced speech segregation as VAD signal

**Implementation** (in protoscribe, not superassp):

```r
#' Draft VAD from TANDEM Voiced Mask
#'
#' Uses TANDEM algorithm's voiced speech segregation for voice activity detection.
#' Complementary to draft_vad_brouhaha() - TANDEM is better for multi-speaker,
#' Brouhaha is better for simple VAD.
#'
#' @param bundle emuR database bundle
#' @param ... Additional arguments
#' @export
draft_vad_tandem <- function(bundle, ...) {
  # Get audio
  audio_path <- get_audio_path(bundle)
  
  # Run TANDEM with mask output
  result <- superassp::trk_tandem(audio_path, return_mask = TRUE, toFile = FALSE)
  
  # Sum voiced mask across channels to get overall voicing strength
  mask_sum <- colSums(result$voiced_mask)
  
  # Threshold to binary VAD
  vad_binary <- mask_sum > threshold
  
  # Convert to Suggestion object
  segments <- binary_to_segments(vad_binary, sample_rate = 100)  # 100 Hz frame rate
  
  create_suggestion(
    bundle = bundle,
    level = "vad",
    type = "SEGMENT",
    segments = segments,
    method = "TANDEM",
    confidence = calculate_confidence(mask_sum),
    cache_key = generate_cache_key("tandem", audio_path)
  )
}
```

**Benefits over Brouhaha**:
- Multi-speaker capable (can separate multiple voices)
- More detailed time-frequency information
- Physics-based (gammatone filterbank + neural networks)

**When to use**:
- Multi-speaker recordings
- Overlapping speech
- Need for detailed speech/non-speech separation

**When NOT to use**:
- Simple VAD (Brouhaha is faster)
- Need SNR/C50 estimates (Brouhaha provides these)

---

## Estimated Timeline

| Phase | Time | Dependencies |
|-------|------|--------------|
| **Phase 1**: Refactor source | 3-4 hours | None |
| **Phase 2**: Update wrapper | 2 hours | Phase 1 complete |
| **Phase 3**: Neural networks | 1 hour | None (parallel) |
| **Phase 4**: Build system | 30 min | Phase 1-2 complete |
| **Phase 5**: Testing | 1-2 hours | All phases complete |
| **Total** | **7.5-9.5 hours** | Sequential + testing |

**Minimum viable**: 6 hours (skip extensive testing)  
**Full integration**: 8-10 hours (with validation)

---

## Current Workaround

**For users who need TANDEM now**:

The placeholder implementation provides:
1. ✅ Function interface (matches final implementation)
2. ✅ Documentation (complete)
3. ✅ Zero-crossing pitch (crude but functional)
4. ✅ Proper data structure (AsspDataObj)

**For users who need accurate TANDEM**:

Use standalone TANDEM:
```bash
cd src/tandem/tandem_64
make
./tandem input.txt output
```

Then import results into R.

---

## Decision

**Status**: ✅ **DEFERRED** (prudent choice)

**Rationale**:
1. ✅ Framework complete - function exists and works
2. ⚠️ Full integration complex - 6-8 hours careful work needed
3. 🎯 Better to do it right than fast
4. 📋 Clear implementation plan exists
5. 🔧 Workarounds available for urgent needs

**When to proceed**:
- Dedicated time block available (full day)
- Test suite expanded
- Validation strategy in place
- User demand justifies effort

---

**Document Version**: 1.0  
**Date**: 2025-11-07  
**Status**: Deferred pending careful implementation
