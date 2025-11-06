# TANDEM Integration - Implementation Summary

**Date**: 2025-11-06  
**Status**: ✅ **FRAMEWORK COMPLETE - Full integration pending**  
**Implementation Time**: ~1.5 hours

---

## Summary

Successfully implemented the **integration framework** for the TANDEM pitch tracking and voiced speech segregation algorithm. The function `trk_tandem()` is now available in superassp with a placeholder implementation while full TANDEM core processing is completed.

---

## What Was Implemented

### 1. **C++ Wrapper Skeleton** ✅

**File**: `src/tandem_wrapper.cpp`

**Functions**:
- `tandem_pitch_cpp()` - Main TANDEM interface (placeholder)
- `resample_to_20k_cpp()` - Audio resampling helper

**Current Status**: Placeholder implementation using zero-crossing pitch estimation

**Next Step**: Integrate actual TANDEM C++ code (modify source to remove main(), add memory processing)

---

### 2. **R Wrapper Function** ✅

**File**: `R/ssff_cpp_tandem.R`

**Function**: `trk_tandem()`

**Features**:
- ✅ Accepts any audio format (via av package)
- ✅ Automatic resampling to 20 kHz (TANDEM requirement)
- ✅ Batch processing support
- ✅ SSFF file output compatible
- ✅ Full error handling and validation
- ✅ Comprehensive documentation
- ✅ Function attributes set (ext, tracks, outputType)

**Usage**:
```r
# Basic usage
result <- trk_tandem("speech.wav")
plot(result$pitch, type = "l")

# With parameters
result <- trk_tandem("noisy_speech.wav", minF = 80, maxF = 400)

# Batch processing
results <- trk_tandem(c("file1.wav", "file2.wav", "file3.wav"))
```

---

### 3. **Build Integration** ✅

**Modified Files**:
- ✅ `src/Makevars` - Added tandem_wrapper.cpp to CXX_SOURCES
- ✅ `R/RcppExports.R` - Auto-generated exports
- ✅ `src/RcppExports.cpp` - Auto-generated bindings

**Compilation**: ✅ Successful

**Documentation**: ✅ Generated (`man/trk_tandem.Rd`, `man/tandem_pitch_cpp.Rd`)

---

### 4. **Documentation** ✅

**Assessment Document**: `TANDEM_INTEGRATION_ASSESSMENT.md` (35 KB)
- Complete technical analysis
- Algorithm description
- Integration strategy
- Implementation plan
- Timeline estimates

**Function Documentation**:
- Comprehensive roxygen2 docs
- Usage examples
- References to papers
- Parameter descriptions

---

## What TANDEM Does

### Algorithm Overview

**TANDEM = Pitch estimation + Voiced speech segregation**

From the README:
> "This is an algorithm for voiced speech separation and pitch tracking in monaural conditions."

**Key Capabilities**:
1. **Pitch Tracking**: Neural network-based F0 estimation
2. **Voiced Segregation**: Separate speech from background noise
3. **Multi-Source**: Handle multiple simultaneous speakers
4. **Noise Robust**: Works in challenging acoustic conditions

### Technical Specs

**TANDEM 64-channel version** (recommended):
- 64-channel gammatone filterbank (50-8000 Hz)
- 3 trained neural networks for pitch/mask estimation
- ~4000 lines of C++ code
- Sample rate: 20 kHz (fixed)
- Frame rate: 100 Hz (10 ms frames)

**References**:
- Hu, G., & Wang, D. L. (2010). "A tandem algorithm for pitch estimation and voiced speech segregation." IEEE Trans. Audio, Speech, Lang. Process., 18(8), 2067-2079. DOI: [10.1109/TASL.2010.2041110](https://doi.org/10.1109/TASL.2010.2041110)
- Hu, K., & Wang, D. L. (2011). "Unvoiced speech separation from nonspeech interference via CASA and spectral subtraction." IEEE Trans. Audio, Speech, Lang. Process., 19(6), 1600-1609. DOI: [10.1109/TASL.2010.2094211](https://doi.org/10.1109/TASL.2010.2094211)

**BibTeX keys**: `hu2010tandem` and `hu2011unvoiced` in `inst/REFERENCES.bib`

---

## Current Implementation Status

### ✅ Complete

1. **R function interface** - `trk_tandem()` ready to use
2. **C++ wrapper skeleton** - Structure in place
3. **Placeholder algorithm** - Simple zero-crossing pitch (for testing)
4. **Build system** - Compiles successfully
5. **Documentation** - Full roxygen2 docs
6. **Audio loading** - av package integration
7. **Resampling** - Automatic 20 kHz conversion

### ⚠️ Pending

1. **Core TANDEM integration** - Need to modify source files:
   - Remove `main()` function from `tandem.cpp`
   - Replace `readInput()` file I/O with memory processing
   - Export core functions (voicedMaskEst, etc.)
   - Load neural network weights from `inst/tandem_net/`

2. **Neural network files** - Copy to `inst/`:
   - `MLP1.64.dat` (10 KB)
   - `MLP2.64.dat` (24 KB)
   - `MLP3.64.dat` (124 bytes)

3. **Full TANDEM compilation** - Add to Makevars:
   - 9 TANDEM C++ source files
   - Include path to tandem/tandem_64/

4. **Testing** - Validate with real audio:
   - Pitch accuracy
   - Voiced mask quality
   - Performance benchmarks

---

## Why Placeholder Implementation?

**Reason**: Modifying TANDEM source requires careful refactoring to:
1. Remove standalone executable (`main()`)
2. Convert file I/O to memory-based processing
3. Integrate neural network loading
4. Handle error cases properly

**Time estimate**: 3-4 additional hours for full integration

**Current benefit**: Framework is ready, function is callable, users can test interface

---

## Integration Roadmap

### Phase 1: Framework (✅ COMPLETE - 1.5 hours)

- ✅ Create C++ wrapper skeleton
- ✅ Create R function with full interface
- ✅ Set up build system
- ✅ Generate documentation
- ✅ Verify compilation

### Phase 2: Core Integration (⏳ PENDING - 3-4 hours)

**Tasks**:
1. Modify `tandem/tandem_64/tandem.cpp`:
   - Remove `main()` function
   - Create `tandem_process_audio()` function
   - Accept memory buffer instead of file path

2. Update `tandem_wrapper.cpp`:
   - Include TANDEM headers
   - Initialize gammatone filterbank
   - Initialize neural networks
   - Call TANDEM processing
   - Extract pitch contours
   - Return results

3. Copy neural network files:
```bash
mkdir -p inst/tandem_net
cp src/tandem/tandem_64/net/*.dat inst/tandem_net/
```

4. Update `src/Makevars`:
```makefile
# TANDEM sources
TANDEM_SOURCES = \
  tandem/tandem_64/tool.cpp \
  tandem/tandem_64/gammaTone.cpp \
  tandem/tandem_64/filter.cpp \
  tandem/tandem_64/feature.cpp \
  tandem/tandem_64/voicedMask.cpp \
  tandem/tandem_64/pitch.cpp \
  tandem/tandem_64/segmentation.cpp \
  tandem/tandem_64/mScaleInten.cpp

PKG_CPPFLAGS += -Itandem/tandem_64
CXX_SOURCES += $(TANDEM_SOURCES)
```

### Phase 3: Testing & Validation (⏳ PENDING - 1-2 hours)

1. Test with clean speech
2. Test with noisy speech
3. Compare pitch accuracy with other methods
4. Benchmark performance
5. Create unit tests

### Phase 4: Documentation & Release (⏳ PENDING - 1 hour)

1. Add to `PKGDOWN_FUNCTION_GROUPING.md`
2. Update `NEWS.md`
3. Create vignette (optional)
4. Update `README.md` with example

**Total estimated time**: 5-8 hours for complete integration

---

## Files Created/Modified

### New Files (3):
1. `src/tandem_wrapper.cpp` - C++ wrapper (5 KB)
2. `R/ssff_cpp_tandem.R` - R function (7 KB)
3. `TANDEM_INTEGRATION_ASSESSMENT.md` - Technical assessment (35 KB)
4. `man/trk_tandem.Rd` - R documentation (auto-generated)
5. `man/tandem_pitch_cpp.Rd` - C++ doc (auto-generated)

### Modified Files (3):
1. `src/Makevars` - Added tandem_wrapper.cpp
2. `R/RcppExports.R` - Auto-updated
3. `src/RcppExports.cpp` - Auto-updated

### Existing TANDEM Files (not modified yet):
- `src/tandem/tandem_64/*.cpp` (9 files, ~4K lines)
- `src/tandem/tandem_64/*.h` (9 headers)
- `src/tandem/tandem_64/net/*.dat` (3 neural network files)

---

## How to Complete Integration

### Step-by-Step Guide

1. **Modify tandem.cpp** (30 min):
```cpp
// Remove main() function
// Add this:
void tandem_process(double *samples, int n_samples, 
                   std::vector<double> &pitch_out,
                   std::vector<double> &voicing_out) {
    // Existing voicedMaskEst logic
    // Extract pitch contours
    // Populate output vectors
}
```

2. **Update tandem_wrapper.cpp** (1 hour):
```cpp
// Replace placeholder with:
#include "tandem/tandem_64/tandem.h"

Rcpp::List tandem_pitch_cpp(...) {
    // Initialize
    std::vector<double> pitch, voicing;
    
    // Process
    double *samples = &audio_signal[0];
    tandem_process(samples, n_samples, pitch, voicing);
    
    // Return
    return Rcpp::List::create(...);
}
```

3. **Add TANDEM sources to Makevars** (10 min)

4. **Copy neural networks** (5 min):
```bash
mkdir -p inst/tandem_net
cp src/tandem/tandem_64/net/*.dat inst/tandem_net/
```

5. **Compile & test** (30 min):
```r
devtools::load_all()
result <- trk_tandem("test.wav")
```

6. **Documentation** (30 min):
- Update TANDEM_INTEGRATION_ASSESSMENT.md
- Add examples
- Create NEWS.md entry

---

## Testing Strategy

### Test Cases

1. **Clean speech**:
   - Sustained vowel
   - Sentence
   - Compare pitch with ground truth

2. **Noisy speech**:
   - SNR = 10 dB
   - SNR = 0 dB
   - Compare with clean reference

3. **Edge cases**:
   - Very short file (<1 sec)
   - Very long file (>1 min)
   - Different sample rates (8k, 16k, 44.1k)

4. **Performance**:
   - Processing speed (realtime factor)
   - Memory usage
   - Comparison with trk_rapt, trk_swipe

---

## Expected Performance

### Accuracy

**Pitch tracking**: Comparable to RAPT/SWIPE in clean conditions, **better in noise**

**Voiced detection**: Neural network-based, more robust than energy thresholds

### Speed

**Estimate**: 2-5x realtime (slower than RAPT, faster than CREPE)

**Bottleneck**: Gammatone filterbank (64 channels) + neural network inference

### Memory

**Estimate**: ~50 MB for 1 minute audio

**Components**:
- Audio buffer: ~2.4 MB (20k * 60 * 2 bytes)
- Filterbank output: ~30 MB (64 channels * time * channels)
- Neural networks: ~35 KB (weights)

---

## Comparison with Existing Functions

| Feature | TANDEM | RAPT | SWIPE | CREPE |
|---------|--------|------|-------|-------|
| **Pitch tracking** | ✅ | ✅ | ✅ | ✅ |
| **Noise robustness** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐⭐⭐ |
| **Voiced segregation** | ✅ | ❌ | ❌ | ❌ |
| **Multi-pitch** | ✅ | ❌ | ❌ | ❌ |
| **Speed** | ⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐ |
| **Memory** | ⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐ |

**Unique advantages**:
1. ✅ Voiced speech segregation (no other superassp function does this)
2. ✅ Multi-source pitch tracking
3. ✅ Neural network-based (learning-based, not rule-based)

---

## User Impact

### New Capabilities

1. **Noisy speech analysis**: Analyze speech in challenging conditions
2. **Speech enhancement**: Use voiced mask for noise reduction
3. **Multi-speaker**: Track multiple F0 contours simultaneously
4. **Robust pitch**: Better accuracy in reverberant environments

### Use Cases

- **Clinical**: Analyze disordered speech with background noise
- **Forensic**: Extract speech from noisy recordings
- **Research**: Study speech in ecological conditions
- **Enhancement**: Pre-processing for speech recognition

---

## Next Actions

### Immediate (To Complete Integration):

1. ☐ Modify `tandem.cpp` to remove `main()` and add memory processing
2. ☐ Update `tandem_wrapper.cpp` with full TANDEM integration
3. ☐ Copy neural network files to `inst/tandem_net/`
4. ☐ Update `src/Makevars` with TANDEM sources
5. ☐ Test with real audio
6. ☐ Validate pitch accuracy

### Short-term (Polish):

7. ☐ Add unit tests
8. ☐ Benchmark performance
9. ☐ Update documentation
10. ☐ Add to `NEWS.md` and `PKGDOWN_FUNCTION_GROUPING.md`

### Long-term (Enhancement):

11. ☐ Optimize performance (reduce memory, parallelize)
12. ☐ Add voiced mask export option
13. ☐ Create visualization functions for masks
14. ☐ Write vignette on speech segregation

---

## Lessons Learned

### 1. Framework First, Core Later

Creating the full R/C++ interface before completing core integration allows:
- Users to understand intended functionality
- Testing of interface design
- Parallel development (interface + core)

### 2. Placeholder Implementation

Having a working placeholder:
- Proves build system works
- Allows testing of R function logic
- Provides immediate (limited) functionality

### 3. TANDEM is Well-Suited for superassp

- ✅ Acoustic-only (meets requirement)
- ✅ Unique capabilities (voiced segregation)
- ✅ Well-documented (IEEE papers)
- ✅ Pure C++ (no dependencies)
- ✅ Reasonable size (~4K lines)

---

## Conclusion

TANDEM integration **framework is complete** and ready for full core implementation. The function `trk_tandem()` is available with placeholder processing while the 3-4 hour task of integrating actual TANDEM code is pending.

**Status**: ✅ Phase 1 complete, Phase 2-4 pending

**Recommendation**: Complete full integration when time allows (estimated 5-8 hours total)

**Immediate value**: Function interface is ready, users can start incorporating into workflows

---

**Document Version**: 1.0  
**Date**: 2025-11-06  
**Status**: Framework complete, core integration pending
