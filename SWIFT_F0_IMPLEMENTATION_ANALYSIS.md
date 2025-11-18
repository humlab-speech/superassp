# Swift-F0 Implementation Analysis for trk_swiftf0()

**Date:** October 20, 2025
**Task:** Implement `trk_swiftf0()` function for Swift-F0 pitch tracker
**Model Location:** `inst/onnx/swift-f0/`

---

## Overview

Swift-F0 is a fast and accurate F0 detector using a CNN-based approach on STFT spectrograms. According to benchmarks, it outperforms CREPE in both speed and accuracy.

**Key Specifications:**
- **Speed:** 132ms for 5 seconds of audio (CPU)
- **Frequency Range:** 46.875 Hz to 2093.75 Hz (G1 to C7)
- **Sample Rate:** 16kHz (required)
- **Hop Length:** 256 samples (16ms frames)
- **Frame Length:** 1024 samples (STFT window)
- **Model Format:** ONNX (model.onnx, ~800KB)

---

## Three Implementation Approaches

### Approach 1: Python via reticulate ✅ RECOMMENDED

**Implementation:**
- Use reticulate to call Swift-F0 Python package
- Load audio via av package → pass to Python
- Return PitchResult → convert to AsspDataObj

**Advantages:**
1. ✅ **Simplest implementation** (~100 lines of R code)
2. ✅ **Zero compilation** - works immediately
3. ✅ **Automatic updates** - pip install swift-f0 updates model
4. ✅ **Proven reliability** - uses official implementation
5. ✅ **Fast development** - reticulate integration already established
6. ✅ **Dependencies met** - Python 3.12 + onnxruntime + numpy already available
7. ✅ **Resampling handled** - librosa in Python does 16kHz conversion
8. ✅ **Cross-platform** - Works on macOS, Linux, Windows

**Disadvantages:**
1. ⚠️ Requires Python + onnxruntime dependency
2. ⚠️ Slightly slower than native C++ (~10-20% overhead from reticulate)
3. ⚠️ Installation: users must run `install_voice_analysis()` or similar

**Performance Estimate:**
- 2-second audio: ~50-80ms (model: 50ms, reticulate overhead: 5-30ms)
- Acceptable for typical speech analysis workflows

**Code Pattern:**
```r
trk_swiftf0 <- function(listOfFiles, minF = 75, maxF = 400,
                        confidence_threshold = 0.9, toFile = TRUE, ...) {
  # Load audio via av
  audio <- av_to_asspDataObj(listOfFiles)

  # Get Python module
  swift_f0 <- reticulate::import("swift_f0")

  # Create detector
  detector <- swift_f0$SwiftF0(fmin = minF, fmax = maxF,
                               confidence_threshold = confidence_threshold)

  # Run detection
  result <- detector$detect_from_array(audio_samples, sample_rate)

  # Convert to AsspDataObj
  obj <- list(
    F0 = result$pitch_hz,
    confidence = result$confidence,
    voicing = result$voicing
  )
  attr(obj, "sampleRate") <- sample_rate
  attr(obj, "startTime") <- result$timestamps[1]
  # ... etc
  class(obj) <- "AsspDataObj"

  return(obj)
}
```

---

### Approach 2: onnx2c Compiler (C Code Generation)

**Implementation:**
- Compile `model.onnx` to C code using `src/onnx2c/`
- Integrate generated C into R via Rcpp
- Implement STFT preprocessing in C++
- Return AsspDataObj from C++

**Advantages:**
1. ✅ **Maximum speed** - Pure C/C++ with no Python overhead
2. ✅ **No Python dependency** - Self-contained
3. ✅ **Small footprint** - Compiled code is compact
4. ✅ **Predictable performance** - No interpreter overhead

**Disadvantages:**
1. ❌ **Complex implementation** - Need to implement STFT preprocessing in C++
2. ❌ **High development time** - Estimated 2-3 days for full implementation
3. ❌ **Model updates difficult** - Requires recompilation for new models
4. ❌ **Limited ONNX operator support** - onnx2c only supports 91/166 operators
5. ❌ **Testing overhead** - Need to verify STFT matches Python exactly
6. ❌ **Maintenance burden** - Must maintain C++ STFT + model wrapper
7. ⚠️ **Risk of incompatibility** - Model might use unsupported operators
8. ⚠️ **Build complexity** - Requires protobuf, cmake during package build

**Performance Estimate:**
- 2-second audio: ~30-50ms (best case)
- Speedup vs Python: ~1.5-2x faster
- But: development time >> performance gain

**Verdict:** Not worth the complexity for this use case

---

### Approach 3: R torch Package

**Implementation:**
- Use torch R package to load ONNX model
- Implement STFT in R/torch
- Run model inference
- Return AsspDataObj

**Advantages:**
1. ✅ **Pure R solution** - No Python dependency
2. ✅ **torch available** - Package already installed
3. ✅ **GPU support** - Could leverage CUDA if available

**Disadvantages:**
1. ❌ **torch ONNX support limited** - May not load all ONNX models correctly
2. ❌ **STFT implementation needed** - Must reimplement preprocessing in R
3. ❌ **Higher memory usage** - torch tensors less efficient than numpy
4. ❌ **Complex debugging** - torch errors harder to diagnose
5. ❌ **Development time** - Estimated 1-2 days
6. ⚠️ **Resampling needed** - Must handle 16kHz conversion in R
7. ⚠️ **Model compatibility uncertain** - Need to test if torch can load model

**Performance Estimate:**
- 2-second audio: ~60-100ms (model: 50ms, R overhead: 10-50ms)
- Similar to Python reticulate approach but more complex

**Verdict:** More complex than Python with no clear benefit

---

## Recommendation: Approach 1 (Python via reticulate)

**Rationale:**

1. **Fastest Development:** ~1-2 hours vs 1-3 days for other approaches
2. **Proven Reliability:** Uses official Swift-F0 implementation
3. **Easy Maintenance:** pip install swift-f0 handles updates
4. **Acceptable Performance:** 50-80ms for 2s audio is fine for speech analysis
5. **Consistency:** Follows existing pattern from voice_sauce, COVAREP integrations
6. **Cross-platform:** Works everywhere Python + onnxruntime work

**Performance is adequate:**
- RAPT (C++): ~100-200ms per 3s audio
- Swift-F0 (Python): ~75-120ms per 3s audio
- Speedup from onnx2c: ~25-40ms (not significant)

**Development cost-benefit:**
- Python approach: 1-2 hours → production ready
- onnx2c approach: 2-3 days → marginal speedup
- torch approach: 1-2 days → similar performance to Python

---

## Implementation Plan (Python/reticulate)

### Phase 1: Core Implementation (1 hour)

1. Create `R/ssff_python_swiftf0.R` following existing patterns:
   - Use `processMediaFiles_LoadAndProcess()` pattern
   - Load audio via `av_to_asspDataObj()`
   - Pass to Python swift_f0 module
   - Convert PitchResult to AsspDataObj
   - Support time windowing (beginTime, endTime)

2. Function signature:
```r
trk_swiftf0(listOfFiles,
           minF = 75,           # Speech default (was 46.875 for music)
           maxF = 400,          # Speech default (was 2093.75 for music)
           confidence_threshold = 0.9,
           beginTime = 0.0,
           endTime = 0.0,
           toFile = TRUE,
           explicitExt = "f0",
           outputDirectory = NULL,
           verbose = TRUE)
```

3. Output tracks:
   - "F0" (Hz)
   - "confidence" (0-1)
   - "voicing" (boolean)

### Phase 2: Installation Helper (30 min)

Add to existing `install_voice_analysis()` or create new helper:
```r
install_swiftf0 <- function(method = "auto", envname = NULL) {
  # Install swift-f0 via pip
  reticulate::py_install("swift-f0", envname = envname, pip = TRUE)
}

swiftf0_available <- function() {
  reticulate::py_module_available("swift_f0")
}
```

### Phase 3: Testing (1 hour)

Create `tests/testthat/test-swiftf0.R`:
- Basic functionality with sustained vowel
- Custom F0 range (minF, maxF)
- Confidence threshold parameter
- Time windowing
- Batch processing
- File I/O modes (toFile=TRUE/FALSE)
- Non-WAV formats
- S7 dispatch (AVAudio input)

### Phase 4: Documentation (30 min)

- Roxygen2 documentation
- Usage examples
- Comparison with other pitch trackers
- Performance characteristics

---

## Expected Performance Characteristics

**Benchmark (3-second sustained vowel @44.1kHz):**
- Audio loading + resampling: ~10-20ms
- Swift-F0 inference: ~75-100ms
- Result conversion: ~5-10ms
- **Total: ~90-130ms**

**Comparison with existing trackers:**
- RAPT (C++): ~100-200ms
- REAPER (C++): ~80-150ms
- DIO (C++): ~50-100ms
- SWIPE (C++): ~150-250ms
- Swift-F0 (Python): ~90-130ms ← **Competitive!**

---

## Dependencies

**Required:**
- Python 3.7+
- onnxruntime (pip)
- numpy (pip)
- swift-f0 (pip)

**Optional:**
- librosa (for resampling, but we handle via av package)

**Installation:**
```bash
pip install swift-f0
# This installs onnxruntime and numpy as dependencies
```

---

## Integration with S7 Dispatch

Function will automatically support S7 dispatch like all other trk_* functions:

```r
# Pattern 1: File path
result <- trk_swiftf0("speech.wav", toFile = FALSE)

# Pattern 2: AVAudio object (automatic S7 dispatch)
audio <- read_avaudio("speech.wav", sample_rate = 16000)
result <- trk_swiftf0(audio, toFile = FALSE)
```

---

## Conclusion

**Recommended Approach:** Python via reticulate (Approach 1)

**Implementation Time:** 2-3 hours total
**Performance:** Competitive with C++ implementations
**Maintenance:** Low (pip install handles updates)
**Risk:** Low (proven implementation)

**Next Steps:**
1. Implement `R/ssff_python_swiftf0.R`
2. Add installation helper
3. Create comprehensive tests
4. Add documentation
5. Commit with proper attribution

---

**Analysis Date:** October 20, 2025
**Status:** Analysis complete, ready for implementation
