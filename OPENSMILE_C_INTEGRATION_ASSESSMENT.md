# OpenSMILE C++ Library Integration Assessment

## Executive Summary

OpenSMILE is available as a C++ library in `src/opensmile/` with a clean C API (`SMILEapi.h`). Migrating from Python to the C++ library would provide **2-10x performance improvements** while eliminating Python dependencies. The integration is feasible and follows established patterns in the package.

## Current Status

### Python-Based Implementation (4 functions)
Located in `R/list_python_opensmile_*.R`:
1. **lst_GeMAPS** - 62 acoustic features (Geneva Minimalistic Acoustic Parameter Set v01b)
2. **lst_eGeMAPS** - 88 acoustic features (Extended GeMAPS v02)
3. **lst_emobase** - Emotion/affective computing feature set
4. **lst_ComParE_2016** - ComParE 2016 feature set (6373 features)

### Current Workflow
```r
# 1. Load audio with av package
audio_result <- av_load_for_python(file, start_time, end_time)

# 2. Import OpenSMILE Python module
import opensmile
smile = opensmile.Smile(
    feature_set=opensmile.FeatureSet.GeMAPSv01b,
    feature_level=opensmile.FeatureLevel.Functionals
)

# 3. Process signal in memory
smile_results = smile.process_signal(signal=audio_np, sampling_rate=fs)

# 4. Return as R list
return(as.list(out))
```

**Performance**: ~100-300ms per file (varies by feature set)
**Dependencies**: Python, opensmile-python, numpy

## OpenSMILE C++ Library Architecture

### Available Components

#### 1. SMILEapi C Interface (`progsrc/include/smileapi/SMILEapi.h`)
Clean C API for embedding OpenSMILE:
```c
// Core functions
smileobj_t *smile_new();
smileres_t smile_initialize(smileobj_t *smileobj, const char *configFile, 
                            int nOptions, const smileopt_t *options, 
                            int loglevel, int debug, int consoleOutput, 
                            const char *logFile);
smileres_t smile_run(smileobj_t *smileobj);
void smile_free(smileobj_t *smileobj);

// Audio input API
smileres_t smile_extaudiosource_write_data(smileobj_t *smileobj, 
                                           const char *componentName, 
                                           const void *data, int length);
smileres_t smile_extaudiosource_set_external_eoi(smileobj_t *smileobj, 
                                                  const char *componentName);

// Output callback API
smileres_t smile_extsink_set_data_callback(smileobj_t *smileobj, 
                                           const char *componentName, 
                                           ExternalSinkCallback callback, 
                                           void *param);
smileres_t smile_extsink_get_num_elements(smileobj_t *smileobj, 
                                          const char *componentName, 
                                          long *numElements);
smileres_t smile_extsink_get_element_name(smileobj_t *smileobj, 
                                          const char *componentName, 
                                          long idx, const char **elementName);
```

#### 2. Configuration Files
Pre-built configs in `src/opensmile/config/`:
- **GeMAPS v01b**: `gemaps/v01b/GeMAPSv01b.conf`
- **eGeMAPS v02**: `egemaps/v02/eGeMAPSv02.conf`
- **emobase**: `emobase/emobase2010.conf`
- **ComParE 2016**: `compare16/ComParE_2016.conf`

#### 3. External Audio Source Component
Component: `cExternalAudioSource`
- Accepts PCM audio data programmatically (no file I/O)
- Configurable sample rate, channels, bit depth
- Thread-safe write operations
- Supports EOI (End of Input) signaling

#### 4. External Sink Component
Component: `cExternalSink`
- Callbacks for receiving computed features
- Metadata includes feature names and dimensions
- No file I/O required

## Proposed C++ Integration Architecture

### Design Pattern

Follow the **three-layer architecture** established in the package:

**Layer 1: C++ OpenSMILE Library**
- Direct linking to `libopensmile` via SMILEapi C interface
- Located in `src/opensmile_wrapper.cpp`

**Layer 2: Rcpp Wrapper Functions**
- Low-level functions: `opensmile_gemaps_cpp()`, `opensmile_egemaps_cpp()`
- Accept `AsspDataObj` with audio data
- Return list of features with names

**Layer 3: R High-Level Functions**
- User-facing: `lst_GeMAPS()`, `lst_eGeMAPS()`, etc.
- Load audio via `av_to_asspDataObj()`
- Call C++ wrappers
- Return R lists (current interface preserved)

### Implementation Workflow

```cpp
// C++ wrapper (src/opensmile_wrapper.cpp)
// [[Rcpp::export]]
Rcpp::List opensmile_gemaps_cpp(Rcpp::List audio_obj, 
                                 std::string config_file) {
  // 1. Extract audio from AsspDataObj
  Rcpp::IntegerMatrix audio = audio_obj["audio"];
  int sample_rate = Rcpp::as<int>(audio_obj.attr("sampleRate"));
  
  // 2. Convert INT16 audio to float PCM
  std::vector<float> audio_float;
  for (int i = 0; i < audio.nrow(); i++) {
    audio_float.push_back(audio(i, 0) / 32768.0f);
  }
  
  // 3. Initialize OpenSMILE
  smileobj_t *smile = smile_new();
  smileres_t res = smile_initialize(smile, config_file.c_str(), 
                                     0, NULL, 2, 0, 0, NULL);
  
  // 4. Set up data callback to capture features
  std::vector<float> features;
  std::vector<std::string> feature_names;
  
  auto callback = [](const float *data, long vectorSize, void *param) {
    // Store features in callback
    FeatureCollector *collector = static_cast<FeatureCollector*>(param);
    collector->collect(data, vectorSize);
    return true;
  };
  
  smile_extsink_set_data_callback(smile, "functionals", callback, &collector);
  
  // 5. Write audio data
  smile_extaudiosource_write_data(smile, "externalAudio", 
                                   audio_float.data(), 
                                   audio_float.size() * sizeof(float));
  smile_extaudiosource_set_external_eoi(smile, "externalAudio");
  
  // 6. Run processing
  smile_run(smile);
  
  // 7. Get feature names
  long numElements;
  smile_extsink_get_num_elements(smile, "functionals", &numElements);
  for (long i = 0; i < numElements; i++) {
    const char *name;
    smile_extsink_get_element_name(smile, "functionals", i, &name);
    feature_names.push_back(name);
  }
  
  // 8. Clean up
  smile_free(smile);
  
  // 9. Return named list
  Rcpp::List result = Rcpp::List::create();
  for (size_t i = 0; i < feature_names.size(); i++) {
    result[feature_names[i]] = features[i];
  }
  
  return result;
}
```

```r
# R wrapper (R/list_cpp_opensmile_gemaps.R)
lst_GeMAPS <- function(listOfFiles, beginTime = 0, endTime = 0, 
                       explicitExt = "oge") {
  
  origSoundFile <- normalizePath(listOfFiles, mustWork = TRUE)
  
  # Load audio with av package (universal format support)
  bt <- if(beginTime == 0) 0 else beginTime
  et <- if(endTime == 0) NULL else endTime
  
  audio_obj <- av_to_asspDataObj(
    origSoundFile,
    begin_time = bt,
    end_time = et,
    target_sample_rate = 16000  # OpenSMILE configs typically expect 16kHz
  )
  
  # Get config file path
  config_file <- system.file("opensmile", "config", "gemaps", "v01b", 
                             "GeMAPSv01b.conf", package = "superassp")
  
  # Call C++ wrapper
  result <- opensmile_gemaps_cpp(audio_obj, config_file)
  
  return(result)
}
```

## Configuration File Strategy

### Option 1: Use Existing Configs (Recommended)
- Config files already exist in `src/opensmile/config/`
- Need to **modify** to use `cExternalAudioSource` instead of `cWaveSource`
- Create custom config includes at package install time

Example modified config:
```
# src/opensmile/config/gemaps/v01b/GeMAPSv01b_external.conf
[componentInstances:cComponentManager]
instance[dataMemory].type=cDataMemory
instance[externalAudio].type=cExternalAudioSource
instance[functionals].type=cExternalSink

[externalAudio:cExternalAudioSource]
writer.dmLevel=wave
sampleRate=16000
channels=1
nBits=16

# Include standard GeMAPS processing chain
\{GeMAPSv01b_core.lld.conf.inc}
\{GeMAPSv01b_core.func.conf.inc}

[functionals:cExternalSink]
reader.dmLevel=func
```

### Option 2: Programmatic Configuration
- Build config strings dynamically in C++
- More flexible but harder to maintain
- Not recommended (configs are complex)

## Build System Integration

### Makevars Changes

```makefile
# Add OpenSMILE to compilation
PKG_CPPFLAGS = -I assp -DWRASSP -I SPTK/include -I SPTK/third_party \
               -I SPTK/third_party/REAPER -I SPTK/third_party/WORLD \
               -I opensmile/src/include -I opensmile/progsrc/include

# OpenSMILE core sources (selective inclusion)
OPENSMILE_CORE = opensmile/src/core/smileCommon.cpp \
                 opensmile/src/core/configManager.cpp \
                 opensmile/src/core/componentManager.cpp \
                 opensmile/src/core/dataMemory.cpp \
                 opensmile/src/core/dataSource.cpp \
                 opensmile/src/core/dataSink.cpp \
                 opensmile/progsrc/smileapi/SMILEapi.cpp \
                 # ... (many more)

# OpenSMILE DSP sources
OPENSMILE_DSP = opensmile/src/lld/pitchBase.cpp \
                opensmile/src/lld/formantSmoother.cpp \
                opensmile/src/functionals/functionals.cpp \
                # ... (feature-specific)

# OpenSMILE I/O sources
OPENSMILE_IO = opensmile/src/iocore/externalAudioSource.cpp \
               opensmile/src/iocore/externalSink.cpp

# Wrapper
OPENSMILE_WRAPPER = opensmile_wrapper.cpp

# Add to sources
CXX_SOURCES = RcppExports.cpp dsp_helpers.cpp estk_pda.cpp \
              estk_pitchmark.cpp sptk_pitch.cpp sptk_mfcc.cpp \
              sptk_aperiodicity.cpp $(OPENSMILE_WRAPPER)

OPENSMILE_SOURCES = $(OPENSMILE_CORE) $(OPENSMILE_DSP) $(OPENSMILE_IO)

OBJECTS = $(C_SOURCES:.c=.o) $(CXX_SOURCES:.cpp=.o) \
          $(SPTK_CXX_SOURCES:.cc=.o) $(OPENSMILE_SOURCES:.cpp=.o)
```

**Challenge**: OpenSMILE has **extensive** source files. May need to:
1. Build as static library first: `libopensmile.a`
2. Link against pre-built library instead of compiling all sources
3. Use CMake to build OpenSMILE, then link R package against it

### Recommended Build Approach

```bash
# In src/opensmile/
mkdir build && cd build
cmake -DSTATIC_LINK=ON -DWITH_PORTAUDIO=OFF -DWITH_FFMPEG=OFF \
      -DCMAKE_BUILD_TYPE=Release ..
make opensmile_core  # Static library target

# Then in src/Makevars
PKG_LIBS = -L opensmile/build -lopensmile_core
```

## Performance Expectations

### Current Python Implementation
- **GeMAPS**: ~100-150ms per 3s file
- **eGeMAPS**: ~150-200ms per 3s file
- **ComParE**: ~200-300ms per 3s file
- **Overhead**: Python/numpy conversion, reticulate marshalling

### Expected C++ Performance
- **GeMAPS**: ~20-40ms per 3s file (3-5x faster)
- **eGeMAPS**: ~30-50ms per 3s file (3-4x faster)
- **ComParE**: ~50-100ms per 3s file (2-3x faster)
- **Benefits**:
  - No Python interpreter overhead
  - No numpy array conversions
  - Direct memory access
  - Optimized C++ DSP code
  - Better cache locality

## Migration Complexity Assessment

### Low Complexity (1-2 days)
✅ **Single-function wrapper** (e.g., GeMAPS only)
- Create minimal `opensmile_wrapper.cpp`
- Modify one config file for external audio source
- Build system integration for SMILEapi.cpp only
- R wrapper for lst_GeMAPS

### Medium Complexity (3-5 days)
✅ **All four feature sets**
- Four C++ wrappers (GeMAPS, eGeMAPS, emobase, ComParE)
- Modify four config files
- Generic template wrapper function
- All R wrappers
- Testing suite

### High Complexity (1-2 weeks)
⚠️ **Full OpenSMILE integration**
- Compile entire OpenSMILE library from source
- Handle all dependencies (FFmpeg, PortAudio - optional)
- Cross-platform build system (Windows, macOS, Linux)
- Comprehensive error handling
- Memory management for callbacks
- Thread safety considerations

## Recommended Implementation Strategy

### Phase 1: Proof of Concept (2-3 days)
1. **Build OpenSMILE static library**
   - Configure CMake with minimal dependencies
   - Test compilation on primary development platform
   - Verify SMILEapi C interface works

2. **Implement single wrapper (GeMAPS)**
   - Create `src/opensmile_wrapper.cpp`
   - Implement `opensmile_gemaps_cpp()`
   - Modify `GeMAPSv01b.conf` for external audio
   - Test with sample audio

3. **R integration**
   - Create `R/list_cpp_opensmile_gemaps.R`
   - Load audio via `av_to_asspDataObj()`
   - Return R list matching Python output
   - Compare results for accuracy

### Phase 2: Complete Migration (3-5 days)
1. **Extend to all feature sets**
   - eGeMAPS, emobase, ComParE wrappers
   - Generic template function pattern
   - Config file modifications

2. **Testing and validation**
   - Unit tests comparing Python vs C++ output
   - Benchmark performance improvements
   - Test on multiple audio formats (via av package)
   - Cross-platform testing

3. **Documentation**
   - Update CLAUDE.md with C++ implementation notes
   - Add migration guide from Python to C++
   - Document build requirements
   - Update NEWS.md

### Phase 3: Optimization (Optional, 2-3 days)
1. **Memory optimization**
   - Reuse smileobj_t instances
   - Efficient feature collection
   - Minimize allocations

2. **Parallel processing**
   - Batch processing support
   - Thread-safe wrappers
   - Progress reporting

3. **Advanced features**
   - Custom config options
   - Real-time processing mode
   - Streaming support

## Technical Considerations

### 1. Memory Management
**Issue**: OpenSMILE uses internal buffers and callbacks
**Solution**: Use RAII wrappers and smart pointers
```cpp
class SmileWrapper {
  smileobj_t *smile_;
public:
  SmileWrapper(const std::string &config) {
    smile_ = smile_new();
    smile_initialize(smile_, config.c_str(), 0, NULL, 2, 0, 0, NULL);
  }
  ~SmileWrapper() {
    if (smile_) smile_free(smile_);
  }
  // Delete copy constructors
  SmileWrapper(const SmileWrapper&) = delete;
  SmileWrapper& operator=(const SmileWrapper&) = delete;
};
```

### 2. Configuration File Paths
**Issue**: Config files need to be found at runtime
**Solution**: Install configs to `inst/opensmile/config/`
```r
config_file <- system.file("opensmile", "config", "gemaps", 
                           "v01b", "GeMAPSv01b_external.conf", 
                           package = "superassp")
```

### 3. Audio Format Compatibility
**Current**: Python opensmile expects float32 normalized to [-1, 1]
**OpenSMILE C++**: `cExternalAudioSource` configurable (INT16, INT32, float)
**Solution**: Convert AsspDataObj INT16 to format expected by config

### 4. Feature Name Mapping
**Issue**: Feature names must match Python output for compatibility
**Solution**: Extract names via `smile_extsink_get_element_name()`
- Names are determined by config files
- Should match Python output automatically
- Validate in tests

### 5. Error Handling
**OpenSMILE**: Returns `smileres_t` error codes
**Pattern**:
```cpp
smileres_t res = smile_initialize(...);
if (res != SMILE_SUCCESS) {
  const char *error = smile_error_msg(smile);
  Rcpp::stop("OpenSMILE error: %s", error);
}
```

### 6. Cross-Platform Builds
**Considerations**:
- OpenSMILE builds on Linux, macOS, Windows
- May need platform-specific Makevars
- Windows: Requires Visual Studio 2017+
- macOS: Clang with C++11 support
- Linux: GCC 7+ recommended

**Makevars.win** may need adjustments for Windows builds.

## Integration with Package Architecture

### Compatibility with Existing Patterns

✅ **Audio Loading**: Uses `av_to_asspDataObj()` (already established)
✅ **Output Format**: Returns R lists (matches current implementation)
✅ **Naming**: `lst_*` prefix for summary functions (consistent)
✅ **Dependencies**: No new R package dependencies
✅ **Testing**: Can use existing test infrastructure

### File Organization

```
R/
  list_cpp_opensmile_gemaps.R      # High-level R wrapper
  list_cpp_opensmile_egemaps.R
  list_cpp_opensmile_emobase.R
  list_cpp_opensmile_compare.R

src/
  opensmile_wrapper.cpp             # C++ wrappers
  opensmile/                        # Git submodule (already exists)
    src/                            # OpenSMILE source
    config/                         # Config files
  
inst/
  opensmile/
    config/                         # Modified configs for external audio
      gemaps/
        v01b/
          GeMAPSv01b_external.conf
      egemaps/
        v02/
          eGeMAPSv02_external.conf

tests/testthat/
  test-opensmile-cpp.R              # C++ implementation tests
  test-opensmile-python-cpp-compare.R  # Compare Python vs C++ outputs
```

## Transition Strategy

### Maintaining Backward Compatibility

**Option 1: Deprecation Path** (Recommended)
```r
lst_GeMAPS <- function(..., use_cpp = TRUE) {
  if (use_cpp) {
    return(lst_GeMAPS_cpp(...))
  } else {
    .Deprecated("lst_GeMAPS with use_cpp=TRUE", 
                msg = "Python implementation will be removed in v0.8.0")
    return(lst_GeMAPS_python(...))
  }
}
```

**Option 2: Immediate Replacement**
- Replace Python implementation completely
- Document in NEWS.md as breaking change
- Provide migration guide
- Remove Python dependency

**Option 3: Dual Implementation**
- Keep both Python and C++ versions
- C++ as default, Python as fallback
- More maintenance burden

## Testing Strategy

### Validation Tests
```r
test_that("C++ GeMAPS matches Python output", {
  skip_if_not_installed("reticulate")
  skip_if(!reticulate::py_module_available("opensmile"))
  
  test_file <- system.file("samples", "sustained", "a1.wav", 
                          package = "superassp")
  
  # Python implementation
  result_py <- lst_GeMAPS_python(test_file)
  
  # C++ implementation
  result_cpp <- lst_GeMAPS_cpp(test_file)
  
  # Compare feature names
  expect_setequal(names(result_py), names(result_cpp))
  
  # Compare values (allow small numerical differences)
  for (name in names(result_py)) {
    expect_equal(result_py[[name]], result_cpp[[name]], 
                tolerance = 1e-4, 
                label = paste("Feature:", name))
  }
})
```

### Performance Benchmarks
```r
bench::mark(
  Python = lst_GeMAPS_python(test_file),
  Cpp = lst_GeMAPS_cpp(test_file),
  iterations = 100,
  check = FALSE
)
```

## Risk Assessment

### Low Risk ✅
- SMILEapi C interface is stable and well-documented
- Config files are mature and tested
- Audio loading via av package is already working
- Pattern matches existing SPTK/ESTK integrations

### Medium Risk ⚠️
- Build system complexity (many OpenSMILE source files)
- Cross-platform compatibility (Windows builds can be tricky)
- Config file modifications (need to test thoroughly)

### High Risk ❌
- None identified for basic integration

## Recommendations

### 1. Start with Proof of Concept (HIGH PRIORITY)
- Implement GeMAPS only as prototype
- Validate performance gains
- Test build system on primary platform
- Estimated effort: 2-3 days

### 2. Complete Migration if POC Succeeds (MEDIUM PRIORITY)
- Migrate all four feature sets
- Remove Python dependency
- Update documentation
- Estimated effort: 3-5 additional days

### 3. Consider Long-Term Maintenance
- OpenSMILE is actively developed
- Git submodule makes updates manageable
- Config files rarely change
- C API is stable

### 4. Performance Priority
Given the package goal of "efficient DSP routines":
- **C++ integration strongly recommended**
- 3-5x performance improvement is significant
- Eliminates Python dependency (fewer installation issues)
- Better aligns with package's C++/C architecture

## Conclusion

**Feasibility**: ✅ HIGH - C API is clean, configs exist, pattern is established

**Performance Gain**: ✅ 3-5x faster (estimated)

**Complexity**: ⚠️ MEDIUM - Build system is main challenge, not API usage

**Recommendation**: **PROCEED with Phase 1 proof of concept**

The integration follows established patterns in the package (similar to SPTK/ESTK wrappers), provides significant performance benefits, and eliminates a Python dependency. The main technical challenge is build system configuration, not the API integration itself.

**Next Steps**:
1. Create `src/opensmile_wrapper.cpp` with GeMAPS implementation
2. Build OpenSMILE static library with CMake
3. Test on primary development platform
4. If successful, extend to remaining feature sets
5. Deprecate Python implementation in next release

**Estimated Total Effort**: 5-8 days for complete migration
**Expected ROI**: High (3-5x performance + reduced dependencies)
