# Strategy: Direct Linking to Praat C Library

## Overview

Instead of reimplementing Praat algorithms in C++ or calling through Python, we can link directly to the Praat C library that's already bundled in the Parselmouth submodule at `src/Parselmouth/praat/`.

## Current Architecture

### Parselmouth Approach
```
Python → pybind11 → Praat C library (compiled as static lib)
```

### Our New Approach
```
R → Rcpp → Praat C library (same compiled code as Parselmouth uses)
```

## Advantages

1. **Best Performance**: Use the exact same optimized Praat C code as Parselmouth
2. **No Python Dependency**: Direct R/Rcpp integration
3. **Maintained Code**: Leverage decades of Praat development
4. **Multiple Functions**: Easy to add more Praat procedures beyond intensity

## Implementation Strategy

### Phase 1: Build Praat Library for R Package

1. **Identify Required Praat Sources**
   - Core: `melder/`, `sys/`, `stat/`
   - DSP: `fon/` (includes Sound, Intensity, Pitch, Formant, etc.)
   - Optional: `dwtools/` for extended analyses

2. **Create Praat Library Makefile**
   - Extract relevant sources from Parselmouth's CMakeLists
   - Compile as static library: `libpraat.a`
   - Include in R package Makevars

3. **Header Management**
   - Expose necessary Praat headers
   - Handle C++ wrappers around C code where needed

### Phase 2: Create Rcpp Wrappers

Pattern for each Praat function:

```cpp
// In src/praat_wrappers.cpp

#include <Rcpp.h>

// Praat headers
extern "C" {
#include "Parselmouth/praat/sys/melder.h"
}
#include "Parselmouth/praat/fon/Sound.h"
#include "Parselmouth/praat/fon/Intensity.h"
#include "Paraat/fon/Sound_to_Intensity.h"

//' Praat Intensity (Direct Library)
// [[Rcpp::export]]
List praat_intensity_direct(SEXP audio_obj, 
                           double minimal_f0_frequency = 50.0,
                           double time_step = 0.0, 
                           bool subtract_mean = true) {
  // Convert AsspDataObj to Praat Sound
  autoSound sound = asspDataObj_to_Sound(audio_obj);
  
  // Call Praat function directly
  autoIntensity intensity = Sound_to_Intensity(sound.get(), 
                                               minimal_f0_frequency, 
                                               time_step, 
                                               subtract_mean);
  
  // Convert back to R format
  return intensity_to_list(intensity.get());
}
```

### Phase 3: Handle Praat's Memory Management

Praat uses its own memory management (`Thing`, `autoThing`):
- Use Praat's `auto*` types for automatic cleanup
- Initialize Melder on first use
- Handle exceptions properly

### Phase 4: Systematic Migration

Migrate existing Parselmouth-based functions:
1. `praat_intensity` → `praat_intensity_direct`
2. `praat_formant_burg` → `praat_formant_direct`
3. `praat_pitch` → `praat_pitch_direct`
4. Other analyses as needed

## Technical Considerations

### Include Paths

Praat expects specific include structure:
```makefile
PKG_CPPFLAGS = -I Parselmouth/praat/melder \
               -I Parselmouth/praat/sys \
               -I Parselmouth/praat/stat \
               -I Parselmouth/praat/fon \
               -I Parselmouth/praat/dwtools
```

### Compilation Flags

From Parselmouth's CMakeLists:
```makefile
PKG_CXXFLAGS = -DNDEBUG -DNO_GUI -DNO_AUDIO -DNO_GRAPHICS -DNO_NETWORK \
               -D_FILE_OFFSET_BITS=64 -DUNIX
```

### Praat Sources Needed

Core (minimal for intensity):
- `melder/*.c` - Memory, strings, utilities
- `sys/melder*.cpp` - System functions  
- `fon/Sound.cpp` - Sound object
- `fon/Intensity.cpp` - Intensity object
- `fon/Sound_to_Intensity.cpp` - Conversion
- `fon/Vector.cpp`, `fon/Matrix.cpp` - Base classes
- `NUM*.cpp` - Numerical utilities (FFT, Bessel, etc.)

### Library Dependencies

Praat uses external libraries:
- `fmt` - String formatting (submodule in extern/)
- Standard C++ library
- Math library (-lm)

## Files to Create/Modify

### New Files
- `src/praat_library/Makefile` - Build Praat static library
- `src/praat_wrappers.cpp` - Rcpp wrappers for Praat functions
- `R/ssff_praat_direct.R` - R wrappers for direct Praat calls
- `PRAAT_DIRECT_LINKING.md` - Documentation

### Modified Files
- `src/Makevars` - Link against libpraat.a, add include paths
- `src/superassp_init.c` - Register new functions
- `DESCRIPTION` - Update SystemRequirements if needed

## Expected Performance

Based on Parselmouth benchmarks:
- **Python/Parselmouth**: ~6 ms (baseline)
- **Direct Rcpp/Praat**: ~6 ms (same optimized code)
- **Pure C++ reimplementation**: ~20 ms (our current approach)

Expected speedup: **3-4x faster** than pure C++ implementation, **no dependency on Python**.

## Risks and Mitigations

### Risk 1: Praat Source Complexity
- **Mitigation**: Start with minimal subset (intensity only)
- **Mitigation**: Use Parselmouth's CMakeLists as guide

### Risk 2: Build System Compatibility
- **Mitigation**: Test on multiple platforms (Mac, Linux, Windows)
- **Mitigation**: Fallback to Python version if build fails

### Risk 3: Praat Memory Model
- **Mitigation**: Study Parselmouth's approach
- **Mitigation**: Use RAII wrappers (autoSound, autoIntensity)

### Risk 4: Praat Updates
- **Mitigation**: Pin to specific Parselmouth/Praat version
- **Mitigation**: Document update procedure

## Next Steps

1. ✅ Assess feasibility (this document)
2. Extract minimal Praat source set for intensity
3. Create Praat library build system
4. Implement first wrapper (intensity)
5. Benchmark and compare
6. Document for other functions
7. Systematically migrate procedures

## References

- Parselmouth: https://github.com/YannickJadoul/Parselmouth
- Praat: https://github.com/praat/praat
- Parselmouth CMake: `src/Parselmouth/CMakeLists.txt`
- Praat bindings: `src/Parselmouth/src/parselmouth/*.cpp`

---

**Status**: Strategy document - ready for implementation  
**Date**: 2025-10-18  
**Priority**: High - enables multiple procedure migrations
