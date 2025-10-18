# Direct Praat Linking: Step-by-Step Implementation Guide

*Use this guide ONLY if you decide to proceed with direct Praat linking after considering the cost-benefit analysis in PRAAT_DIRECT_LINKING_ANALYSIS.md*

## Prerequisites

- [ ] Confirmed 10+ Praat procedures need migration
- [ ] Python deployment is confirmed blocker
- [ ] 100+ hours development time allocated
- [ ] Team familiar with C/C++, CMake, and R package development

## Phase 1: Minimal Proof of Concept (2 weeks)

### Step 1.1: Identify Minimal Source Set

Create a list of required Praat sources for intensity only:

```bash
cd src/Parselmouth/praat
# Core numerical
find melder -name "*.cpp" -o -name "*.c" | grep -E "(melder\.c|NUMmath|NUMspecfunc|NUM\.cpp|VEC\.cpp|MAT\.cpp)"

# Core objects  
find sys -name "*.cpp" | grep -E "(Thing|Data|Collection|Daata)"

# Sound and Intensity
find fon -name "*.cpp" | grep -E "(Function|Sampled|Matrix|Vector|Sound|Intensity)"
```

Document each file and its dependencies in `PRAAT_SOURCES.txt`.

### Step 1.2: Create Standalone Build

Create `src/praat_build/Makefile`:

```makefile
# Praat minimal library for intensity
PRAAT_DIR = ../Parselmouth/praat

# Compiler flags
CXXFLAGS = -std=c++17 -O3 -DNDEBUG \
           -DNO_GUI -DNO_AUDIO -DNO_GRAPHICS -DNO_NETWORK \
           -D_FILE_OFFSET_BITS=64 -DUNIX

# Include paths
INCLUDES = -I$(PRAAT_DIR)/melder \
           -I$(PRAAT_DIR)/sys \
           -I$(PRAAT_DIR)/stat \
           -I$(PRAAT_DIR)/fon \
           -I../Parselmouth/extern/fmt/include

# Source files (extracted from CMakeLists)
MELDER_SRCS = $(PRAAT_DIR)/melder/melder.cpp \
              $(PRAAT_DIR)/melder/NUMmath.cpp \
              $(PRAAT_DIR)/melder/NUM.cpp \
              # ... add all needed sources

SYS_SRCS = $(PRAAT_DIR)/sys/Thing.cpp \
           $(PRAAT_DIR)/sys/Data.cpp \
           # ... add all needed sources

FON_SRCS = $(PRAAT_DIR)/fon/Function.cpp \
           $(PRAAT_DIR)/fon/Sampled.cpp \
           $(PRAAT_DIR)/fon/Vector.cpp \
           $(PRAAT_DIR)/fon/Matrix.cpp \
           $(PRAAT_DIR)/fon/Sound.cpp \
           $(PRAAT_DIR)/fon/Intensity.cpp \
           $(PRAAT_DIR)/fon/Sound_to_Intensity.cpp \
           # ... add all needed sources

ALL_SRCS = $(MELDER_SRCS) $(SYS_SRCS) $(FON_SRCS)
OBJS = $(ALL_SRCS:.cpp=.o)

# Target
libpraat_minimal.a: $(OBJS)
	ar rcs $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJS) libpraat_minimal.a
```

Test build:
```bash
cd src/praat_build
make
# Should produce libpraat_minimal.a
```

### Step 1.3: Create Minimal Rcpp Wrapper

Create `src/praat_intensity_direct.cpp`:

```cpp
#include <Rcpp.h>

// Praat headers
extern "C" {
#include "Parselmouth/praat/sys/melder.h"
}
#include "Parselmouth/praat/fon/Sound.h"
#include "Parselmouth/praat/fon/Intensity.h"
#include "Parselmouth/praat/fon/Sound_to_Intensity.h"

using namespace Rcpp;

// Helper: Convert AsspDataObj to Praat Sound
static autoSound asspobj_to_sound(SEXP audio_obj) {
  List obj(audio_obj);
  NumericMatrix audio = obj["audio"];
  double sampleRate = as<double>(obj.attr("sampleRate"));
  int nSamples = audio.nrow();
  int nChannels = audio.ncol();
  
  double duration = nSamples / sampleRate;
  double dx = 1.0 / sampleRate;
  double x1 = 0.5 * dx;
  
  autoSound sound = Sound_create(nChannels, 0.0, duration, nSamples, dx, x1);
  
  // Normalize INT16 to float
  const double norm = (std::abs(audio(0, 0)) > 10.0) ? 32768.0 : 1.0;
  
  for (int ch = 0; ch < nChannels; ch++) {
    for (int i = 0; i < nSamples; i++) {
      sound->z[ch + 1][i + 1] = audio(i, ch) / norm;
    }
  }
  
  return sound;
}

// Helper: Convert Praat Intensity to R list
static List intensity_to_list(Intensity intensity) {
  int nFrames = intensity->nx;
  double dt = intensity->dx;
  double t1 = intensity->x1;
  
  NumericVector values(nFrames);
  NumericVector times(nFrames);
  
  for (int i = 0; i < nFrames; i++) {
    values[i] = intensity->z[1][i + 1];
    times[i] = t1 + i * dt;
  }
  
  return List::create(
    Named("intensity") = values,
    Named("times") = times,
    Named("sample_rate") = 1.0 / dt,
    Named("n_frames") = nFrames,
    Named("dt") = dt,
    Named("t1") = t1
  );
}

//' Direct Praat Intensity (via Praat C library)
// [[Rcpp::export]]
List praat_intensity_direct(SEXP audio_obj,
                            double minimal_f0_frequency = 50.0,
                            double time_step = 0.0,
                            bool subtract_mean = true) {
  try {
    // Initialize Melder if needed
    static bool melder_initialized = false;
    if (!melder_initialized) {
      Melder_init();
      melder_initialized = true;
    }
    
    // Convert to Praat Sound
    autoSound sound = asspobj_to_sound(audio_obj);
    
    // Call Praat function
    autoIntensity intensity = Sound_to_Intensity(
      sound.get(),
      minimal_f0_frequency,
      time_step,
      subtract_mean
    );
    
    // Convert back to R
    return intensity_to_list(intensity.get());
    
  } catch (MelderError) {
    Melder_clearError();
    stop("Praat intensity extraction failed");
  } catch (...) {
    stop("Unknown error during intensity extraction");
  }
}
```

### Step 1.4: Update Build Configuration

Update `src/Makevars`:

```makefile
# Add to existing configuration
PRAAT_DIR = Parselmouth/praat
PRAAT_BUILD = praat_build

PKG_CPPFLAGS = -I $(PRAAT_DIR)/melder \
               -I $(PRAAT_DIR)/sys \
               -I $(PRAAT_DIR)/stat \
               -I $(PRAAT_DIR)/fon \
               -I Parselmouth/extern/fmt/include \
               -DNO_GUI -DNO_AUDIO -DNO_GRAPHICS -DNO_NETWORK

PKG_LIBS = $(PRAAT_BUILD)/libpraat_minimal.a

# Add praat_intensity_direct.cpp to CXX_SOURCES
CXX_SOURCES = RcppExports.cpp dsp_helpers.cpp ... praat_intensity_direct.cpp

# Build Praat library before R package
$(SHLIB): $(PRAAT_BUILD)/libpraat_minimal.a

$(PRAAT_BUILD)/libpraat_minimal.a:
	cd $(PRAAT_BUILD) && $(MAKE)
```

### Step 1.5: Test and Benchmark

```r
# Test
devtools::load_all()
test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
audio_obj <- av_to_asspDataObj(test_wav)

# Compare implementations
library(microbenchmark)
bench <- microbenchmark(
  python = praat_intensity(test_wav, toFile = FALSE),
  cpp_reimpl = praat_intensity_cpp_wrapper(test_wav, toFile = FALSE),
  direct = {
    obj <- av_to_asspDataObj(test_wav)
    praat_intensity_direct(obj)
  },
  times = 50
)
print(bench)
```

**Go/No-Go Decision**: If direct linking isn't faster than Python, STOP HERE.

## Phase 2: Production Integration (2 weeks)

### Step 2.1: Error Handling

Add comprehensive error handling and memory cleanup:
- Try-catch blocks around all Praat calls
- Ensure Melder errors are properly cleared
- Test with invalid inputs

### Step 2.2: R Wrapper Functions

Create `R/ssff_praat_direct.R`:

```r
praat_intensity_direct_wrapper <- function(listOfFiles, ...) {
  # Follow pattern from ssff_cpp_sptk_rapt.R
  # Load audio, call praat_intensity_direct, create AsspDataObj
}
```

### Step 2.3: Testing

Create `tests/testthat/test-praat-direct.R`:
- Compare output with Python version
- Test error cases
- Test edge cases (very short audio, very long audio)
- Test different audio formats

### Step 2.4: Documentation

- Update package documentation
- Add migration guide for new procedures
- Document build requirements

## Phase 3: Scale to Other Procedures (Ongoing)

### Adding New Procedures

For each additional Praat procedure:

1. Identify required Praat sources (check CMakeLists.txt)
2. Add sources to `praat_build/Makefile` if needed
3. Create Rcpp wrapper in `src/praat_*.cpp`
4. Create R wrapper in `R/ssff_praat_*.R`
5. Add tests
6. Benchmark vs Python version

### Maintenance

- Monitor Parselmouth updates
- Test on all target platforms
- Keep documentation current

## Platform-Specific Notes

### macOS
- Use clang++ from Xcode Command Line Tools
- May need `-stdlib=libc++`
- Test on both Intel and ARM architectures

### Linux
- Use g++ or clang++
- Ensure glibc compatibility
- Test on Ubuntu, CentOS/RHEL

### Windows
- Use Rtools compiler
- May need `/EHsc` flag
- Handle path separators correctly

## Troubleshooting

### Build fails with undefined references
- Check all required source files are included
- Verify link order (Praat lib must come before system libs)
- Check for missing Praat dependencies

### Memory leaks
- Use valgrind: `R -d valgrind -f test_script.R`
- Ensure all `auto*` objects are properly scoped
- Check that Praat errors don't leak memory

### Different results from Python
- Check audio normalization
- Verify parameter mapping
- Compare intermediate values
- Check floating-point precision

## Success Criteria

- [ ] Build succeeds on all platforms
- [ ] Results match Python/Parselmouth (< 0.2 dB difference)
- [ ] Performance equal or better than Python
- [ ] No memory leaks (valgrind clean)
- [ ] Comprehensive tests pass
- [ ] Documentation complete

## Abort Criteria

Stop and revert if:
- Build takes >2 weeks to stabilize
- Performance not better than Python
- Platform-specific issues can't be resolved
- Memory management proves problematic

---

**Note**: This is a complex undertaking. Budget 100+ hours and expect challenges. Keep Python/Parselmouth as fallback throughout development.

**Alternative**: Consider contributing to Parselmouth to add needed features rather than reimplementing.
