# SPTK C++ Integration Status

## Summary
Integration of SPTK C++ pitch tracking algorithms (RAPT, SWIPE, REAPER, DIO/Harvest) was attempted but encountered significant dependency complexity.

## Work Completed

### 1. Files Created/Copied

#### C++ Wrapper
- Created `/src/sptk_pitch.cpp` (484 lines)
  - `rapt_cpp()` - RAPT pitch tracking
  - `swipe_cpp()` - SWIPE pitch tracking
  - `reaper_cpp()` - REAPER pitch tracking
  - `dio_cpp()` - DIO pitch tracking
  - All functions accept `AsspDataObj` and return R lists with F0 contours

#### SPTK Source Files
- Analysis: 6 pitch extraction implementation files (`.cc`)
- Analysis headers: 7 pitch extraction headers (`.h`)
- Math: 37 mathematical operation files copied
- Utils: 1 utility file (sptk_utils.cc)
- Compression: ~20 compression/quantization files
- Third-party: 25 files (REAPER, Snack, SWIPE, WORLD)

#### Build Configuration
- Updated `/src/Makevars` to include:
  - SPTK include paths
  - REAPER and WORLD subdirectory includes
  - Math, utils, and compression source lists

### 2. Issues Encountered

#### Deep Dependency Chain
SPTK pitch extraction has cascading dependencies:

```
PitchExtraction classes
  ├── RealValuedFastFourierTransform
  │     ├── FastFourierTransform
  │     └── sptk_utils (IsPowerOfTwo, etc.)
  ├── SymmetricMatrix
  │     └── Matrix operations
  ├── LindeBuzoGrayAlgorithm (for clustering)
  │     └── Compression/VQ modules
  └── (many more...)
```

Each resolved symbol led to 1-3 new missing symbols.

#### Build Errors Progression
1. ✅ Fixed: Missing pitch extraction `.cc` files
2. ✅ Fixed: REAPER include paths (`-I SPTK/third_party/REAPER`)
3. ✅ Fixed: WORLD include paths (`-I SPTK/third_party/WORLD`)
4. ✅ Fixed: Missing `RealValuedFastFourierTransform` (copied FFT files)
5. ✅ Fixed: Missing `IsPowerOfTwo` (copied sptk_utils.cc)
6. ✅ Fixed: Missing `SymmetricMatrix` (copied all 37 math `.cc` files)
7. ❌ **BLOCKED**: Missing `LindeBuzoGrayAlgorithm` (requires compression module)
8. ❌ **Unknown**: Likely more missing dependencies after compression

The official SPTK repository contains **327 `.cc` source files**, suggesting that pitch extraction requires a substantial portion of the library.

## Architectural Challenge

The SPTK library is **monolithic** - it was designed as a complete toolkit where components share infrastructure. The pitch extraction classes are not standalone modules.

### Alternative Approaches

#### Option 1: Include Entire SPTK Library
- **Pros**: Guaranteed to resolve all dependencies
- **Cons**:
  - Adds ~300 source files to package
  - Significantly increases compilation time (currently ~2-3 minutes, would likely be 5-10+ minutes)
  - Increases package size substantially
  - Maintenance burden (tracking SPTK upstream changes)

#### Option 2: Continue Python Integration
- **Pros**:
  - Already working (rapt, swipe, reaper via pysptk)
  - No dependency management issues
  - Easier to maintain
- **Cons**:
  - Requires Python environment
  - Performance overhead (though Python implementations use C underneath)
  - Complex setup for users

#### Option 3: Use ESTK Instead
- **Pros**:
  - ✅ **Already successfully integrated** (estk_pda_cpp, estk_pitchmark_cpp)
  - Self-contained, minimal dependencies
  - Fast native C++ performance
- **Cons**:
  - Different algorithms than SPTK
  - May have different accuracy characteristics

#### Option 4: Minimal SPTK Subset (Time-Intensive)
- Manually trace and copy only required files
- **Estimate**: 50-100 source files still needed based on dependency patterns observed
- **Time required**: 4-8 hours of dependency chasing + testing
- **Risk**: High chance of discovering more dependencies during testing

## Recommendation

**Prioritize Option 3 (ESTK)** combined with keeping existing Python SPTK integration:

1. **ESTK provides**:
   - Pitch tracking (Period Detection Algorithm)
   - Pitch marking (Laryngograph-style epoch detection)
   - ✅ Already integrated and working
   - Fast native C++ performance

2. **Python/pysptk provides**:
   - RAPT, SWIPE, REAPER for users who want those specific algorithms
   - Works without additional integration effort

3. **Future work** (if needed):
   - Consider full SPTK integration only if there's strong user demand
   - Alternatively, create lightweight standalone implementations of specific algorithms

## Files for Reference

### Main wrapper
- `/src/sptk_pitch.cpp` - Ready to use if dependencies resolved

### SPTK sources copied
- `/src/SPTK/src/analysis/*.cc` (6 files)
- `/src/SPTK/src/math/*.cc` (37 files)
- `/src/SPTK/src/utils/*.cc` (1 file)
- `/src/SPTK/src/compression/*.cc` (~20 files)
- `/src/SPTK/third_party/` (REAPER, Snack, SWIPE, WORLD - 25 files)

### Headers
- `/src/SPTK/include/SPTK/analysis/*.h` (7 files)
- `/src/SPTK/include/SPTK/math/*.h` (37 files)
- `/src/SPTK/include/SPTK/utils/*.h` (3 files)
- `/src/SPTK/include/SPTK/compression/*.h` (~20 files)

## Conclusion

SPTK C++ integration requires essentially the entire SPTK library due to deep interdependencies. The current package already has:
- ✅ Working ESTK C++ integration (fast, native)
- ✅ Working Python SPTK integration (functional, complete)

The cost-benefit analysis suggests **not** proceeding with full SPTK C++ integration at this time unless there's a compelling performance or functionality requirement that cannot be met by ESTK or the Python integration.
