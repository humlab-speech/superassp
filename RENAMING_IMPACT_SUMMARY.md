# Function Renaming Impact Analysis

## Overview

Implementing the new naming convention with `trk_` and `lst_` prefixes will affect **55 functions** across the superassp package:
- **45 SSFF functions** to be renamed with `trk_` prefix
- **8 list functions** to be renamed with `lst_` prefix  
- **1 function kept as-is** (`praat_sauce`)
- **7 functions to deprecate** (slower Python implementations)

## Statistics

### Functions by Action
- **Rename**: 53 functions
- **Keep as-is**: 1 function (`praat_sauce`)
- **Deprecate**: 7 functions

### Functions by Output Type
- **SSFF (trk_)**: 45 functions + 5 to deprecate = 50 total SSFF functions
- **list (lst_)**: 8 functions + 2 to deprecate = 10 total list functions

### Functions by Implementation
- **C (ASSP)**: 14 functions → all renamed with `trk_` prefix
- **C++ (SPTK)**: 7 functions → all renamed with `trk_` prefix
- **C++ (ESTK)**: 1 function → renamed with `trk_` prefix
- **C (Snack)**: 2 functions → renamed with `trk_` prefix
- **Python (Parselmouth)**: 11 functions (6 SSFF + 4 list + 1 exception)
- **Python (Other)**: 18 functions (various pitch/feature extraction libraries)

## Breaking Changes

### User Code Impact

All user code calling these functions will break unless deprecated aliases are provided. Example:

```r
# Old code (will break)
result <- rapt(wav_file)
formants <- praat_formant_burg(wav_file)
features <- eGeMAPS(wav_file)

# New code (required)
result <- trk_rapt(wav_file)
formants <- trk_formant_burgp(wav_file)
features <- lst_egemaps(wav_file)
```

### Package Files to Update

#### 1. R Function Files (~55 files)
- All `R/ssff_*.R` files
- All `R/list_*.R` files
- Each file needs function name change + documentation updates

#### 2. Documentation (~55 files)
- All `man/*.Rd` files corresponding to renamed functions
- Examples in each help file
- Cross-references between help files

#### 3. NAMESPACE (1 file)
- Export statements for all new names
- Keep old names exported during deprecation period

#### 4. Tests (~20-30 files)
- `tests/testthat/test-*.R` files
- All function calls in tests
- Test descriptions and expectations

#### 5. Benchmarks (~10 files)
- `benchmarking/benchmark_*.R` scripts
- `benchmarking/benchmark_suite.R`
- Result analysis scripts

#### 6. README and Vignettes
- `README.md` - all code examples
- All vignette files
- Quick start guides

#### 7. Supporting Documentation
- `CLAUDE.md` - development guide
- Example code snippets
- Migration guides

## Deprecated Functions

These 7 functions will show deprecation warnings and eventually be removed:

| Function | Reason | Replacement |
|----------|--------|-------------|
| `nonopt_rapt` | 2-3x slower than C++ | `trk_rapt` |
| `nonopt_swipe` | 2-3x slower than C++ | `trk_swipe` |
| `nonopt_reaper` | 2-3x slower than C++ | `trk_reaper` |
| `aperiodicities` | Slower Python implementation | `trk_d4c` |
| `mfcc` | Slower Torch implementation | `trk_mfcc` |
| `vat` | Will be replaced with new implementation | TBD |
| `voice_analysis_toolkit` | Will be replaced with new implementation | TBD |

## Key Naming Pattern Changes

### ASSP Functions
- Remove `ana` suffix where present
- Add `trk_` prefix
- Example: `rmsana` → `trk_rms`, `acfana` → `trk_acf`

### SPTK Functions
- Remove `sptk_` prefix (redundant)
- Add `trk_` prefix
- Example: `sptk_mfcc` → `trk_mfcc`

### ESTK Functions
- Remove `estk_` prefix (redundant)
- Add `trk_` prefix
- Example: `estk_pitchmark` → `trk_pitchmark`

### Parselmouth Functions
- Remove `praat_` prefix
- Add `trk_` or `lst_` prefix
- Add `p` suffix to indicate Parselmouth
- Examples:
  - `praat_pitch` → `trk_pitchp`
  - `praat_voice_report` → `lst_voice_reportp`
- Exception: `praat_sauce` stays as-is

### OpenSmile Functions
- Add `lst_` prefix
- Normalize casing (all lowercase)
- Examples:
  - `eGeMAPS` → `lst_egemaps`
  - `ComParE_2016` → `lst_compare_2016`

### Snack Functions
- Keep `snack` in name (algorithm identifier)
- Add `trk_` prefix
- Examples:
  - `snack_pitch` → `trk_snack_pitch`
  - `snack_formant` → `trk_snack_formant`

## Migration Strategy

### Phase 1: Create Deprecation Wrappers (Low Risk)
```r
# In each old function file, add deprecation wrapper
rapt <- function(...) {
  .Deprecated("trk_rapt", package = "superassp",
              msg = "rapt() is deprecated. Please use trk_rapt() instead.")
  trk_rapt(...)
}
```

### Phase 2: Create New Functions
- Copy function files to new names
- Update function definitions
- Update internal calls
- Update documentation

### Phase 3: Update Tests
- Create new test files with new names
- Update all function calls
- Keep old tests working via deprecation wrappers

### Phase 4: Update Documentation
- Update all examples
- Update README
- Create migration vignette
- Update CLAUDE.md

### Phase 5: Update Benchmarks
- Update all benchmark scripts
- Regenerate benchmark results
- Update README with new names

### Phase 6: Deprecation Period (1-2 releases)
- Old names work but warn
- Give users time to update

### Phase 7: Remove Old Names (Future release)
- Remove deprecation wrappers
- Remove old function definitions
- Clean up NAMESPACE

## Estimated Effort

### By Task
- **Function renaming**: ~8 hours (systematic but tedious)
- **Documentation updates**: ~6 hours (all .Rd files + examples)
- **Test updates**: ~4 hours (update all test calls)
- **Benchmark updates**: ~2 hours (update scripts + regenerate)
- **README/vignettes**: ~2 hours (update examples)
- **Testing/validation**: ~4 hours (ensure everything works)

**Total estimate**: ~26 hours of focused work

### Parallelization Opportunities
- Function renaming can be done in batches by type (ASSP, SPTK, etc.)
- Documentation updates can be semi-automated with regex
- Tests can be updated in parallel with functions

## Risk Assessment

### High Risk
- Breaking all user code without deprecation period
- Missing function references in documentation
- Breaking cross-package dependencies (if any)

### Medium Risk
- Inconsistent renaming (missing some references)
- Broken examples in documentation
- Failed tests due to missed updates

### Low Risk
- Performance changes (none - only names changing)
- Functionality changes (none - only names changing)

## Recommendation

**Recommended Approach**: Implement with 2-release deprecation period

1. **Version 2.0.0**: Introduce new names + deprecation wrappers
   - All new `trk_*` and `lst_*` functions work
   - Old names work but show warnings
   - Update all documentation to use new names
   - Provide migration guide

2. **Version 2.1.0**: Continue deprecation warnings
   - Monitor for issues
   - Help users migrate
   - Fix any problems discovered

3. **Version 3.0.0**: Remove old names
   - Clean break
   - Remove all deprecation wrappers
   - Smaller, cleaner codebase

## Benefits

1. **Clarity**: Immediate visual indication of return type (`trk_` vs `lst_`)
2. **Consistency**: All functions follow same pattern
3. **Discoverability**: Tab completion in IDE groups related functions
4. **Documentation**: Easier to organize and find functions
5. **Performance**: Deprecating slower Python versions guides users to faster C++ implementations

## Compatibility Notes

### Backward Compatibility
- Maintain via `.Deprecated()` wrappers
- No performance penalty (just wraps new function)
- Clear migration path for users

### Forward Compatibility
- New code uses new names exclusively
- Cleaner namespace
- Better long-term maintainability

## Communication Plan

### Announcement
- NEWS.md with complete change list
- README with prominent notice
- Package startup message (temporary)
- Email to known users/maintainers

### Documentation
- Migration vignette with examples
- Table of old→new name mappings
- FAQ section for common issues

### Support
- Monitor GitHub issues closely
- Quick responses during transition
- Consider blog post or tutorial
