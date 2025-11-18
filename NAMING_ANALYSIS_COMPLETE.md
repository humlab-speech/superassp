# Function Naming Convention Analysis - COMPLETE

## Analysis Completed

A comprehensive analysis of all DSP functions in the superassp package has been completed according to the new naming convention requirements.

## New Naming Convention

### Rules
1. **SSFF-returning functions**: Prefix `trk_` (for "track-based data")
2. **List-returning functions**: Prefix `lst_` (for "feature lists")
3. **Framework identifiers removed**: `praat`, `sptk`, `estk` removed from names
4. **Exception**: `praat_sauce` keeps its original name
5. **Parselmouth implementations**: Add `p` suffix (e.g., `trk_pitchp`)
6. **Deprecate slow implementations**: Python versions superseded by C++ versions
7. **Deprecate VAT**: `vat` and `voice_analysis_toolkit` will be replaced

## Analysis Results

### Total Functions: 61
- **53 functions to rename**
- **7 functions to deprecate**  
- **1 function to keep as-is** (`praat_sauce`)

### By Implementation Type
- **ASSP (C)**: 14 functions → rename with `trk_` prefix
- **SPTK (C++)**: 7 functions → rename with `trk_` prefix
- **ESTK (C++)**: 1 function → rename with `trk_` prefix
- **Snack (C)**: 2 functions → rename with `trk_` prefix
- **Parselmouth (Python)**: 10 functions → 6 SSFF + 4 list
- **OpenSmile (Python)**: 4 functions → rename with `lst_` prefix
- **Other Python**: 15 functions → various pitch/feature extractors

### Aliases (No Action Required)
The following are aliases and don't require separate renaming:
- `fo`, `foana`, `fo_ksv` (aliases for `ksvfo` → `trk_ksvf0`)
- `pitch`, `pitch_mhs` (aliases for `mhspitch` → `trk_mhsf0`)

## Key Examples

### Most Common Renames
```r
# Pitch tracking
rapt()               → trk_rapt()
swipe()              → trk_swipe()
reaper()             → trk_reaper()
dio()                → trk_dio()
harvest()            → trk_harvest()
praat_pitch()        → trk_pitchp()
ksvfo()              → trk_ksvf0()
mhspitch()           → trk_mhsf0()

# Formant analysis
praat_formant_burg() → trk_formant_burgp()
snack_formant()      → trk_snack_formant()
forest()             → trk_forest()

# Spectral features
sptk_mfcc()          → trk_mfcc()
praat_intensity()    → trk_intensityp()
cepstrum()           → trk_cepstrum()

# Voice quality features
praat_voice_report() → lst_voice_reportp()
praat_avqi()         → lst_avqip()
eGeMAPS()            → lst_egemaps()
d4c()                → trk_d4c()
```

### Functions to Deprecate
```r
nonopt_rapt()           → DEPRECATED, use trk_rapt()
nonopt_swipe()          → DEPRECATED, use trk_swipe()
nonopt_reaper()         → DEPRECATED, use trk_reaper()
aperiodicities()        → DEPRECATED, use trk_d4c()
mfcc()                  → DEPRECATED, use trk_mfcc()
vat()                   → DEPRECATED, will be replaced
voice_analysis_toolkit() → DEPRECATED, will be replaced
```

## Documentation Created

The following comprehensive documentation has been created:

1. **NAMING_CONVENTION_ANALYSIS.md** (315 lines)
   - Detailed analysis with complete function inventory
   - Implementation patterns and examples
   - Phase-by-phase implementation checklist

2. **NAMING_CONVENTION_MAPPING.csv** (55 rows)
   - Machine-readable mapping of old → new names
   - Includes action (rename/deprecate/keep)
   - Includes implementation type and library

3. **RENAMING_IMPACT_SUMMARY.md** (7,895 characters)
   - Impact analysis and risk assessment
   - Estimated effort (26 hours)
   - Migration strategy
   - Communication plan

4. **RENAMING_QUICK_REFERENCE.md** (4,000+ characters)
   - User-friendly migration guide
   - Alphabetical lookup table
   - Search & replace patterns for automated migration
   - Common usage examples

## Files Requiring Updates

### R Code (~55 files)
- All `R/ssff_*.R` files
- All `R/list_*.R` files

### Documentation (~55 files)
- All corresponding `man/*.Rd` files

### Tests (~30 files)
- All `tests/testthat/test-*.R` files

### Benchmarks (~10 files)
- All `benchmarking/*.R` scripts

### Other Documentation
- README.md
- CLAUDE.md
- Vignettes
- NAMESPACE

**Total: ~200+ files**

## Implementation Impact

### Effort Estimate: 26 hours
- Function renaming: 8 hours
- Documentation updates: 6 hours
- Test updates: 4 hours
- Benchmark updates: 2 hours
- README/vignettes: 2 hours
- Testing/validation: 4 hours

### Breaking Changes
All user code calling renamed functions will break unless deprecation wrappers are provided.

### Recommended Strategy
Implement with **3-phase approach** over multiple releases:
1. **v2.0.0**: New names + deprecation wrappers
2. **v2.1.0**: Continue deprecation period
3. **v3.0.0**: Remove old names

## Benefits

1. **Type Clarity**: Return type immediately visible (`trk_` vs `lst_`)
2. **Consistency**: Uniform naming across all DSP functions
3. **Discoverability**: Related functions grouped in IDE autocomplete
4. **Performance Guidance**: Deprecating slow versions guides users to fast C++ implementations
5. **Maintainability**: Cleaner, more organized codebase

## Next Steps

To proceed with implementation:

1. Review the analysis documents
2. Confirm the naming convention approach
3. Begin Phase 1 implementation (create new functions + deprecation wrappers)
4. Update tests and documentation
5. Run full test suite
6. Update README and benchmarks
7. Release as version 2.0.0 with deprecation warnings

## Questions to Confirm

Before proceeding, please confirm:

1. ✓ Use `trk_` and `lst_` prefixes (with underscores)
2. ✓ Keep `praat_sauce` unchanged
3. ✓ Add `p` suffix for Parselmouth implementations
4. ✓ Deprecate slower Python implementations
5. ✓ Deprecate VAT functions
6. ? Preferred deprecation period (1-2 releases recommended)
7. ? Target version number for breaking changes (3.0.0 suggested)

## Status

**✅ ANALYSIS COMPLETE**

Ready to proceed with implementation when approved.
