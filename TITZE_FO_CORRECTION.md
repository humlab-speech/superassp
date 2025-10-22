# CRITICAL CORRECTION: Titze 2015 Uses f_o (NOT f_0)

**Date**: 2025-10-22
**Status**: ✅ All Phase 1 documents updated

## The Correction

Titze et al. (2015) JASA paper uses **f_o** (subscript letter 'o') for fundamental frequency, NOT **f_0** (subscript digit zero).

### Notation Meaning

- **f_o** = **f**requency of **o**scillation (semantic meaning)
- **NOT f_0** = arbitrary zero subscript

### Why This Matters

1. **Semantic clarity**: 'o' explicitly means "oscillation"
2. **Distinguishes from harmonics**: Harmonics use numeric subscripts (H₁, H₂, H₃)
3. **Titze 2015 rationale**: Meaningful subscripts preferred over arbitrary numbering

## Impact on superassp

### Good News: ASSP Functions Already Correct! ✅

The ASSP library functions were **already using correct Titze notation**:

```r
attr(trk_ksvfo, "tracks") <- "fo[Hz]"  # ✅ Correct!
attr(fo, "tracks") <- "fo[Hz]"         # ✅ Correct!
attr(foana, "tracks") <- "fo[Hz]"      # ✅ Correct!
attr(fo_ksv, "tracks") <- "fo[Hz]"     # ✅ Correct!
```

### What Needs Updating

Only Python/C++ pitch tracking functions need correction:

**Before (incorrect)**:
```r
attr(trk_rapt, "tracks") <- "f0"       # ✗ Wrong
attr(trk_swipe, "tracks") <- "f0"      # ✗ Wrong
attr(trk_swiftf0, "tracks") <- "F0"    # ✗ Wrong
```

**After (correct)**:
```r
attr(trk_rapt, "tracks") <- "fo[Hz]"   # ✓ Correct
attr(trk_swipe, "tracks") <- "fo[Hz]"  # ✓ Correct
attr(trk_swiftf0, "tracks") <- "fo[Hz]" # ✓ Correct
```

## Updated Standards

### Fundamental Frequency

| Current | SSFF/Attr (Titze 2015) | Data Frame (R-friendly) | Status |
|---------|------------------------|-------------------------|--------|
| `fo[Hz]` | `fo[Hz]` | `fo_Hz` | ✅ Already correct (ASSP) |
| `f0` | `fo[Hz]` | `fo_Hz` | ⚠️ Needs update (Python/C++) |
| `F0` | `fo[Hz]` | `fo_Hz` | ⚠️ Needs update (Swift) |

### Complete Titze 2015 Notation System

- **f_o**: Frequency of oscillation (fundamental)
- **F_i** (i=1,2,3...): Frequencies of resonance (formants)
- **B_i** (i=1,2,3...): Bandwidths of resonances  
- **H_i** (i=1,2,3...): Harmonics (integer multiples of f_o)
- **A_i** (i=1,2,3...): Amplitude of harmonic nearest F_i

## Functions Affected

### Already Correct (4 functions - no change needed)
- `trk_ksvfo` → `fo[Hz]` ✅
- `fo` → `fo[Hz]` ✅
- `foana` → `fo[Hz]` ✅
- `fo_ksv` → `fo[Hz]` ✅
- `trk_mhspitch` → `fo[Hz]` ✅

### Need Update (20 functions - f0 → fo[Hz])
Python implementations:
- `trk_rapt`, `trk_swipe`, `trk_dio`, `trk_reaper`, `trk_harvest`
- `trk_yin`, `trk_pyin`, `trk_crepe`, `trk_yaapt`, `trk_torch_pitch`
- `trk_kaldi_pitch`, `trk_snackp`, `trk_npy_import`
- `trk_world_harvest`, `trk_world_dio`
- `nonopt_rapt`, `nonopt_swipe`, `nonopt_reaper`

C++ implementations:
- Various SPTK wrappers

Other:
- `trk_swiftf0` (F0 → fo[Hz])
- `trk_covarep_srh` (F0[Hz] → fo[Hz])
- `trk_praat_sauce` (f0 → fo[Hz])

## Benefits of Correct Notation

1. **Semantic meaning**: "oscillation" not arbitrary zero
2. **Distinguishes roles**: 
   - f_o = oscillation (fundamental)
   - H_0 would conflict (first harmonic IS f_o)
   - H_1, H_2 = integer multiples of f_o
3. **Less work**: Most speech functions already correct
4. **Standards compliant**: Exact Titze 2015 notation

## Example: Titze Notation System

```
f_o = 100 Hz        (frequency of oscillation)
H_1 = 100 Hz        (first harmonic = f_o)
H_2 = 200 Hz        (second harmonic = 2×f_o)
H_3 = 300 Hz        (third harmonic = 3×f_o)

F_1 = 800 Hz        (first formant)
F_2 = 1200 Hz       (second formant)
F_3 = 2400 Hz       (third formant)

A_1 = amplitude of H near F_1  (harmonic amplitude at ~800 Hz)
A_2 = amplitude of H near F_2  (harmonic amplitude at ~1200 Hz)
A_3 = amplitude of H near F_3  (harmonic amplitude at ~2400 Hz)
```

## What Changed in Phase 1 Documents

All Phase 1 documents have been updated:

1. ✅ `TRACK_NAMING_MIGRATION_PLAN.md` - f0 → fo throughout
2. ✅ `TRACK_NAMES_MAPPING.csv` - Updated 24 entries  
3. ✅ `TRACK_INVENTORY.md` - Regenerated with fo
4. ✅ `VOICESAUCE_EQUIVALENTS.md` - Updated all F0 references
5. ✅ `NAMING_AESTHETICS_ANALYSIS.md` - Examples use fo
6. ✅ `PHASE1_COMPLETE.md` - Summary corrected
7. ✅ `scripts/validate_track_names.R` - Validation rules updated
8. ✅ `scripts/generate_track_inventory.R` - Regenerated

## Migration Impact

### Reduced Complexity

**Original plan** (incorrect):
- Change ASSP functions: `fo` → `f0` (4 functions)
- Change Python/C++: `f0` → `f0` (no change, 20 functions)
- Total changes: 4 functions

**Corrected plan**:
- Keep ASSP functions: `fo[Hz]` (no change, 4 functions) ✅
- Change Python/C++: `f0` → `fo[Hz]` (20 functions)
- Total changes: 20 functions

**But**: Speech-specific ASSP functions are already compliant!

## References

Titze, I. R., et al. (2015). "Toward a consensus on symbolic notation of harmonics, resonances, and formants in vocalization." *The Journal of the Acoustical Society of America*, 137(5), 3005–3007. https://doi.org/10.1121/1.4919349

**Key quote**: "The fundamental frequency is denoted f_o (frequency of oscillation)"

## Conclusion

This correction **simplifies** the migration:
- ASSP functions already comply with Titze 2015
- Only Python/C++ wrappers need updating
- Semantic meaning is clearer (oscillation vs arbitrary zero)
- Total work is similar but better motivated

**All Phase 1 documents have been corrected to use `fo` (not `f0`).**
