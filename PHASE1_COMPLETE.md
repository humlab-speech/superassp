# Phase 1 Complete: Track Name Migration Documentation

**Status**: ✅ Complete
**Date**: 2025-10-22
**Phase**: Documentation & Mapping (No Code Changes)

## Summary

Phase 1 of the track name migration to Titze 2015 / Nylén 2024 standards is complete. All documentation, mapping, and validation tools have been created **without any code changes**.

## Deliverables

### 1. Master Planning Document
**File**: `TRACK_NAMING_MIGRATION_PLAN.md` (45 pages)

**Contents**:
- Complete migration strategy (3 phases)
- Standards from Titze 2015 and Nylén 2024
- Analysis of current state (149 track definitions)
- Detailed mapping tables (8 categories)
- Implementation checklist
- Risk assessment and mitigation
- Timeline proposal

**Key Statistics**:
- 149 track definitions analyzed
- ~50+ functions affected
- 8 categories identified
- 3-phase migration strategy proposed

### 2. Track Name Mapping
**File**: `TRACK_NAMES_MAPPING.csv` (112 entries)

**Contents**:
- Complete mapping: old name → new name
- Function associations
- Category classifications
- Unit specifications
- Change notes

**Categories Covered**:
- F0 / frequency of oscillation (24 tracks)
- Formants (2 tracks)
- Bandwidths (2 tracks)
- Harmonics (VoiceSauce)
- Harmonic differences (corrected/uncorrected)
- Voice quality measures (11 tracks)
- Spectral measures (10 tracks)
- Metadata (16 tracks)
- Other (45 tracks)

**Key Findings**:
- 40 tracks need changes (36.0%)
- 71 tracks unchanged (64.0%)
- Most common changes: adding units, case normalization

### 3. Comprehensive Inventory
**File**: `TRACK_INVENTORY.md`
**Script**: `scripts/generate_track_inventory.R`

**Contents**:
- Summary statistics by category
- Function breakdown (top 20)
- Detailed category tables
- Complete mapping appendix

**Top Functions** (by track count):
1. `lst_voice_reportp` - 26 tracks
2. `trk_spectral_momentsp` - 4 tracks
3. `arfana`, `larana`, `lpcana` - 3 tracks each

### 4. Validation System
**Script**: `scripts/validate_track_names.R`

**Features**:
- Validates proposed names against standards
- Checks Titze 2015 compliance
- Checks Nylén 2024 extensions
- Verifies unit consistency
- Category-specific validation rules

**Results**:
- ✅ All 111 proposed track names pass validation (100%)
- 0 issues found
- Standards checked: Titze 2015, Nylén 2024, explicit units

**Validation Rules**:
- F0: lowercase 'f0', with [Hz] unit
- Formants: uppercase 'F' + number + [Hz]
- Bandwidths: uppercase 'B' + number + [Hz]
- Harmonics: Hi or Ai format + [dB]
- Harmonic differences: hyphen separator + suffix + [dB]
- Voice quality: appropriate units ([dB], [%], [us])

### 5. VoiceSauce Equivalents
**File**: `VOICESAUCE_EQUIVALENTS.md`

**Contents**:
- Complete mapping of VoiceSauce parameters
- 40+ parameter equivalents documented
- Migration guide for VoiceSauce users
- Backwards compatibility strategy
- Code examples (old vs new workflow)

**Key Mappings**:
- fo variants → `fo[Hz]`
- Formants F1-F4 → `F1[Hz]`-`F4[Hz]`
- Harmonics H1, A1, etc. → `H1[dB]`, `A1[dB]`, etc.
- Differences H1H2c → `H1-H2c[dB]` (with hyphen)
- Quality measures → all get appropriate units

### 6. emuR Compatibility Analysis
**File**: `EMUR_COMPATIBILITY_ANALYSIS.md`

**Contents**:
- Impact analysis on emuR integration
- SSFF file format compatibility
- Migration scenarios (3 types)
- Proposed migration tools
- Risk assessment
- Testing plan
- Implementation checklist

**Key Findings**:
- ✅ SSFF format is flexible - supports any track names
- ✅ No structural breaking changes
- ⚠️ Existing databases require migration
- ✅ Migration tools will ease transition

**Proposed Tools**:
- `check_emur_database()` - scan for old names
- `migrate_emur_database()` - update track names
- `verify_emur_database()` - confirm migration

**Risk Level**: Low - migration is fully backwards compatible

## Standards Compliance

### Titze 2015 (JASA)
**"Toward a consensus on symbolic notation..."**

✅ Compliant:
- `fo` (lowercase) for fundamental frequency
- `Fi` (uppercase F, subscript i) for formants
- `Bi` (uppercase B, subscript i) for bandwidths
- `Hi` for harmonics, `Ai` for harmonic at Fi
- Hyphen separator for differences (H1-H2)

### Nylén 2024 (JASA)
**"Acoustic cues to femininity and masculinity..."**

✅ Compliant:
- Suffix `c` for corrected measures (H1-H2c)
- Suffix `u` for uncorrected measures (H1-H2u)
- Enables gender/voice quality research

### superassp Conventions

✅ Maintained:
- Explicit units in brackets: `[Hz]`, `[dB]`, `[%]`
- Integration with R `units` package
- Consistent notation across all functions

## Example Transformations

### Before (Old Names)
```r
attr(trk_rapt, "tracks") <- c("fo")
attr(trk_forest, "tracks") <- c("F[Hz]", "B[Hz]")
attr(trk_praat_sauce, "tracks") <- c("fo", "F1", "F2", "H1H2c")
```

### After (New Names)
```r
attr(trk_rapt, "tracks") <- c("fo[Hz]")
attr(trk_forest, "tracks") <- c("F1[Hz]", "F2[Hz]", "F3[Hz]", "F4[Hz]",
                                "B1[Hz]", "B2[Hz]", "B3[Hz]", "B4[Hz]")
attr(trk_praat_sauce, "tracks") <- c("fo[Hz]", "F1[Hz]", "F2[Hz]", "H1-H2c[dB]")
```

## Key Benefits

1. **Scientific Compliance**: Published JASA standards
2. **Cross-Software Compatibility**: Matches Praat, VoiceSauce
3. **Clarity**: Explicit units prevent confusion
4. **Automatic Unit Conversion**: R `units` package integration
5. **Research Support**: Corrected formants for gender/voice analysis
6. **Consistency**: Unified notation across 149 definitions
7. **Citability**: Reference Titze 2015 and Nylén 2024 in papers

## Next Steps (Phase 2)

**Goal**: Implement backwards-compatible migration

**Tasks**:
1. Implement `.track_name_aliases()` function
2. Update `[[.AsspDataObj` with alias support
3. Update all `attr(*, "tracks")` definitions
4. Add deprecation warnings
5. Create migration utilities for emuR
6. Update tests and documentation
7. Release v0.7.0 with deprecation period

**Timeline**: 1-2 months

## Statistics Summary

| Metric | Count | Percentage |
|--------|-------|------------|
| Total track definitions | 111 | 100% |
| Tracks needing changes | 40 | 36.0% |
| Tracks unchanged | 71 | 64.0% |
| Functions affected | 53 | — |
| Categories | 8 | — |
| VoiceSauce parameters | 40+ | — |
| Validation pass rate | 111/111 | 100% |

## Files Created (Phase 1)

1. `TRACK_NAMING_MIGRATION_PLAN.md` - 45-page master plan
2. `TRACK_NAMES_MAPPING.csv` - 112 entries
3. `TRACK_INVENTORY.md` - Complete inventory
4. `VOICESAUCE_EQUIVALENTS.md` - VoiceSauce mapping
5. `EMUR_COMPATIBILITY_ANALYSIS.md` - emuR analysis
6. `scripts/generate_track_inventory.R` - Inventory generator
7. `scripts/validate_track_names.R` - Validation tool
8. `PHASE1_COMPLETE.md` - This summary

**Total**: 8 files, 0 code changes, 100% documentation

## Validation Results

✅ **All proposed track names validated successfully**

- Titze 2015 compliance: 100%
- Nylén 2024 compliance: 100%
- Unit specification: 100%
- Naming consistency: 100%

## Impact Assessment

### Low Risk
- ✅ No immediate breaking changes
- ✅ Long deprecation period planned
- ✅ Alias system for backwards compatibility
- ✅ Migration tools for emuR databases

### Medium Priority
- ⚠️ Requires user communication
- ⚠️ Documentation updates needed
- ⚠️ Coordination with emuR maintainers

### High Value
- ✅ Scientific standard compliance
- ✅ Improved cross-platform compatibility
- ✅ Better research support
- ✅ Future-proof notation

## References

1. Titze, I. R., et al. (2015). "Toward a consensus on symbolic notation of harmonics, resonances, and formants in vocalization." *JASA*, 137(5), 3005-3007.

2. Nylén, F., et al. (2024). "Acoustic cues to femininity and masculinity in spontaneous speech." *JASA*, 155(2), 1373-1387.

3. VoiceSauce Documentation: https://www.phonetics.ucla.edu/voicesauce/documentation/parameters.html

4. emuR Manual: https://ips-lmu.github.io/The-EMU-SDMS-Manual/

## Approval for Phase 2

Phase 1 is complete and ready for review. Approval to proceed to Phase 2 (implementation) requires:

- [ ] Review of all documentation
- [ ] Approval of migration strategy
- [ ] Decision on timeline
- [ ] Coordination with emuR maintainers (optional but recommended)

**Status**: ✅ Ready for Phase 2 implementation
