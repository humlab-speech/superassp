# Pladdrr Integration Migration Status

**Date**: 2026-02-06  
**Version**: 0.11.0  
**Branch**: pladdrr-integration

## Overview

Complete migration from Python's parselmouth to R's pladdrr for all Praat-based functionality.

## Phase 1: Setup ✅ COMPLETE

- [x] Updated DESCRIPTION to v0.11.0
- [x] Added pladdrr (>= 4.8.16) to Imports
- [x] Created `R/install_pladdrr.R` with install/info helpers
- [x] Verified pladdrr v4.8.16 available

## Phase 2: Infrastructure ✅ COMPLETE

- [x] Created `R/pladdrr_helpers.R`
  - `av_to_pladdrr_sound()` - convert av audio to pladdrr Sound
  - `av_load_for_pladdrr()` - complete workflow
  - `get_pladdrr_ptr()` - extract C pointer
  - `pladdrr_df_to_superassp()` - format conversion

## Phase 3: Migrate Existing Functions ⏳ IN PROGRESS

### Track Functions (6 files)

1. **trk_pitchp** - `R/ssff_python_pm_ppitch.R`
   - Source: `plabench/R_implementations/pitch.R`
   - Status: 🔄 TODO
   - Priority: HIGH
   - Notes: CC/AC methods, SPINET/SHS unavailable

2. **trk_formantp** - `R/ssff_python_pm_pformantb.R`
   - Source: `plabench/R_implementations/formant.R`
   - Status: 🔄 TODO
   - Priority: HIGH
   - Notes: Formant bug fixed in v4.8.16, verify

3. **trk_formantpathp** - `R/ssff_python_pm_pformantpathb.R`
   - Source: Merge into trk_formantp with track_formants parameter
   - Status: 🔄 TODO
   - Priority: MEDIUM
   - Notes: HMM tracking may have limitations

4. **trk_praatsaucep** - `R/ssff_python_pm_psauce.R`
   - Source: `plabench/R_implementations/praatsauce.R`
   - Status: 🔄 TODO
   - Priority: HIGH
   - Notes: Complex - F0, formants, harmonics, HNR, CPP

5. **trk_spectral_momentsp** - `R/ssff_python_pm_pspectral_moments.R`
   - Source: `plabench/R_implementations/spectral_moments.R`
   - Status: 🔄 TODO
   - Priority: LOW
   - Notes: Direct port

6. **trk_intensityp** - `R/ssff_python_pm_pintensity.R`
   - Source: `plabench/R_implementations/intensity.R`
   - Status: 🔄 TODO
   - Priority: MEDIUM
   - Notes: Simple, use to_intensity_direct()

### Summary Functions (4 files)

7. **lst_avqip** - `R/list_python_pm_pavqi.R`
   - Source: `plabench/R_implementations/avqi.R`
   - Status: 🔄 TODO
   - Priority: HIGH
   - Notes: 3x faster than Python! v2.03 & v3.01

8. **lst_dsip** - `R/list_python_pm_pdsi.R`
   - Source: `plabench/R_implementations/dsi.R`
   - Status: 🔄 TODO
   - Priority: HIGH
   - Notes: Production ready

9. **lst_voice_tremorp** - `R/list_python_pm_pvoice_tremor.R`
   - Source: `plabench/R_implementations/tremor.R`
   - Status: 🔄 TODO
   - Priority: MEDIUM
   - Notes: 18 measures, validated

10. **lst_voice_reportp** - `R/list_python_pm_pvoice_report.R`
    - Source: `plabench/R_implementations/voice_report.R`
    - Status: 🔄 TODO
    - Priority: MEDIUM
    - Notes: Jitter, shimmer, HNR

## Phase 4: New Functions 📋 PENDING

11. **trk_cpps** - NEW from `plabench/R_implementations/cpp.R`
    - Status: 📋 TODO
    - Notes: Cepstral Peak Prominence

12. **trk_vuv** - NEW from `plabench/R_implementations/vuv.R`
    - Status: 📋 TODO
    - Notes: Voice/Unvoiced detection

13. **lst_vq** - NEW from `plabench/R_implementations/vq.R`
    - Status: 📋 TODO
    - Notes: Voice quality measures

14. **lst_pharyngeal** - NEW from `plabench/R_implementations/pharyngeal.R`
    - Status: 📋 TODO
    - Notes: Pharyngeal analysis

15. **utils_momel** - NEW from `plabench/R_implementations/momel_pure_r.R`
    - Status: 📋 TODO
    - Notes: Pure R, check if already present

16. **utils_intsint** - NEW from `plabench/R_implementations/intsint_pure_r.R`
    - Status: 📋 TODO
    - Notes: Pure R, check if already present

17. **lst_dysprosody** - UPDATE from `plabench/R_implementations/dysprosody.R`
    - Status: 📋 TODO
    - Notes: Migrate parselmouth → pladdrr

18. **Additional utils** - From plabench shared code
    - Status: 📋 TODO

## Phase 5: Cleanup 🧹 PENDING

- [ ] Delete `R/parselmouth_helpers.R`
- [ ] Delete `R/utils_av_parselmouth_helpers.R`
- [ ] Delete `tests/test_parselmouth_equivalence.R`
- [ ] Delete Python scripts from `inst/python/`:
  - praat_pitch.py
  - praat_formant_burg.py
  - praat_formant_path_burg.py
  - praat_intensity.py
  - praat_spectral_moments.py
  - praat_sauce.py
  - AVQI/DSI/tremor/voice_report Python scripts
- [ ] Update imports in affected files
- [ ] Search/replace parselmouth references in docs

## Phase 6: Testing 🧪 PENDING

- [ ] Create `tests/testthat/test-pladdrr-integration.R`
- [ ] Port 3-way validation from plabench
- [ ] Test each function:
  - Single file
  - Multiple files
  - Time windowing
  - toFile=TRUE/FALSE
  - Various media formats
  - Parameter edge cases
- [ ] Test emuR compatibility
- [ ] Test S7 AVAudio dispatch
- [ ] Test parallel processing

## Phase 7: Performance 📊 PENDING

- [ ] Create `benchmarking/benchmark_pladdrr_vs_parselmouth.R`
- [ ] Benchmark all 18 functions
- [ ] Profile bottlenecks
- [ ] Create `PLADDRR_PERFORMANCE_ANALYSIS.md`
- [ ] Create `PLADDRR_OPTIMIZATION_REQUESTS.md` for developers

## Phase 8: Documentation 📝 PENDING

- [ ] Update CLAUDE.md
- [ ] Update README.md
- [ ] Create `MIGRATION_GUIDE.md` for users
- [ ] Create `PLADDRR_API_REFERENCE.md`
- [ ] Create `PLADDRR_VS_PARSELMOUTH.md`
- [ ] Update `PKGDOWN_FUNCTION_GROUPING.md`
- [ ] Update all function man pages
- [ ] Update NEWS.md

## Phase 9: Polish ✨ PENDING

- [ ] Run `devtools::check()` - zero errors
- [ ] Fix all warnings
- [ ] Run `lintr::lint_package()`
- [ ] Update build configuration
- [ ] Final testing
- [ ] Prepare release notes

## Phase 10: Transition 🚀 PENDING

- [ ] Create GitHub release
- [ ] Update website
- [ ] Announce changes
- [ ] Archive parselmouth code in branch

## Known Issues

### Formant Bug (v4.8.16)
- **Status**: To be verified
- **Issue**: Previous versions (v4.6.4) had 35-55% underestimation
- **Resolution**: Fixed in v4.8.16 per release notes
- **Action**: Test and document if still present

### Performance Concerns
- Monitor for any slower implementations
- Document and report to pladdrr developers
- Keep optimization list

## Success Metrics

- [ ] All 10 existing functions migrated
- [ ] All 8 new functions added
- [ ] Zero test failures
- [ ] Package passes `devtools::check()`
- [ ] All functions documented
- [ ] SSFF output unchanged (emuR compatible)
- [ ] Performance report complete

## Timeline

- **Started**: 2026-02-06
- **Target**: 2026-02-27 (3 weeks)
- **Current Phase**: 3 (Migration)
