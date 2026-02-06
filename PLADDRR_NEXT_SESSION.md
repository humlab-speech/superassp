# Pladdrr Integration - Next Session Quick Start

**Current State**: Branch `pladdrr-integration`, 2/18 functions complete  
**Last Commit**: `daaaa88` - docs: Add session 2 summary

## Immediate Next Task: trk_formantp

### Quick Commands

```bash
# Check status
cd /Users/frkkan96/Documents/src/superassp
git log --oneline -5
git status

# Load in R
R
devtools::load_all()

# Test current functions
test_file <- system.file("samples/sustained/a1.wav", package = "superassp")
trk_pitchp(test_file, toFile = FALSE)  # Should work
```

### Implementation Steps for trk_formantp

1. **Read source**: `../plabench/R_implementations/formant.R` (238 lines)

2. **Create/edit**: `R/ssff_python_pm_pformantb.R`

3. **Key pladdrr calls**:
```r
formant_ptr <- pladdrr::to_formant_direct(
  sound,
  time_step = 0.005,
  max_n_formants = 5,
  maximum_formant = 5500,
  window_length = 0.025,
  pre_emphasis_from = 50
)
formant_obj <- pladdrr::Formant(.xptr = formant_ptr)
df <- formant_obj$as_data_frame()  # Long format: (time, formant_number, freq, bw)
```

4. **Convert to wide format**:
```r
# Use helper from pladdrr_helpers.R
wide_df <- pladdrr_df_to_superassp(df, type = "formant", n_formants = 5)
# Result: (time, fm1, fm2, ..., fm5, bw1, bw2, ..., bw5)
```

5. **Build AsspDataObj**:
- 10 tracks: fm1-fm5 (formant frequencies), bw1-bw5 (bandwidths)
- SSFF format
- Standard superassp interface (toFile, batch, time windowing)

6. **CRITICAL TEST**: Verify formant bug fix
```r
# Known issue in v4.6.4: F1/F2/F3 were 35-55% too low
# v4.8.16 should be fixed - test with known audio
# Compare with praat/other reference implementations
```

### If Formant Bug Persists

- Add warning to function docs
- Create `PLADDRR_OPTIMIZATION_REQUESTS.md`
- Document issue for pladdrr maintainer
- Continue with migration (still better than Python latency)

### After trk_formantp

Move to Batch 2 (summary functions):
1. `lst_avqip` - AVQI (3x faster!)
2. `lst_dsip` - DSI
3. `lst_voice_tremorp` - 18 measures
4. `lst_voice_reportp` - Voice report

### Reference Files

- Template: `R/ssff_python_pm_pitchp.R` (just completed)
- Source: `../plabench/R_implementations/formant.R`
- Helpers: `R/pladdrr_helpers.R`
- Plan: `PLADDRR_IMPLEMENTATION_PLAN.md`
- Status: `PLADDRR_MIGRATION_STATUS.md`

### Timeline

- **Day 3** (today): trk_formantp + formant testing
- **Days 4-5**: Batch 2 (4 summary functions)
- **Target completion**: 2026-02-20 (12 days remaining)

---

**Start here**: Read `../plabench/R_implementations/formant.R` and begin implementation.
