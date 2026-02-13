# Plan: Assess and Update superassp pladdrr API Compatibility

## Context

pladdrr has been updated to v4.8.23. superassp requires `pladdrr >= 4.8.20` and uses 14 functions built on the pladdrr API. This plan assesses API compatibility and identifies needed changes.

## Assessment Summary

**Overall verdict: No breaking changes detected. Minor improvements recommended.**

The superassp codebase uses pladdrr correctly and is fully compatible with v4.8.23. The API surface used by superassp has not changed in breaking ways. Below is the detailed analysis.

---

## 1. API Calls Audited (all OK)

### Tier 1: Object Construction (OK)
- `pladdrr::Sound(file_path)` - Used in 10+ files. Still the primary constructor.
- `pladdrr::Pitch(.xptr = ptr)` - Function factory pattern. Unchanged.
- `pladdrr::PointProcess(.xptr = ptr)` - Unchanged.
- `pladdrr::TextGrid(path)` - Now properly exported (was missing before v4.8.8). No change needed.

### Tier 2: Direct API (OK)
- `to_pitch_cc_direct()` - Full-parameter version (v4.0.2+). Unchanged.
- `to_pitch_ac_direct()` - Full-parameter version. Unchanged.
- `to_formant_direct()` - Unchanged.
- `to_intensity_direct()` - Unchanged.
- `to_harmonicity_direct()` - Unchanged.
- `to_point_process_from_sound_and_pitch()` - Unchanged.

### Tier 4: Ultra API (OK)
- `two_pass_adaptive_pitch()` - Unchanged.
- `get_jitter_shimmer_batch()` - Unchanged.
- `calculate_multiband_hnr_ultra()` - Bug fixes in v4.8.x (CC method now correct). **Results may improve** but API unchanged.
- `get_durations_batch()` - Unchanged.
- `calculate_f0_stats_ultra()` - Unchanged.
- `calculate_minimum_intensity_ultra()` - Bug fix verified. API unchanged.
- `get_voice_quality_ultra()` - Unchanged.
- `sound_concatenate_all()` - Unchanged.
- `extract_voiced_segments_ultra()` - Bug fixes (ZCR accuracy). API unchanged.
- `calculate_cpps_ultra()` - Bug fixes (NA return, defaults). API unchanged.

### R6 Methods (OK)
- `sound$extract_part()` - Unchanged.
- `sound$to_spectrum()` - Unchanged.
- `sound$to_spectrogram()` - Segfault fixed in v4.8.20. API unchanged.
- `sound$to_powercepstrogram()` - Unchanged.
- `sound$to_intensity()` - Unchanged.
- `sound$to_formant_burg()` - Unchanged.
- `sound$to_harmonicity_gne()` - Unchanged.
- `sound$filter_pass_hann_band()` - Unchanged.
- `sound$resample()` - Unchanged.
- `sound$.cpp$duration` / `.cpp$xmin` / `.cpp$xmax` / `.cpp$sampling_frequency` - Property access unchanged.
- `spectrum$to_ltas_1to1()` - Unchanged. (`to_ltas()` without args now also works as alias.)
- `spectrum$get_band_energy()` - Unchanged.
- `spectrum$pass_hann_band()` - Unchanged.
- `spectrum$formula()` - Unchanged.
- `spectrum$to_power_cepstrum()` - **Correct name used** (not the deprecated `to_powercepstrum()`).
- `ltas$get_slope()` - LTAS unit fix in v4.0.4 but superassp passes `"dB"` which is correct.
- `ltas$compute_trend_line()` - Unchanged.
- `ltas$get_maximum()` - Unchanged.
- `formant$get_value_at_time()` - Unchanged.
- `formant$get_bandwidth_at_time()` - Unchanged.
- `pitch$get_value_at_time()` - Unchanged.
- `intensity$get_time_of_maximum()` - Unchanged.
- `intensity$get_value_at_time()` - Unchanged.
- `point_process$get_mean_period()` - Unchanged.
- `point_process$get_stdev_period()` - Unchanged.
- `point_process$to_textgrid_vuv()` - Unchanged.
- `textgrid$get_number_of_intervals()` - Unchanged.
- `textgrid$get_interval_text()` - Unchanged.
- `textgrid$get_interval_start_time()` / `get_interval_end_time()` - Unchanged.
- `power_cepstrogram$get_number_of_frames()` - Unchanged.
- `power_cepstrogram$get_time_from_frame_number()` - Unchanged.
- `power_cepstrogram$to_powercepstrum_slice()` - Unchanged.

### Static Methods (OK)
- `pladdrr::Sound$from_values(data, sampling_rate)` - Used in `lst_voice_tremorp`. Unchanged.

### Internal API (OK)
- `ns$.powercepstrum_get_peak_prominence()` - Used in `trk_cpps` and `lst_vq`. Unchanged.
- `ns$.ltas_get_frequency_of_maximum()` - Used in `lst_pharyngeal`. Unchanged.

---

## 2. Deprecated API Usage Check

| Deprecated Function | Used in superassp? | Status |
|---|---|---|
| `to_powercepstrum()` | **NO** - uses `to_power_cepstrum()` | OK |
| `sound$to_pitch()` | Only in doc comments, not code | OK |
| `sound$to_formant()` | **NO** - uses `to_formant_burg()` | OK |
| `create_sound()` | **NO** - uses `Sound$from_values()` | OK |
| `read_sound()` | **NO** - uses `Sound()` constructor | OK |
| `sound_as_matrix_zerocopy_impl` (removed) | **NO** | OK |

**No deprecated API usage found in superassp code.**

---

## 3. Recommended Improvements (Optional, Non-Breaking)

### 3a. Bump minimum pladdrr version to 4.8.23

**File**: `DESCRIPTION` line 37
**Change**: `pladdrr (>= 4.8.20)` -> `pladdrr (>= 4.8.23)`
**Rationale**: v4.8.23 includes important fixes:
- CPPS default parameter alignment
- NaN/NA input guards on all query methods (prevents crashes)
- Spectrogram segfault fix (v4.8.20)
- XPtr memory corruption fix (v4.8.15)
- `calculate_multiband_hnr_ultra()` CC method fix (correct HNR values)
- `extract_voiced_segments_ultra()` ZCR accuracy fix

### 3b. Use new batch LTAS methods in `lst_pharyngeal` (performance)

**File**: `R/list_pladdrr_pharyngeal.R` function `extract_ltas_peaks()`
**Current**: Loop with `ltas$get_maximum()` + internal `ns$.ltas_get_frequency_of_maximum()` per band
**Available**: `ltas$get_peaks_batch(fmins, fmaxs)` (18x faster per AGENT_GUIDE)
**Impact**: Minor perf improvement for pharyngeal analysis
**Risk**: Low - new API is stable since v4.4.0

### 3c. Use `calculate_cpps_fast()` in `trk_cpps` instead of frame-by-frame loop

**File**: `R/ssff_pladdrr_cpps.R`
**Current**: Manual frame-by-frame loop through PowerCepstrogram slices calling internal API
**Available**: `calculate_cpps_fast(sound)` gives a single CPPS value (not per-frame)
**Note**: `trk_cpps` produces a **time-series** track, while `calculate_cpps_fast()` gives a single summary value. The frame-by-frame approach is **intentionally different** and correct for the `trk_*` use case. No change needed here.

### 3d. Update roxygen doc comments referencing old API names

**File**: `R/pladdrr_helpers.R` lines 33, 131
**Current**: Doc examples use `sound$to_pitch()` (deprecated)
**Should be**: `sound$to_pitch_cc()` or `sound$to_pitch_ac()`
**Impact**: Doc quality only, no functional change

---

## 4. Implementation Plan

### Step 1: Update DESCRIPTION version requirement
- `DESCRIPTION:37` — change `pladdrr (>= 4.8.20)` to `pladdrr (>= 4.8.23)`

### Step 2: Optimize `extract_ltas_peaks()` with batch API
- `R/list_pladdrr_pharyngeal.R:645-666` — replace loop with `ltas$get_peaks_batch()`

### Step 3: Fix doc examples
- `R/pladdrr_helpers.R:33` — `to_pitch()` -> `to_pitch_cc()`
- `R/pladdrr_helpers.R:131` — `to_pitch()` -> `to_pitch_cc()`

### Step 4: Regenerate docs
- `devtools::document()`

---

## 5. Verification

1. `devtools::document()` — regenerate NAMESPACE and man pages
2. `devtools::load_all()` — verify package loads
3. Run pladdrr-related tests: `devtools::test(filter = "pladdrr|vuv|cpps|vq|pharyngeal|voice_report|tremor|avqi|dsi|pitch|formant|intensity|spectral_moments|sauce")`
4. Spot-check a pladdrr function: `trk_cpps("tests/testthat/test_audio.wav", toFile = FALSE)`

---

## 6. Summary

| Category | Finding |
|---|---|
| Breaking changes | **None** |
| Deprecated usage | **None found** |
| Required changes | **0** (all optional) |
| Recommended improvements | **3** (version bump, batch LTAS, doc fix) |
| Risk level | **Very low** |
