# superassp — Comprehensive Package Assessment

*Generated 2026-09-06. Read-only audit (4 parallel deep-dive agents); no code changed. File:line references verified against current `master`.*

---

## 1. Runtime power-consumption reduction opportunities

Ordered by likely CPU-cycle/power impact.

1. **Five SPTK pitch/formant wrappers run sequentially with zero parallelization**, despite CLAUDE.md's "auto-parallelize for 2+ files" contract and a working template already in the codebase.
   - Sequential `for` loops, no `parallel::`/`mclapply`/`pbmcapply`: `R/ssff_cpp_sptk_rapt.R:99`, `ssff_cpp_sptk_reaper.R:85`, `ssff_cpp_sptk_swipe.R:84`, `ssff_cpp_sptk_dio.R:74`, `ssff_cpp_sptk_harvest.R:93`.
   - Working counter-example already in-tree: `R/ssff_cpp_estk_pitchmark.R:182-194,423-460` auto-selects `mclapply`/`pbmcapply` (fork, Unix) or `parLapply` (PSOCK, Windows) for >1 files.
   - **Fix**: extract that dispatch block into a shared helper (e.g. `run_parallel_files()` in `R/track_helpers.R`) and adopt it in the 5 SPTK wrappers. These are the most heavily-used pitch trackers — batch jobs over N files currently burn N× wall-clock/idle-core power where cores are available.

2. **Duplicated, un-vectorized autocorrelation + Levinson-Durbin implemented three separate times**, only one of which uses the mandated SIMD path.
   - `src/vat_iaif_lpc.cpp:19-22` (runs once per GCI — potentially thousands of calls/file) and `src/vat_srh_pitch.cpp:25-27` (once per 25ms frame) both hand-roll a scalar `for(k) for(i) r(k)+=s(i)*s(i+k)` loop, ignoring `src/simd_utils.hpp`'s `sasp::simd_dot`/`simd_energy`.
   - `src/srh_variant.cpp` implements the *same* algorithm correctly using SIMD — proving the inconsistency is avoidable, not architectural.
   - A fourth, pure-R copy of the same math exists in `R/list_r_polarity.R:266-293` (`.polarity_lpc()`), run in interpreted R — far slower than any of the C++ copies.
   - **Fix**: one `sasp::simd_autocorr(x, order)` primitive in `simd_utils.hpp` (double precision, `#ifdef RCPPXSIMD_AVAILABLE` guarded), consumed by all three C++ sites; expose a thin internal Rcpp binding for `.polarity_lpc()` to call instead of reimplementing in R.

3. **Quadratic vector growth in `R/list_r_polarity.R:214,248`** — `residuals <- c(residuals, ...)` inside a per-frame loop reallocates/copies the whole accumulated vector every iteration (O(n_frames²)). Preallocate a list/vector and `unlist()` once.

4. **Window-function formulas (Hamming/Hann) hand-written independently in ~11 files** (`list_cpp_covarep_gci.R`, `list_r_polarity.R`, `ssff_cpp_covarep_{creak,env_te,gfmiaif,hmpd,vad_drugman,vq_gci}.R`, `vat_internal_{creak,dsp,voice_quality}.R`). Low per-call cost but wasted repeated cycles plus a fidelity-drift risk (any subtle formula divergence between the 11 copies is itself a correctness bug). Consolidate into one `R/window_helpers.R` (+ a C++ equivalent in `simd_utils.hpp`), verified numerically identical before switching call sites.

5. **`build_media_manifest()` (`R/audio_loader.R:154-169`) is dead code** — designed to do one upfront existence/extension pass across a batch instead of N redundant per-file checks, but called nowhere except its own test. Low individual cost (the expensive ffprobe call is already deduped via `R/cache_media_info.R`), but either wire it in or delete it.

**Confirmed non-issues**: `assp_load_audio_for_dsp`/`av_to_asspDataObj`/`media_info()` do **not** redundantly reload or re-transcode — `media_info()` is properly memoized and shared across call sites.

---

## 2. Runtime performance opportunities (non-power-specific)

Largely the same list as §1 (power and CPU-cycle savings coincide here) — the additional performance-only findings:

6. **`R/ssff_cpp_tvwlp.R:542-552`** — `.tvwlp_build_system()` fills the LP regression matrix with a scalar R triple-nested loop, once per ~200ms analysis frame, despite the file name implying a C++ backend (`src/tvwlp_core.cpp` exists and is used elsewhere in the same pipeline). Vectorize with `embed()`/matrix slicing, or move the fill into the existing C++ file as an Armadillo routine.

7. **Rcpp-boundary check**: no anti-pattern found of per-frame R→C++ calls inside a hot loop; the one candidate (`R/ssff_pladdrr_spectral_moments.R:154-167` looping into pladdrr's R6 Praat objects per-frame) is inherent to the Praat object model, not a fixable boundary issue.

8. **Sibling-algorithm duplication is real but scoped**: multiple pitch trackers (rapt/swipe/reaper/harvest/dio/yin/pyin/srh/ac/cc/ksv/mhs/vat/snack/pda/crepe) intentionally offer different algorithms (not redundant), but do redundantly reimplement common DSP primitives (windowing, autocorrelation — see §1.2/1.4) rather than sharing helpers.

---

## 3. Standards compliance, code quality, maintainability

Ordered by impact.

1. **Dead validation helpers vs. 20+ file duplicated inline boilerplate.** `R/validation_helpers.R` defines `validate_file_paths()`, `validate_jstf_output_params()`, `validate_time_range()` explicitly as the shared parameter-checking layer — **zero callers** anywhere in `R/`. Meanwhile the same error strings are duplicated by hand across the `ssff_*.R` family: "No input files specified" (23 files), "Some files do not exist" (21 files), "Error processing {.file...}" (25 files), "Successfully processed" (14 files). Cheapest, highest-value single fix: start calling the existing helpers instead of the inline copies.

2. **Unused Imports in `DESCRIPTION`**: `tidyr`, `tidyselect`, `R.utils` have zero call sites anywhere in `R/`. Remove them (reduces install footprint/attack surface).

3. **CLAUDE.md ↔ DESCRIPTION drift**: CLAUDE.md's dependency list places `R.matlab` under Required/Imports; `DESCRIPTION` actually lists it under Suggests. Fix the doc.

4. **Undocumented naming-scheme origin tag**: 7 files use `voiceanalysis` as an implementation-origin token (`ssff_voiceanalysis_{creak,gci,iaif,mdq,peakslope,pitch}.R`, `list_voiceanalysis_{lf,vq}.R`), which isn't in CLAUDE.md's enumerated origin list (`c_assp, cpp, cpp_sptk, cpp_estk, cpp_snack, cpp_covarep, cpp_opensmile, pladdrr, r`). These front the pure-R `vat_internal_*` pipeline — likely should be `r_vat` or the scheme doc should be extended to name `voiceanalysis` explicitly.

5. **14 of 62 exported `trk_*` functions have no dedicated test**: `trk_afdiff`, `trk_affilter`, `trk_css_spectrum`, `trk_formant_snack`, `trk_ksvfo`, `trk_lps_spectrum`, `trk_pitch_ac`, `trk_pitch_crepe`, `trk_pitch_mhs`, `trk_pitch_pda`, `trk_pitch_shs`, `trk_pitch_snack`, `trk_pitch_spinet`, `trk_pitchmark_estk`, `trk_praatsauce`. All `lst_*` exports have coverage.

6. **6 legacy `test_<name>.R` (underscore) files** alongside ~65 `test-<name>.R` (hyphen) files — both run under testthat's discovery, but the underscore files (`test_aaa_initDemoDatabase.R`, `test_parallel_processing.R`, `test_praat.R`, `test_praat_slicefunctions.R`, `test_ssff.R`, `test_wrassp.R`) are stale wrassp-era stubs, several mostly commented out. Candidates for removal (overlaps with §6 below).

7. **Minor naming inconsistencies**: `lst_GeMAPS`/`lst_eGeMAPS`/`lst_ComParE_2016` use camelCase segments (only camelCase cluster among 19 `lst_*` exports — likely intentional to preserve official openSMILE names, worth an explicit documented exception). Separately, four internal (correctly unexported) Layer-2 functions reuse the `lst_` prefix (`lst_GeMAPS_cpp`, `lst_eGeMAPS_cpp`, `lst_ComParE_2016_cpp`, `lst_emobase_cpp`) instead of the documented `<algo>_cpp` pattern — confusing next to the real exported `lst_*` functions, though not a policy violation.

8. **Fully compliant**: Export Policy (all 114 exports match policy regex, `test-export-policy.R` passes conceptually) and error-handling standard (zero bare `warning()`/`stop()` outside one doc-comment example) — both are clean.

---

## 4. Vignettes (highest documentation-effort/value ratio)

Only one vignette exists: `vignettes/getting_started.Rmd` (137 lines) — covers `trk_*`/`lst_*` conventions, `AsspDataObj`/`JsonTrackObj` shape, and a pitch/formant comparison table. Given 114 exported functions, this drastically under-represents the package.

**Recommended new vignettes, highest value first:**

1. **Voice quality and creak** — `trk_covarep_iaif`/`trk_gfmiaif` → `lst_vq`/`lst_avqi`/`lst_dsi` pipeline; when to pick the VAT- vs COVAREP-backed variant.
2. **Extracting standardized feature sets with openSMILE** — `lst_GeMAPS`/`lst_eGeMAPS`/`lst_emobase`/`lst_ComParE_2016`, output-column meaning, batch extraction across a corpus.
3. **SSFF and JSTF: reading, writing, round-tripping** — `read_ssff`/`write_ssff` vs `read_jstf`/`write_jstf`, JSTF slice-manipulation helpers, and loading into `emuR` (a headline interoperability claim in both the README and the existing vignette that is never actually demonstrated anywhere).
4. **Unit conversions for psychoacoustics** — one cookbook page for all 14 `ucnv_*` functions, currently zero vignette mentions.
5. **Neural and Praat-backed trackers: setup and caveats** — ONNX auto-download behavior (`trk_pitch_crepe`, `trk_formant_deepformants`, `trk_formant_formantnet`) and the Praat/pladdrr dependency scope (`trk_formant_burg`, `trk_praatsauce`, `lst_pharyngeal`, `lst_voice_tremor`, `lst_voice_report`) — this is also where the README's misleading blanket Praat requirement (§5) should actually be clarified for users.
6. *(nice-to-have)* **From superassp to emuR** — demonstrate the interoperability claim end-to-end.

---

## 5. Documentation quality (roxygen + README)

**Roxygen** (12 functions sampled across families): the large majority are good — usage-first descriptions, clear `@return` schemas, runnable `@examples` (`trk_pitch_rapt`, `trk_formant_burg`, `lst_voice_report`, `read_ssff`/`write_ssff`, `trk_covarep_iaif`, `trk_pitchmark_reaper`, `ucnv_hz_to_bark`, `trk_pitch_crepe`). Two clear 90/10-rule violations:
- `R/list_cpp_opensmile_gemaps.R:8-41` — `@details` spends ~34 lines enumerating every LLD/functional/statistic before `@param` even starts; should be trimmed to a short summary + pointer to a reference vignette.
- `R/assp_dataobj_methods.R:1-13` — `read.AsspDataObj` (legacy, unexported) has a thin `@return` ("list object containing file data" — doesn't even name it as `AsspDataObj`) and **no `@examples`**; either flesh it out or mark `@keywords internal @noRd` since it's already unexported.
- Minor: `R/ssff_voiceanalysis_creak.R:8-13` front-loads "Bit-faithful Rcpp port of the Kane-Drugman + Ishi 36-feature pipeline..." ahead of usage — borderline but short enough to be low-priority.

**README.md** (179 lines, no `README.Rmd` — risk of example drift since nothing re-executes it):
- **Accuracy bug**: line 46 states "The package requires the Praat program to be installed in the user's PATH" as a blanket requirement, contradicting both the vignette ("No external Praat or wrassp installation is required") and reality (Praat/pladdrr is only needed for a specific function subset). This is a newcomer-facing correctness bug — scope the wording.
- Typos: "succuessor", "inporporates", "audil feature extractor", "assocaiated".
- Quick Start only demonstrates pitch/formant tracking — no `lst_*`, `ucnv_*`, or JSTF example, under-representing breadth relative to the 114-function export surface.

---

## 6. Superseded functionality — recommend removing

Ordered by value; all confirmed by direct usage-grep before flagging, nothing here is speculative.

1. **15 tracked dead scripts directly under `tests/`** (not `tests/testthat/`) — these run as standalone R CMD check scripts and reference things that no longer exist: `reticulate::*`, `praat_intensity_opt()`/`praat_formant_burg_opt()` (Python removed April 2026, never existed in current `R/`), `av_load_for_python()` (explicitly noted as removed in `R/wav_helpers.R:109`), `lst_phonet()`/`phonet_available()` (undefined), old wrassp-era `rmsana()` calls, and `source("../R/av_helpers.R")` (won't resolve from an installed package). Files: `debug_file_writing.R`, `debug_formant_burg.R`, `debug_intensity.R`, `test_python_dsp_code_verification.R`, `test_python_memory.R`, `test_python_memory_based_dsp.R`, `test_opensmile_code_verification.R`, `test_opensmile_memory.R`, `test_ftrack_tvwlp.R`, `test_phonet_integration.R`, `test_load_and_process.R`, `test_memory_dsp.R`, `test_av_debug.R`, `test_av_integration.R`, `test_non_native_formats.R`. **`git rm` all 15** — single highest-value cleanup in this whole assessment; they will fail if `R CMD check --as-cran` ever actually executes them.

2. **4 orphaned submodules contributing nothing to the compiled library**: `src/ESTK`, `src/tcl-snack`, `src/Yin-Pitch-Tracking`, `src/pyin` — none appear in `src/Makevars`' include paths or source lists; the equivalent functionality (`estk_pda.cpp`/`estk_pitchmark.cpp`, YIN/pYIN in `src/yin_wrapper.cpp`, Snack-derived code from `SPTK/third_party/Snack/`) is actually implemented as standalone code or vendored inside the SPTK submodule instead. Same category of waste already fixed once in commit `591460f` for Parselmouth/praat.github.io/swift-f0. `inst/SUBMODULES.md` incorrectly still lists all 4 as "must be pinned." **Investigate `src/ESTK` specifically before removing** — see next item.

3. **Broken ESTK-binary fallback + already-deprecated wrapper**: `R/ssff_cpp_estk_pitchmark.R` has a fallback to an external `pitchmark` binary via a hardcoded developer-machine path (`/Users/frkkan96/Documents/src/superassp/src/ESTK/bin/pitchmark`), which is unreachable on any other machine since ESTK isn't compiled. The wrapper itself is already `@note **DEPRECATED**. Use protoscribe::draft_pitchmark() instead`. Delete the dead fallback branch now; drop the wrapper entirely once the protoscribe migration completes — this also confirms ESTK has no live consumer, clearing the way for #2's submodule removal.

4. **`inst/onnx/swift-f0/` has been reintroduced as an untracked plain directory** — this is the *exact* thing commit `591460f` removed as an orphan submodule, now back on disk as a full unmanaged Python package (`swift_f0/core.py`, `model.onnx`, `pyproject.toml`, `__pycache__`), directly undoing both the orphan-submodule cleanup and the Python-removal migration. A stale `[submodule "inst/onnx/swift-f0"]` block also still lingers in `.git/config`. **Delete the directory, clear the stale `.git/config` entry, do not commit it.**

5. **`tests/testthat/_problems/`** (87 untracked files, `test-*-phaseN-N.R`-style auto-extracted debug snippets, headed "# Extracted from test-X.R:LINE") — scratch tooling output, not source. **Add to `.gitignore` or delete.**

6. **`graphify-out/`** (42MB, untracked analysis-tool output) — **add to `.gitignore`**, don't commit.

7. **Not actually broken** (contradicts an initial hypothesis worth recording so it isn't re-investigated): `src/opensmile` shows as untracked/modified only because of in-tree build byproducts (`src/opensmile/build_r/`, a stray `.o`) — the submodule link itself is correctly pinned as of commit `ba9cc6f`. No action beyond confirming `.gitignore` covers the nested `.o`.

8. **Deprecated accessor aliases** (`rate`, `numRecs`, `dur`, `startTime`, `tracks`): zero internal call sites remain (not exercised by the package's own tests/examples), but NEWS.md sets no target removal version. Not overdue by the stated policy, but worth pinning removal to a named version (e.g. 3.0.0) now that internal usage is fully gone.

9. **`tests/testthat/test-new-asspdataobj.R`** — checked and it is **not** dead/duplicate work; it's a legitimate new faithfulness test for `new_asspdataobj()`, already consumed by several existing wrappers. Keep/commit as in-progress work (matches memory note on the deferred "Task 6" constructor work).

---

## 7. pkgdown site improvements

1. **`trk_ksvfo` is exported but missing from `_pkgdown.yml`** entirely — will render as an "Undocumented/missing" entry. Add it (likely alongside `trk_pitch_ksv` in the Pitch & F0 section, since it's the deprecated alias per §6.8's sibling finding).
2. **Legacy unexported functions surfaced in a primary section**: `_pkgdown.yml:143-144` lists `read.AsspDataObj`/`write.AsspDataObj` under "I/O — Audio & SSFF" even though neither is in `NAMESPACE`'s exports — move both to the existing "Legacy / Internal Reference" section (currently only holds `harmonics`) so the primary I/O section only surfaces the current `read_ssff`/`write_ssff`.
3. **No `home:`/`navbar:` customization** — the file is just `url:` + `template: bootstrap: 5` + `reference:`. With only one vignette today, pkgdown's auto-generated single-item "Articles" dropdown is easy to miss; once §4's new vignettes land, add an explicit `navbar: articles:` menu and a `home:` block (tagline/links) so guides get first-class nav treatment instead of being buried.
4. Section groupings themselves (15 titled sections mapping cleanly onto the CLAUDE.md function-family taxonomy) are otherwise coherent — no restructuring needed beyond items 1–3.

---

## Suggested execution order

Given effort/value:
1. `git rm` the 15 dead `tests/*.R` scripts (§6.1) — zero risk, immediate hygiene win.
2. Gitignore `graphify-out/`, `tests/testthat/_problems/`; delete reintroduced `inst/onnx/swift-f0/` (§6.4-6.6).
3. Wire the 5 SPTK wrappers into the existing parallel-dispatch pattern (§1.1) — highest perf/power win for the effort.
4. Route `ssff_*.R` error paths through the already-written `R/validation_helpers.R` (§3.1) — highest maintainability win for the effort.
5. Add the SIMD autocorrelation primitive and consolidate the 3 duplicate implementations (§1.2).
6. Start the two highest-value vignettes (voice quality, SSFF/JSTF I/O) (§4).
7. Fix README's Praat-requirement overstatement + typos, trim `lst_GeMAPS`'s `@details` (§5).
8. Remove the 4 orphaned submodules once ESTK's live-consumer status is confirmed void (§6.2-6.3).
9. pkgdown fixes (§7) — cheap, do alongside the vignette work.
