# addTrack() format bookkeeping and the `unlockBinding()` check NOTE — 2026-09-20

**Status:** plan. Nothing here is implemented yet.

Two independent defects, both surfaced by the cross-validation work in
`3783415` and by the current `R CMD check` output:

* **A** — `superassp::addTrack()` silently drops the storage format of a *new*
  track, and every one of its 20 call sites in the package works around it.
* **B** — `R/s7_methods.R` calls `unlockBinding()`, which `R CMD check` reports
  as *"possibly unsafe calls"* (the second NOTE in `cran-comments.md`). The
  justification recorded there ("intrinsic to the design") is measurably wrong:
  the call is unnecessary.

## Part A — `addTrack()` drops the format of a new track

### A1 Root cause (measured 2026-09-20)

`wrassp::addTrack()` and ours are the same function except for one assignment
(`deparse()` of the installed wrassp vs `R/assp_dataobj_methods.R:172-175`):

```r
# wrassp
else attr(dobj, "trackFormats") = append(attr(dobj, "trackFormats"), format)
# superassp
else append(attr(dobj, "trackFormats"), format)          # result discarded
```

R's copy-on-modify semantics make the second line a no-op, so the new track's
format never registers. Everything in the package that adds a track compensates
by hand (§A2), and objects built that way cannot be written:

| observation | superassp | wrassp |
|---|---|---|
| `trackFormats` after `addTrack(x, "extra", m, "REAL32")` on a 2-track object | `INT16,INT16` | `INT16,INT16,REAL32` |
| `names()` / `as.data.frame()` / `as_tibble()` | include the new track | include the new track |
| `track_formats()` | reports 2 tracks for 3 | consistent |
| `write_ssff()` | **error**: "Not enough format specifiers for the data tracks." | writes fine |

The bug is ours, not inherited: the comments blaming wrassp
(`R/ssff_pladdrr_intensity.R:189`, `tests/PARSELMOUTH_FIXES_SUMMARY.md:85,94`)
are wrong and must be corrected. Impact today is a loud failure at write time
(`write_ssff`) plus silent metadata inconsistency for in-memory consumers.

### A2 Call-site inventory (20 sites, 8 files)

Every site manages the format vector by hand. Three patterns:

| pattern | sites | post-fix action |
|---|---|---|
| **pre-set before the calls** — `attr(obj,"trackFormats") <- c(...)` / `character(0)` with the final values, then `addTrack()` whose formats are dropped | `R/ssff_cpp_swiftf0.R:193` → `:203,:204`; `R/ssff_cpp_crepe.R:225` → `:236,:238`; `R/ssff_cpp_formantnet.R:188` → `:198,:199`; `R/ssff_cpp_deepformants.R:183` → `:193`; `R/ssff_cpp_covarep_gfmiaif.R:198` → `:216,:224,:232`; `R/ssff.R:206` → `:221`; `R/ssff_cpp_sptk_mfcc.R:193` → `:202` | delete the pre-set line; `addTrack()` supplies the format |
| **append by hand after the call** | `R/ssff_cpp_sptk_mfcc.R:203`; `R/ssff_pladdrr_formant.R:351,:361,:373`; `R/ssff_pladdrr_intensity.R:190`; `R/ssff_pladdrr_pitch.R:189,:356,:586,:764` | delete the manual append (it would duplicate the entry) |
| **replace path** (`deleteExisting = TRUE`) | `tests/testthat/test-ssff-wrassp-golden.R:200` | none — already correct |

The listed line numbers are a map, not a licence: the fix commit must confirm
per site that the pre-set/append belongs to the same object as the `addTrack()`
call (the constructors differ — `list()` + attributes, or an object passed in;
none of the eight files uses `new_asspdataobj()`).

### A3 Fix

1. `R/assp_dataobj_methods.R`: `attr(dobj, 'trackFormats') <- append(attr(dobj, 'trackFormats'), format)`.
2. Delete the compensating lines at the 16 sites mapped in §A2 (7 pre-set
   sites, 9 manual-append sites; the gfmiaif pre-set spans several lines),
   keeping every `addTrack()` call.
3. Correct the misattributed comments; leave `tests/PARSELMOUTH_FIXES_SUMMARY.md`
   in place but annotate it (it is a historical record).
4. `NEWS.md`: bug fix entry (objects built with `addTrack()` now carry complete
   metadata and can be written; the change is a fix, not a behaviour break).

### A4 Verification

| gate | how |
|---|---|
| unit: new-track path | extend `test-ssff-wrassp-golden.R`: `addTrack(x, "extra", m, "REAL32")` → `names`, `track_formats()` and `length(attr(,"trackFormats"))` consistent, write→read round trip equal, values survive |
| unit: parity with the reference | new assertion in `test-ssff-wrassp-interop.R`: the wrassp subprocess applies `wrassp::addTrack()` to a gold fixture and returns it; the parent applies ours to the same input and compares names, `trackFormats` and values |
| invariant for wrappers | assert `length(track_formats(obj)) == length(names(obj))` for the objects produced by the changed wrappers that have tests here: `test-swiftf0.R`, `test-crepe.R`, `test-deepformants.R`, `test-formantnet.R`, `test-gfmiaif.R`. `test-smoke-trk.R:101` already covers `trk_mfcc`; `harmonics()` (internal, `R/ssff.R:204-221`) has no test — add one smoke call so every changed call site is exercised at least once |
| no output drift | for the wrappers above that write files, the written file must keep the same `trackFormats`/values as before the change (their existing round-trip tests cover this; where a fixture exists, compare bytes) |
| suite | `devtools::test()` green, including the pladdrr-gated files when pladdrr is present |

### A5 Risks

* **Pattern-A sites are only safe to simplify if *all* tracks of the object
  arrive through `addTrack()`.** If a site also assigns a track directly
  (`obj[[name]] <- m`), its format must stay in the pre-set; the per-site review
  in step 2 decides. (`ssff_cpp_swiftf0.R:193` is confirmed clean: the object is
  an empty `list()`.)
* **Three sites are pladdrr-gated** (`ssff_pladdrr_formant.R`,
  `ssff_pladdrr_intensity.R`, `ssff_pladdrr_pitch.R`) and cannot be executed in
  this environment. They are covered by the unit-level parity test and by code
  inspection; the residual risk is a typo in a deleted line, which the wide
  `track_formats()` invariant would catch for anyone with pladdrr installed.
* After the fix, previously "metadata-short" objects become consistent — a
  user-visible improvement; nothing that used to work stops working.

## Part B — the `unlockBinding()` NOTE

### B1 What the check reports and what the code does

`R CMD check` → *"checking R code for possible problems ... NOTE / Found the
following possibly unsafe calls: File 'superassp/R/s7_methods.R':
unlockBinding(fn_name, ns)"*. `cran-comments.md` note 2 explains it as
intrinsic to the design: `.onLoad()` → `.setup_s7_methods()` converts every
exported `lst_*`/`trk_*` binding into an S7 generic and writes it back with

```r
unlockBinding(fn_name, ns); assign(fn_name, generic_fn, envir = ns); lockBinding(fn_name, ns)
```

### B2 The unlock/lock pair is unnecessary (measured)

Minimal-package probe (install + `loadNamespace`, and `pkgload::load_all`):

| question | result |
|---|---|
| `bindingIsLocked("f", ns)` during `.onLoad` | **FALSE** |
| `environmentIsLocked(ns)` during `.onLoad` | **FALSE** |
| plain `assign()` during `.onLoad` | works |
| is the replacement effective after loading? | yes — `pkg::f()` returns the replaced value and the binding *is* locked once the namespace is sealed |
| `assign()` in `.onAttach` instead | fails: "cannot change value of locked binding" |

Check-detection probe (one package containing both patterns): `R CMD check`
flags `unlockBinding()` (NOTE) and does **not** flag a plain `assign()` inside
`.onLoad`. `.setup_s7_methods()` has exactly one caller (`R/zzz.R:49`, inside
`.onLoad`), so the unlocked window is the only context it ever runs in.

### B3 Fix

1. Delete the `unlockBinding()`/`lockBinding()` pair; keep the `assign()`.
2. Document the invariant where it now matters — the helper's roxygen already
   says "called during `.onLoad()`"; add a comment that the assignment relies on
   the namespace not yet being sealed, and that `assign()` fails (rather than
   silently no-ops) once it is. A runtime guard is *not* recommended: the
   conversion loop wraps each function in `tryCatch()`, so an error would be
   swallowed as "skipped" instead of failing loudly.
3. `cran-comments.md`: rewrite note 2 — the NOTE is removed, and the "intrinsic
   to the design" justification is replaced by the measurement above.
4. `NEWS.md`: one line under "Notes"/developer notes that the check NOTE is gone
   (no user-facing change).

### B4 Verification

| gate | how |
|---|---|
| the NOTE is gone | `R CMD build` + `R CMD check --no-manual --no-vignettes --no-build-vignettes` on the tarball; "checking R code for possible problems" must be OK and the status must not regress elsewhere |
| no other flagged calls | grep `R/` for `unlockBinding`, `assignInNamespace` (the roxygen `superassp:::` mentions are documentation text, not calls) |
| S7 behaviour unchanged | `test-s7-dispatch.R` (12 assertions), `test-s7-avaudio.R` (22), `test-smoke-trk.R`, `test-smoke-lst.R`, `test-export-policy.R`; manual spot check that `class(trk_rms)` is `S7_generic` with `ext`/`tracks`/`outputType` attributes intact and that both `character` and `AVAudio` inputs dispatch |
| dev workflow | `pkgload::load_all()` path verified by the probe; the suite runs under it |

### B5 Risks

* A future refactor that moves the conversion out of `.onLoad` would make
  `assign()` fail; because of the surrounding `tryCatch()` that shows up as a
  *silent* skip. Mitigation: the comment in `.convert_to_s7_generic`, plus the
  existing dispatch tests, which fail loudly when a function stops being a
  generic. (Optional hardening: assert in `test-s7-dispatch.R` that the known
  converted set is still a set of `S7_generic`s.)
* Removing `lockBinding()` does not leave the binding unlocked: the namespace is
  sealed after `.onLoad`, which locks it (probe).

## Part C — order and combined gates

1. Part B first (two-line change, independent of Part A): edit → suite → check.
2. Part A per file, one commit-sized step per group (grep-verified each site).
3. `NEWS.md` + `cran-comments.md` updates.
4. Full `devtools::test()`, then `R CMD check` on a built tarball for both
   effects (NOTE gone, no new warnings), then `graphify update .`.

## Appendix — reproducing the measurements

```bash
# A: the writer/reference divergence
Rscript -e 'suppressMessages(library(wrassp)); cat(paste(deparse(wrassp::addTrack), collapse="\n"))'
# B: binding state during .onLoad (see /tmp/bindprobe in the session that produced this plan)
R CMD INSTALL --library=<lib> <pkg>      # .onLoad prints bindingIsLocked()/assign() results
# B: what the check flags
R CMD build <probe pkg> && R CMD check --no-manual <tarball>   # unlockBinding -> NOTE, assign -> clean
```

---

# Implementation record (2026-09-20)

Status: **implemented**. Part B and Part A both landed; verification below.

## What changed

| File | Change |
| --- | --- |
| `R/assp_dataobj_methods.R` | `addTrack()`: `attr(dobj, 'trackFormats') <- append(...)` (the new-track branch assigned nothing before) |
| `R/s7_methods.R` | `.convert_to_s7_generic()`: `unlockBinding()`/`lockBinding()` pair removed; comment records the `.onLoad()` invariant |
| 10 wrapper files | 20 compensating lines deleted; each verified to be a site that only extends metadata the fixed `addTrack()` now sets itself |

Pure-deletion check (difflib, HEAD vs working tree): every wrapper file has
`+0` added lines. `ssff.R` 1, `ssff_pladdrr_pitch.R` 4, `ssff_pladdrr_formant.R` 3,
`ssff_pladdrr_intensity.R` 2, `ssff_cpp_sptk_mfcc.R`/`swiftf0`/`crepe`/
`deepformants`/`formantnet` 1 each, `ssff_cpp_covarep_gfmiaif.R` 5.
One indentation artefact introduced by the deletion in
`ssff_pladdrr_formant.R` was found and fixed (closing brace 12 → 8 spaces).

## Evidence

**No output drift** (`/tmp/ssff_impl/drift_check.R` + `gfmiaif_snapshot.R`
against a detached worktree at `HEAD` with the current DLL copied in — R-only
changes, so the worktree reproduces pre-change behaviour exactly):

| Artifact | before vs after |
| --- | --- |
| `trk_gfmiaif(nv = 12)` | md5 `101752f840ccc8f0f01d5efb85c3405a` both, 19 tracks, all `REAL64` |
| `trk_mfcc()` | byte-identical, 22948 bytes, 14 tracks, all `REAL32` |
| `harmonics(n = 3)` | byte-identical, 2556 bytes, 1 track, `INT16` |
| `trk_pitch_cc` / `trk_pitch_ac` / `trk_pitch_shs` / `trk_pitch_spinet` | byte-identical, 1 track, `REAL32` each |
| `trk_intensity` | byte-identical, 1 track, `REAL32` |
| `trk_formant_burg` | byte-identical, 15 tracks, all `REAL32` (15 formats) |

**Tests added/changed**

| Test | Result |
| --- | --- |
| `test-addtrack.R` (new) | 24 assertions: new track registers format, replace in place, single-track branch relaxes the row check but still needs `deleteExisting`, duplicate/mismatched-rows/non-numeric/bad-name/`AsspDataObj` errors, `harmonics()` smoke + round trip |
| `test-ssff-wrassp-golden.R` | +6: new track survives write/read with `c("INT16","INT16","REAL32")` |
| `test-ssff-wrassp-interop.R` | +5: `addTrack()` parity against `wrassp::addTrack()` (names, formats, values) |
| 5 wrapper tests | +1 invariant each: `length(track_formats(result)) == length(names(result))` |

**Coverage gap (honest):** `test-swiftf0.R`, `test-crepe.R`,
`test-deepformants.R`, `test-formantnet.R` first tests `skip_*` in this
environment (ONNX runtime/models), so the invariant lines added there execute
only where those models are present. The four sites were verified by
inspection instead: no track is ever assigned directly (`outDataObj[[<-`, `$<-`)
and every track goes through `addTrack()`.

**Binding reduction (Part B, argument + probe, not a test):** removing the
`unlockBinding()`/`lockBinding()` pair is valid because `.onLoad()` runs before
the namespace is sealed. Kept out of the test suite: any test that reads
`bindingIsLocked()` passes under both variants and so asserts nothing about the
change. The real gate is `R CMD check`, which is the reason the change exists.

## Gates (final)

| Gate | Result |
| --- | --- |
| Full suite (`test_dir`, summary reporter) | **3271 passed, 0 failed, 0 error, 40 skipped** (previous run: 3235/0/0/40 — the +36 are exactly the new assertions minus the four that live in ONNX-gated tests) |
| `R CMD check` on `superassp_3.1.0.tar.gz` (`--no-manual --no-vignettes --no-build-vignettes --no-tests`, tests run separately above) | **Status: OK** — zero NOTE/WARNING/ERROR, including `checking R code for possible problems ... OK`, the section that flagged `unlockBinding()` before |
| `graphify update .` | 13904 nodes, 40623 edges, 673 communities |

## Out of scope

`read_ssff`/`read_track` (the SSFF reading path) is untouched by this change;
the `0 → NA` and one-pass-read work is recorded in
`planning/2026-09-20-ssff-read-performance.md`.
