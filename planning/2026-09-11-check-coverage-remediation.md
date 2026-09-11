# superassp — `R CMD check --as-cran` + CI coverage remediation plan

*Generated 2026-09-11. Every finding below is backed by an observed artifact (local check log, local tarball, or GitHub Actions run log); sources are named inline.*

## Execution status (2026-09-11, later the same day)

| Tier | Commit | State |
|---|---|---|
| T1 coverage diagnostics | `9b32e4f`, `280acde` | done — second iteration needed `clean = FALSE`, see below |
| T2 `.Rbuildignore` | `7f897a3` | done — tarball 88 MB → 17 MB, 87 object files → 0 |
| T3 DESCRIPTION deps | `f358932` | done — pladdrr optional, verified with the installed package hidden |
| T4 bibliography DOIs | `5ea3b58` | done — all three dead DOIs fixed/removed |
| T5 coverage failure | `8fa3663` | done — root cause: cli line-wrapping; confirmed green by run [34613963214](https://github.com/humlab-speech/superassp/actions/runs/34613963214) |
| T6 local `--as-cran` loop | `538eb25` | done — local check went from `1 ERROR` to `0 ERROR`, 7 WARNINGs (all allowlisted on CI), 5 NOTEs |
| (follow-up) coverage job deps | `9cd31bc` | done — `covr::to_cobertura()` needs `xml2`, which the workflow did not install |
| (residue) Rd line widths | `29f6aa5`, `3307b6d` | done — 76 `\usage` blocks + 8 `\examples` lines wrapped; NOTE cleared |
| (residue) NEWS titles | — | no change — unversioned historical section, retitling would fabricate a release |

**Coverage job status: green.** Run 34613963214 finished `success` with every step passing (`Test coverage` success, `codecov/codecov-action` success, the failure-artifact step correctly skipped) — the first green `test-coverage` run in the workflow's recorded history (29 prior runs: 22 failures, 7 cancelled, 0 successes). The diagnostic scratch branch `ci/coverage-diagnostics` has been deleted again (remote and local).

**R-CMD-check status: green on all four legs.** Run 34627852582 (`ubuntu-latest release`, `ubuntu-latest devel`, `macos-latest release`, `windows-latest release`) finished `success`. The workflow has no `workflow_dispatch`, so it was run from a throwaway `ci/check-verify` branch carrying an extra `workflow_dispatch:` trigger — that commit and the branch were deleted after the run; `master` never contained them. Together with the local `--as-cran` run (which adds the PDF manual and drops the warning allowlist), this covers every platform the project builds on.

Post-execution local `R CMD check --as-cran` (R 4.6.1, full manual): `Status: 7 WARNINGs, 4 NOTEs`, `checking PDF version of manual ... OK`. The remaining NOTEs are the accepted set: the `Remotes` field and `pladdrr` in `Suggests` (inherent to a GitHub-only optional dependency), two unversioned historical `NEWS.md` headings, a missing recent HTML Tidy, plus a single intentional `unlockBinding()` in `R/s7_methods.R`.

### Residue tier (executed after the six main tiers)

* **Rd line widths NOTE: cleared.** 76 `\usage` blocks exceeded the 90-column limit (the longest was 553 columns) because explicit `@usage` is used throughout so the docs show the real parameter list rather than the runtime S7 generic's `(listOfFiles, ...)` signature. All 76 were re-wrapped in roxygen's own generated-usage style (one argument per line) — `29f6aa5`. Verified by snapshotting all 82 `@usage` blocks before and after: byte-identical modulo whitespace. A follow-up (`3307b6d`) wrapped the eight 132-column `path2wav <- list.files(...)` example lines that the same check flags at the 100-column limit. `tools::checkRd()` on all 280 Rd files now reports zero over-wide `\usage` and zero over-wide `\examples` lines.
* **NEWS section titles: left as-is, deliberately.** The `# superassp (development) — Consistency Refactor` section is unique content (not a duplicate of any versioned section — none of its entries appear under 2.0.0) describing an unreleased API-tightening pass, and its version cannot be established from the file or from `git log`. Retitling it would fabricate release history, so the NOTE stands; R's complaint is informational.
* Remaining NOTEs after this tier: 4 (incoming feasibility, NEWS titles, `unlockBinding`, HTML Tidy).

T5 root cause (named by CI run [34597843212](https://github.com/humlab-speech/superassp/actions/runs/34597843212) once the diagnostics landed): `tests/testthat/test-edge-cases.R:59` asserted `regexp = "not a directory"` against a cli-formatted error. cli wraps that message to the session width, and in covr's test process the wrap falls inside the phrase — the log shows `The path '/tmp/RtmpcinUtK/file2c077addc07f' exists but is not a\ndirectory.` — so the assertion failed under covr while passing under `R CMD check` (whose temp path lays out differently). It was never a DSP or logic failure; it is the only assertion in the suite that a wrapped cli message could break. Fixed by matching `not a\\s+directory` (`8fa3663`).

Two diagnostics iterations were needed: the first (`9b32e4f`) found no `Rout.fail` because `covr::package_coverage()` unlinks its temp install dir on exit; `280acde` passes `clean = FALSE` and reports when the file is absent.

**Sources**
- `build.log` (repo root, `R CMD build` output, 2026-09-11 08:27)
- `superassp.Rcheck/00check.log` (local `R CMD check --as-cran`, R 4.6.1 / aarch64-apple-darwin23)
- `superassp_2.9.5.tar.gz` (local build, 88 MB)
- GitHub Actions run [34322427814](https://github.com/humlab-speech/superassp/actions/runs/34322427814) (`test-coverage`, FAIL) and run 34322427728 (`R-CMD-check`, PASS, 7 WARNINGs / 6 NOTEs — all on ubuntu/macOS/windows legs)
- Local reproduction runs (this session): full `testthat` suite (`FAIL 0`), `NOT_CRAN=true` re-run of the 7 `skip_on_cran` files (`FAIL 0`), and an instrumented `covr::package_coverage()` run (`FAIL 0 | WARN 8 | SKIP 13 | PASS 3276`)

---

## Findings

### F1 — `R CMD check --as-cran` ERROR: `Packages suggested but not available` (local only)

`00check.log`: `ERROR … Packages suggested but not available: 'pbapply', 'pbmcapply', 'geometry', 'gsignal', 'wavelets', 'tuneR'`, followed by `DONE / Status: 1 ERROR, 1 WARNING`.

All six are **current CRAN packages** (verified against `crandb.r-pkg.org`: `200` for each) — they are simply not installed in this machine's library. The check aborted at that ERROR, which is why `00check.log` never reaches the NOTES/WARNINGs that CI does see. CI's own check reports `* checking package dependencies ... OK`.

*Conclusion:* not a package defect; a local-library gap. Two sub-issues must still be handled:
- `gsignal` is in `Suggests` but **referenced nowhere** in `R/`, `src/`, `vignettes/`, or `tests/` (grep for `gsignal` returns no code hits) — dead entry.
- Until the missing packages are installed (or `_R_CHECK_FORCE_SUGGESTS_=false` is set), every local `--as-cran` run aborts before the interesting checks.

### F2 — as-cran NOTE: the tarball ships 87 compiled object files

`R CMD check` on CI: `* checking if this is a source package ... NOTE / Found the following apparent object files/libraries: src/assp/acf.o … src/world/stonemask.o`.

*Verified locally:* `superassp_2.9.5.tar.gz` contains exactly those files — `87` entries matching `\.(o|so|a)$` (the CI list is a truncation of the same set). The `.o` mtimes inside the tarball are `Sep 8 10:31–10:34`, i.e. leftovers of an earlier in-tree compile.

**Root cause (R, not superassp):** `R CMD build` runs `tools:::.build_packages`' `temp_install_pkg(pkgdir, libdir)`, which does `R CMD INSTALL` **on the source directory** so that help pages/vignettes can be processed. That compiles the package in place — `src/assp/*.o`, `src/SPTK/**/*.o`, `src/world/*.o`, `src/tandem/tandem_64/*.o`, `src/opensmile/progsrc/smileapi/*.o`. The subsequent "* cleaning src" step (`tools:::.build_packages` → `cleanup_pkg`) runs `make -f …/clean.mk … clean` and then `unlink(Sys.glob(c("*.o","*.so","*.dylib","*.mod")))` — a **non-recursive** glob in `src/`, so only top-level `src/*.o` are removed. Everything below `src/` survives into the tarball. This is why a pristine CI checkout still shows the NOTE: the build itself creates them.

`tools:::inRbuildignore("src/assp/acf.o", ".")` is currently `FALSE` — the exclusion has to be added by us.

**Local 92 MB tarball:** the object files are 8.1 MB of it. The rest is dominated by `inst/onnx/crepe/full.onnx` (84.9 MB) + `tiny.onnx` (1.9 MB) — files that are **`.gitignore`d** (`.gitignore:45 inst/onnx/crepe/*.onnx`, not in git, not on CI) and **not read by any code path**: `trk_pitch_crepe` downloads its model to `tools::R_user_dir("superassp","cache")/onnx/crepe` via `.hf_get_cached_model()` (`R/ssff_cpp_crepe.R:109`, `R/hf_model_cache.R`). They are a stale local copy of the HF cache, and `R CMD build` (which ignores `.gitignore`) sweeps them in because `.Rbuildignore` says nothing about them. Only `inst/onnx/formantnet/normstats.txt` is tracked and used (`system.file("onnx","formantnet","normstats.txt")`).

### F3 — as-cran NOTE ×2: hidden files and non-standard top-level files

From the CI check (the tarball, not the repo):
```
Found the following hidden files and directories:
  .lintr  src/pyin/.git  src/SPTK/.clang-format  src/SPTK/.git
  src/SPTK/doc/_static/.gitkeep  src/tandem/.git  src/Yin-Pitch-Tracking/.git
  .opencode  src/SPTK/.github

Non-standard files/directories found at top level:
  'AGENTS.md' 'CLAUDE.md' 'codecov.yml' 'planning'
```
The submodule `.git` directories are files (gitlink pointers) present in the tarball; `build.log` (the artifact this plan was written from) is also shipped, though R did not name it.

### F4 — CRAN incoming feasibility WARNING (the only WARNING in `00check.log`)

```
Unknown, possibly misspelled, fields in DESCRIPTION: 'Remotes'
Strong dependencies not in the CRAN or BioC software repositories: pladdrr
Found the following (possibly) invalid DOIs:
  DOI: 10.1109/TASL.2012.2188377  From: man/trk_covarep_vad_drugman.Rd  Status: 404
Size of tarball: 92016058 bytes
```
- **`Remotes`** is not a CRAN field. Its first entry, `github::curso-r/torchaudio`, is **vestigial**: `torchaudio` appears nowhere else in the package (`DESCRIPTION` + a 2023 `NEWS.md` line only); it survived the Python purge.
- **`pladdrr` (≥ 4.8.34) is in `Imports`** but is a GitHub-only package (`github::humlab-speech/pladdrr`). A hard dependency on a non-CRAN package is a submission blocker. Every one of the 92 `pladdrr::` call sites is already behind `pladdrr_available()` / `requireNamespace("pladdrr")` guards (`R/pladdrr_helpers.R`, `R/ssff_pladdrr_*.R`, `R/vat_internal_pitch.R`, `R/list_pladdrr_*.R`) and `NAMESPACE` imports nothing from it — so the dependency is *de facto* optional and the guard is already implemented.
- **DOI `10.1109/TASL.2012.2188377`** resolves `404` at `doi.org` and is absent from Crossref ("Resource not found"). The whole bib entry is suspect: the title *"A comparative study of different feature sets for acoustic voice quality assessment"* has no Crossref record either. The citing function (`R/ssff_cpp_covarep_vad_drugman.R:5`) describes Drugman's multi-branch VAD, whose actual paper is *Voice Activity Detection: Merging Source and Filter-based Information*, IEEE SPL 23(2), 2016, `10.1109/LSP.2015.2495219` (Drugman, Stylianou, Kida, Akamine — verified in Crossref).
- **Tarball size 92 MB** is the F2/F3 payload; CRAN submissions above ~5 MB attract reviewer pushback.
- **CI never sees any of this:** the `R-CMD-check` legs run with `--no-manual --as-cran` but produce no `* checking CRAN incoming feasibility` line (the local log has it between "package encoding" and "namespace information"), and no `_R_CHECK_CRAN_INCOMING*` variable in the runner environment. The CI check reports `Status: 7 WARNINGs, 6 NOTEs` and all 7 WARNINGs are on the workflow's allowlist. F4 is therefore only observable in a real local `--as-cran` run — or on CRAN.

Bulk validation of all **58** DOIs in `inst/REFERENCES.bib` against `doi.org` found **3** dead ones:

| key | current DOI | status | resolution |
|---|---|---|---|
| `Drugman2012VAD` (line 1337) | `10.1109/TASL.2012.2188377` | 404, no Crossref record, **cited** in `R/ssff_cpp_covarep_vad_drugman.R` | replace entry with `10.1109/LSP.2015.2495219` (Drugman/Stylianou/Kida/Akamine, IEEE SPL 2016) |
| `Ishi2008` (line 767) | `10.1109/TASL.2007.907343` | 404 | typo — correct DOI is `10.1109/TASL.2007.910791` (same title/journal/vol 16/pp 47–56/2008, verified in Crossref). Uncited, so as-cran never saw it. |
| `White2022` (line 763) | `10.1515/phon-2022-2011` | 404, no Crossref record found | uncited; verify against the publisher or delete the entry |

### F5 — `R CMD build` WARNING: `Meta/` in the source tree

`build.log:21`: `WARNING: Removing directory 'superassp/Meta' which should only occur in an installed package`. `Meta/vignette.rds` is a leftover of an in-tree `R CMD INSTALL`. `.gitignore` has `/Meta/`, but `.Rbuildignore` does not — so anyone who has ever run `R CMD INSTALL .` in place produces this warning on every build.

### F6 — CI `test-coverage` job fails, and hides why

Run 34322427814, step "Test coverage":
```
Error: Failure in `/tmp/Rtmpioj9Dh/R_LIBS.../superassp/superassp-tests/testthat.Rout.fail`
[ FAIL 1 | WARN 8 | SKIP 30 | PASS 3200 ]
```
Established facts:
- The same commit passes `R CMD check` (which *also* runs the suite, with `NOT_CRAN=true`, through `tools::testInstalledPackage`) on ubuntu-release, ubuntu-devel, macOS and windows.
- Locally (macOS, R 4.6.1, same version as CI) the whole suite passes: `FAIL 0`; and re-running the 7 `skip_on_cran` files with `NOT_CRAN=true` also passes (`FAIL 0`).
- The coverage job runs `covr::package_coverage()`, not `R CMD check`. covr installs with `--example --install-tests --with-keep.source --keep.parse.data --no-staged-install` and `+=-O0 --coverage` makevars, sets `R_COVR=true`, then runs the tests via `tools::testInstalledPackage`.
- So the failing assertion is a **covr-instrumentation-specific** failure. Linux/OpenMP is ruled out as the mechanism: there is no `#pragma omp` anywhere under `src/` (own or vendored sources), so thread-scheduling nondeterminism cannot explain it.
- `covr` prints only the last few lines of `Rout.fail`, so the failing test's name never reaches the CI log — and covr deletes its temp install dir on exit, so the file is unrecoverable after the run. This is why the last three coverage runs (34322427814, 34261846972, 34254224970) are all opaque.
- **This is not a fresh regression.** `gh run list --workflow=test-coverage` returns 29 runs going back to 2026-08-20: 22 failures, 7 cancelled, **zero successes**. The job has never been green in its recorded history.

**Local reproduction (completed).** `covr::package_coverage()` with `options(covr.flags = list())` on this machine (macOS, R 4.6.1, testthat 3.3.2, covr 3.6.5 — all identical to CI) gives:

| | FAIL | WARN | SKIP | PASS | total |
|---|---|---|---|---|---|
| CI run 34322427814 | 1 | 8 | 30 | 3200 | 3231 |
| local covr (instrumented) | **0** | 8 | 13 | 3276 | **3289** |

The WARN count matches exactly, and the local 13 skips are all accounted for (10 × `{voiceanalysis}` not installed, 1 × OpenMP placeholder, 1 × reaper no-voiced-region, 1 × `huggingfaceR`-installed path). But CI skips **17 more** tests and runs **58 fewer** tests overall. Both workflows install `dependencies: all` (r-lib default) and CI's pak lockfile lists every `Suggests` (`huggingfaceR`, `wrassp`, `pladdrr`, …), so this is not a plainly missing-package difference.

Ruled out so far: plain logic error (suite is green locally and green in `R CMD check` on all four CI OS legs); thread nondeterminism (no `#pragma omp` under `src/`); testthat/covr version drift; the untracked `tests/testthat/_problems/` snippets (never collected — confirmed absent from the instrumented `Rout`).

Not yet reproducible locally: covr's **native** instrumentation (`+=-O0 --coverage` on C/C++). Attempting it on this Mac fails at load — `symbol not found in flat namespace '_llvm_gcda_emit_arcs'` — because the bundled sources are built with a mixed clang/libomp toolchain. That flag change is the one remaining uncontrolled variable between the local run and CI, and it is the prime suspect (it changes float code generation for every compiled DSP kernel).

### F7 — Non-fatal check NOTEs worth clearing (CI, ubuntu release)

- `checking R code for possible problems … NOTE`: `unlockBinding(fn_name, ns)` flagged in `R/s7_methods.R`; `.creak_extract_features: no visible global function definition for 'sd'`.
- `checking Rd files … NOTE`: escaped LaTeX specials `\_` in `trk_formant_deepformants.Rd:42`, `trk_formant_formantnet.Rd:42,44`, `trk_pitch_crepe.Rd:64,66`, `trk_pitch_swiftf0.Rd:44,46`.
- `checking Rd line widths … NOTE`: `\usage` lines > 90 chars in `lst_ComParE_2016.Rd`, `lst_GeMAPS.Rd`, `lst_avqi.Rd`, ….
- Compiler warning during the coverage build: `dataobj.c:818: format '%d' expects 'int', argument has type 'long int'` (`dop->numRecords`).

---

## Plan

Each task is self-contained and leaves the package working. Order is by leverage: the two CI-diagnostics/infrastructure tasks first (they unblock the root-cause work), then the packaging fixes, then the coverage failure.

**Global constraints** (from `CLAUDE.md` / existing workflows): no change to DSP numeric output; `devtools::document()` after roxygen edits; the `R-CMD-check.yaml` WARNING allowlist may only be extended with a written justification; the coverage failure must be fixed in the test/tooling, not by deleting the test.

---

### Task 1 — Make CI coverage failures self-diagnosing

**Why:** F6. Currently every failure costs a 35-minute CI run and yields a single unnamed assertion.

**Files:** `.github/workflows/test-coverage.yaml`

- [ ] **Step 1.1 — keep the failure text.** Wrap the coverage call so that on error the `testthat.Rout.fail` tail is dumped into the log *before* covr cleans up:

```yaml
      - name: Test coverage
        run: |
          run <- function() {
            cov <- covr::package_coverage(quiet = FALSE)
            covr::to_cobertura(cov)
          }
          tryCatch(run(), error = function(e) {
            fails <- list.files(tempdir(), pattern = "testthat\\.Rout\\.fail$",
                                recursive = TRUE, full.names = TRUE)
            for (f in fails) {
              cat("\n===== ", f, " =====\n", sep = "")
              cat(readLines(f, warn = FALSE), sep = "\n")
            }
            stop(e)
          })
        shell: Rscript {0}
```
  (covr's install step also needs `quiet = FALSE` or an equivalent: the first local reproduction died with `Package installation did not succeed` and no compiler output.)

- [ ] **Step 1.2 — keep the artifact too** (belt and braces): add `if: failure()` + `actions/upload-artifact@v4` for `**/testthat.Rout.fail` so future failures stay inspectable after the runner is gone.

**Verification:** dispatch `test-coverage` on a scratch branch (`gh workflow run test-coverage.yaml --ref <branch>`); confirm that when the suite fails the log contains the `Failure ('test-….R:NN:NN')` block naming the test. If the suite passes, confirm the added code path is inert.

**Risk:** none to package code; touches CI only.

---

### Task 2 — `.Rbuildignore`: stop shipping build artifacts and non-package files

**Why:** F2 (87 `.o` → as-cran NOTE + 92 MB tarball), F3 (hidden files / top-level NOTE), F5 (`Meta/` build WARNING).

**Files:** `.Rbuildignore`, `cleanup`

- [ ] **Step 2.1 — add the exclusions** (append to `.Rbuildignore`; keep the existing entries):

```
^Meta$
^build\.log$
^AGENTS\.md$
^CLAUDE\.md$
^codecov\.yml$
^planning$
^\.opencode$
^\.worktrees$
^\.lintr$
^inst/onnx/crepe/.*\.onnx$
\.(o|so|a|dll|dylib|mod)$
^src/.*/\.git$
^src/.*/\.github$
^src/.*/\.clang-format$
^src/.*/\.gitkeep$
```

  Rationale for the last four: they are submodule plumbing (`src/SPTK/.git` is a *file* containing a gitdir path, and shipping it confuses tooling), not package content. `^\.lintr$` and `^\.opencode$` are dev-tooling files that the check explicitly flags.

  **Verification already performed** (scratch copy of `.Rbuildignore`, `tools:::inRbuildignore()` — the exact predicate `R CMD build` uses): `src/assp/acf.o`, `src/RcppExports.o`, `src/SPTK/src/analysis/pitch_extraction.o`, `src/SPTK/.git`, `src/SPTK/.github`, `src/SPTK/doc/_static/.gitkeep`, `build.log`, `AGENTS.md`, `CLAUDE.md`, `codecov.yml`, `.lintr` all return `TRUE`; `R/foo.R`, `DESCRIPTION`, `src/Makevars` return `FALSE`. Directory patterns (`^Meta$`, `^planning$`, `^\.opencode$`) match the directory *entry* (R prunes the subtree before descending) — confirmed `TRUE` for `Meta`, `planning`, `.opencode`. One caveat: `\.(o|so|a|dll|dylib|mod)$` also matches `src/tcl-snack/unix/pkgIndex.tcl.dll`, which is a Tcl script, not a shared library — harmless because `^src/tcl-snack$` is already excluded.

- [ ] **Step 2.2 — clean the working tree once** (untracked, `.gitignore`d build/model artifacts):
```bash
find src -name '*.o' -not -path 'src/opensmile/build_r/*' -delete
rm -rf Meta
rm -f inst/onnx/crepe/*.onnx        # stale copy of the HF cache; re-downloaded on demand
```
  Do **not** rely on this alone — `R CMD build` recreates the sub-directory objects (F2 root cause) and the `.Rbuildignore` entries are what keep them out of the tarball.

- [ ] **Step 2.3 — top-level `.*\.o$` is already handled by R** (its `cleanup_pkg` globs `src/*.o`), but confirm no stale `src/*.o`/`*.so` remain after Step 2.2.

**Verification:**
1. `Rscript -e 'cat(tools:::inRbuildignore("src/assp/acf.o", "."), tools:::inRbuildignore("Meta/vignette.rds", "."), tools:::inRbuildignore("planning/x.md", "."), tools:::inRbuildignore("inst/onnx/crepe/full.onnx", "."), "\n")'` → note the first entry is the *file-inside-directory* form R never queries, so check the entries R actually prunes: `inRbuildignore("Meta", ".")` and `inRbuildignore("planning", ".")` must be `TRUE`, and `inRbuildignore("inst/onnx/crepe/full.onnx", ".")` must be `TRUE`.
2. `R CMD build .` → expect: no `WARNING: Removing directory 'Meta'`; tarball drops from 88 MB to roughly 4–6 MB; `tar tzf superassp_*.tar.gz | grep -cE '\.(o|so|a)$|onnx/crepe'` → `0`; `tar tzf … | grep -E 'AGENTS.md|CLAUDE.md|codecov.yml|planning/|build\.log'` → empty.
3. `R CMD check --as-cran` on that tarball → the "apparent object files/libraries", "hidden files and directories" and "Non-standard files/directories" NOTEs are gone, and `Size of tarball` is no longer a talking point.

---

### Task 3 — DESCRIPTION: CRAN-compatibility of dependencies

**Why:** F4. `Remotes` is flagged as an unknown field; `pladdrr` in `Imports` is a non-CRAN hard dependency (a blocker if this is a CRAN submission); `gsignal` is dead in `Suggests` (F1).

**Files:** `DESCRIPTION`

- [ ] **Step 3.1 — drop the vestigial remote.** Remove the `github::curso-r/torchaudio,` line from `Remotes` (nothing in the package uses `torchaudio`).
- [ ] **Step 3.2 — move `pladdrr (>= 4.8.34)` from `Imports` to `Suggests`** and keep `github::humlab-speech/pladdrr` in `Remotes`. All call sites are already guarded (`pladdrr_available()`), and `NAMESPACE` has no `pladdrr` import, so no code change is required — but this must be *proved*, see verification.
- [ ] **Step 3.3 — drop `gsignal` from `Suggests`** (no reference in code).
- [ ] **Step 3.4 — decide on `Remotes` itself.** It cannot stay in a CRAN submission (`R CMD check --as-cran` reports it as an unknown field). Recommended: keep it in the repo for `pak`/`devtools`-based development, and strip it in the submission copy (`R CMD build` then `sed '/^Remotes:/,+1d' DESCRIPTION` on the tarball, or handle it via `devtools::submit_cran()`'s workflow). Document the choice in `cran-comments.md` (already `.Rbuildignore`d, so it is not currently in the repo).

**Verification (required before accepting 3.2):**
1. `R CMD check` with `_R_CHECK_FORCE_SUGGESTS_=false` in an environment **without pladdrr** → plumbing tests still pass, `test-pladdrr-fallback.R` exercises the graceful-abort path.
2. `Rscript -e 'pak::local_install_dev_deps()'` then `devtools::test()` → `FAIL 0`.
3. `R CMD check --as-cran` → the `Strong dependencies not in the CRAN …` line and the `Unknown … 'Remotes'` line are gone from CRAN incoming feasibility (Remotes line returns only if Step 3.4 keeps the field).

---

### Task 4 — Repair the three dead bibliography DOIs

**Why:** F4. One of them is the only citation in `man/trk_covarep_vad_drugman.Rd` and is the sole DOI as-cran complains about; the other two are dead entries that will bite the next time they are cited.

**Files:** `inst/REFERENCES.bib`; regenerate `man/trk_covarep_vad_drugman.Rd` via `Rcpp::compileAttributes()` (not needed) + `devtools::document()`.

- [ ] **Step 4.1 — `Ishi2008`** (line 767): `doi = {10.1109/TASL.2007.907343}` → `doi = {10.1109/TASL.2007.910791}` (verified: same title, journal, vol 16, pp 47–56, 2008).
- [ ] **Step 4.2 — `Drugman2012VAD`** (line 1337): the entry as written does not exist. Replace the whole record with the real paper the code implements:
```bibtex
@article{Drugman2016VAD,
  author  = {Drugman, Thomas and Stylianou, Yannis and Kida, Yusuke and Akamine, Masami},
  title   = {Voice Activity Detection: Merging Source and Filter-based Information},
  journal = {IEEE Signal Processing Letters},
  year    = {2016},
  volume  = {23},
  number  = {2},
  pages   = {252--256},
  doi     = {10.1109/LSP.2015.2495219}
}
```
  then update `\insertCite{Drugman2012VAD}` → `\insertCite{Drugman2016VAD}` in `R/ssff_cpp_covarep_vad_drugman.R:5` and regenerate the Rd. (If the project prefers to keep the key stable, keep `Drugman2012VAD` as the key — the key is internal; the DOI and title are what matter. Recommended: keep the key, fix the fields, note the correction in `NEWS.md`.)
- [ ] **Step 4.3 — `White2022`** (line 763): verify against the publisher; if no DOI can be confirmed, delete the entry (it is uncited).
- [ ] **Step 4.4 — guard against recurrence**: add a small script (`inst/scripts/check_dois.R` or similar) that extracts every `doi = {}` from `inst/REFERENCES.bib` and HEAD-requests `https://doi.org/<doi>`, failing on `404`. Run it before releases. *(The 3 bad DOIs were found exactly this way in this session — 58 DOIs, 3 × 404.)*

**Verification:** the bulk DOI check returns `0` non-resolving entries; `R CMD check --as-cran` no longer reports `Found the following (possibly) invalid DOIs`.

---

### Task 5 — Root-cause and fix the coverage-job test failure

**Why:** F6. The failing assertion cannot be named from the artifacts that survive a CI run; the job has never passed.

**Files:** `tests/testthat/test-<the failing file>.R`, possibly `.github/workflows/test-coverage.yaml`

- [ ] **Step 5.1 — get the failing test name (blocking).** Land Task 1, push a scratch branch, `gh workflow run test-coverage.yaml --ref <branch>`, and read the `Failure ('test-….R:NN:NN')` block out of the log. Nothing else in this plan should be attempted for F6 before this step: the failing test set in CI is measurably different from the local one (F6 table — 58 tests that run locally never run in CI), so guessing wastes CI cycles.
- [ ] **Step 5.2 — close the environment gap.** Two candidates remain, both testable from that same CI run:
  - **native coverage flags.** CI builds every C/C++ source with `+= -O0 --coverage` (covr's default), which changes code generation for all DSP kernels; the local run had to disable them (macOS link failure). To test on Linux: `docker run --rm -v "$PWD:/pkg" rocker/r-ver:4.6.1` and run `covr::package_coverage()` there — or simply compare against a local `devtools::install()` with `withr::with_makevars(list(CXXFLAGS="-O0"), assignment="+=")`.
  - **the 17 extra skips.** Once Step 5.1 names the file, check its `skip_if_not_installed()`/`skip_if()` guards against what CI actually has installed (pak lockfile in the log lists every resolved package).
- [ ] **Step 5.3 — fix the assertion, not the harness.** Local instrumented run shows the usual suspects all pass here — `test_parallel_processing.R:168` (memory ceiling), `:135/:140` (parallel ≡ sequential), `:79/:204` (`expect_silent`), and the `skip_on_cran` files — so do not "fix" them speculatively. Whatever the real cause, the fix must assert real behaviour (identical DSP output within a defensible tolerance, or a memory ceiling that is meaningful for the user's workload) — not relax a comparison until it passes, and not delete the test.
- [ ] **Step 5.4 — caveat for local covr runs.** `covr::package_coverage()`'s `pre_clean` calls `covr:::clean_objects()`, which recursively unlinks `src/**` files ending `.o .sl .so .dylib .a .dll`. In this repo that deletes the tracked submodule file `src/tcl-snack/unix/pkgIndex.tcl.dll` (observed here). Restore with `git -C src/tcl-snack checkout -- unix/pkgIndex.tcl.dll`, or pass `pre_clean = FALSE`, after any local run. Note this is the same extension set proposed for `.Rbuildignore` in Task 2 — there the file is already covered by `^src/tcl-snack$`, so it is unaffected.
- [ ] **Step 5.5 — confirm** with a fresh `test-coverage` run: green, with the Codecov upload step succeeding.

**Verification:** the CI `test-coverage` run on the fixed commit is green, `[ FAIL 0 | ... ]`, and `codecov/codecov-action` reports a successful upload.

---

### Task 6 — Restore a usable local `--as-cran` loop

**Why:** F1. The local check aborts at the dependency step, so local results are not comparable with CI.

**Files:** none (developer environment); optionally a `Makefile`/`inst/scripts` recipe.

- [ ] **Step 6.1 — install the CRAN Suggests** into the dev library: `install.packages(c("pbapply","pbmcapply","geometry","wavelets","tuneR"))` (5 — `gsignal` is dropped in Task 3), plus the GitHub-only optional deps via `pak::local_install_dev_deps()`.
- [ ] **Step 6.2 — re-run** `R CMD build .` + `R CMD check --as-cran superassp_*.tar.gz` and record the baseline (expected after Tasks 2–4: `0 ERROR`, `0 WARNING`, ~4 NOTEs).
- [ ] **Step 6.3 — decide on the remaining NOTEs** (F7): either fix them (Rd `\_` escapes, `\usage` line widths, `stats::sd` in `.creak_extract_features`, the `unlockBinding` in `s7_methods.R`, the `%d`/`long` printf in `dataobj.c:818`) or document them as accepted. They are NOTEs, so they do not fail the check — but they are the difference between a check a reviewer reads and one they skip.

**Verification:** `00check.log` shows `Status: 0 ERROR, 0 WARNING` with only explicitly documented NOTEs, and it is reproducible on a clean library.

---

## Open items / decisions needed

1. **Is this a CRAN submission, or GitHub-only?**
   - GitHub-only → Task 3.2 (moving `pladdrr` to `Suggests`) and most of Task 3.4 are optional; `Remotes` warning can be accepted and allow-listed.
   - CRAN submission → Task 3 is mandatory (non-CRAN `Imports` is rejected), as is an explanation for the 92 MB→? tarball size.
2. **Failing coverage test — one CI run away.** Local reproduction came back green (`FAIL 0 | WARN 8 | SKIP 13 | PASS 3276`), so the exact assertion cannot be named from this machine; CI runs a measurably different test set (30 skips / 3200 passes vs 13 / 3276). Task 1 + one `workflow_dispatch` run on a scratch branch will name it. This is the only item in the plan that needs something other than this repository's local tooling (a branch push to `humlab-speech/superassp`).
3. **The 92 MB local tarball is mostly local, not shipped content:** `inst/onnx/crepe/{full,tiny}.onnx` (86.8 MB, `.gitignore`d, unused — models are fetched into the user cache), plus 8.1 MB of object files. Neither is in git, so CI's tarball is small; both still need the `.Rbuildignore` entries (Task 2) because `R CMD build` runs from the working tree.
4. **`inst/onnx/deepformants`, `inst/include`, `tests/signalfiles/AVQI`, `tests/signalfiles/DSI`** are empty (and untracked) locally; `R CMD build` strips empty directories, so any test or doc that expects them must be fixed before it is relied on. Only `tests/TRACK_NAME_NORMALIZATION_FIX.md` mentions `AVQI` today — no test depends on them.
