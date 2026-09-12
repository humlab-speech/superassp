# Win-builder amendment plan — superassp 2.9.5 → CRAN

## Artifacts under review

| Artifact | Provenance |
|---|---|
| `00install.out.txt` (1521 lines) | win-builder `R CMD INSTALL`, R Under development (unstable) r90519 ucrt, gcc/g++ 14.3.0, C++17 |
| `00check.log` (6943 lines) | win-builder `R CMD check --as-cran` on `superassp_2.9.5.tar.gz` |
| `superassp_2.9.5.tar.gz` | the actually-submitted tarball (DESCRIPTION line 37 still carries `Remotes:`) |

**Baseline: `Status: 2 ERRORs, 7 WARNINGs, 4 NOTEs`** (`00check.log:6943`).

`cran-comments.md` currently claims `0 errors | 7 warnings | 4 notes` and that `Remotes` "has been removed from the DESCRIPTION submitted here". Both statements are false against these artifacts: the errors come from a machine without `pladdrr`, and the submitted tarball still contains `Remotes`. `cran-comments.md` must not be reused until corrected.

## Baseline inventory

| # | Check | Level | Root cause | Class |
|---|---|---|---|---|
| E1 | `checking tests ... [10m]` | ERROR | `[ FAIL 15 \| WARN 7 \| SKIP 71 \| PASS 2774 ]` — all 15 are `pladdrr` absent | F1 |
| E2 | `re-building of vignette outputs` | ERROR | `ssff-jstf-io.Rmd`, `voice-quality.Rmd` call `lst_voice_report()`/`lst_vq()` | F1 |
| W1 | `whether package can be installed` | WARNING | 30 lines of "significant warnings" from vendored C/C++ (not a timeout — `cran-comments.md` misattributes this) | F5 |
| W2 | `code/documentation mismatches` | WARNING | 78 Rd files | F2 |
| W3 | `compilation flags in Makevars` | WARNING | `-Wno-register -Wno-deprecated-register` | F6 |
| W4 | `for GNU extensions in Makefiles` | WARNING | `src/Makevars`, `src/SPTK/Makefile`, `src/SPTK/tools/Makefile` | F6 |
| W5 | `pragmas in C/C++ headers and code` | WARNING | `dataPrintSink.hpp`, `csvSink.cpp`, `catch.hpp` | F6 |
| W6 | `compilation flags used` | WARNING | same flags as W3 | F6 |
| W7 | `compiled code` | WARNING | `std::cerr`, `std::cout`, `exit`, `putchar`, `puts`, `rand`, `srand` | F7 |
| N1 | `CRAN incoming feasibility` | NOTE | `Remotes` field; `pladdrr` not in mainstream repos; 11 misspellings; no `inst/WORDLIST` | F3 |
| N2 | `package subdirectories` | NOTE | `NEWS.md` headings `# superassp (development) — Consistency Refactor` (line 248) and `# Earlier Versions` (line 2188) | F8 |
| N3 | `R code for possible problems` | NOTE | `unlockBinding(fn_name, ns)` in `R/s7_methods.R` | F9 |
| N4 | `HTML version of manual` | NOTE | `superassp.html:4412 (format_apply_msg.Rd:5): Error: <fun> is not recognized!` | F4 |

---

## Findings

### F1 — Both ERRORs have one root cause: `pladdrr` is a hard dependency in tests and vignettes

`pladdrr` is in `Suggests` (`DESCRIPTION:55`, floor `>= 4.8.34`) and is GitHub-only (`Remotes: github::humlab-speech/pladdrr`). R's runtime code is correctly defensive — 20 `cli_abort()` guard sites across `R/*.R`. The test suite and two vignettes are not.

15 test failures, all `pladdrr`:

| File | Lines |
|---|---|
| `tests/testthat/test-smoke-lst.R` | 71, 83, 99 |
| `tests/testthat/test-smoke-trk.R` | 124, 138, 143, 149 |
| `tests/testthat/test_praat.R` | 14 (×4 functions) |
| `tests/testthat/test_praat_slicefunctions.R` | 15, 27 (×2 functions) |

Two vignettes fail on the same call: `ssff-jstf-io.Rmd` chunk 4 (lines 57-61, `lst_voice_report`) and `voice-quality.Rmd` (line 63, `lst_vq`). `getting_started.Rmd:80` already gates the identical call with `eval = FALSE`, so the pattern exists in-tree.

Three sub-defects make this worse than a missing skip:

- **`install_pladdrr()` does not exist.** 13 references across 10 files (`R/pladdrr_helpers.R:53`, `R/list_pladdrr_*.R`, `R/ssff_pladdrr_*.R`) tell the user to call a function that is neither defined nor exported. There is no `install_pladdrr` in `R/`, `NAMESPACE`, or `man/`.
- **Version floor is inconsistent.** `DESCRIPTION:55` requires `>= 4.8.34`; `R/list_pladdrr_pharyngeal.R:145` reports `Required version: >= 4.8.16`; the same file's doc says `>= 4.8.16` at line 72.
- **Message text is inconsistent** across the 20 guard sites — five distinct spellings/constructions of the same "not available" abort, which is why the log shows four different message shapes for the same condition.

### F2 — `code/documentation mismatches`: 78 Rd files, one mechanism

`.setup_s7_methods()` (`R/s7_methods.R`, called from `.onLoad`) converts **every** `lst_*`/`trk_*` binding in the namespace into an S7 generic:

```r
generic_fn <- S7::new_generic(name = fn_name, dispatch_args = "listOfFiles")
```

`S7::new_generic()` with no `fun` argument synthesises formals `(listOfFiles, ...)`. The 75 R files carrying hand-written `@usage` blocks (82 blocks) still document the pre-conversion signature, so the *installed* function is `function(listOfFiles, ...)` while the Rd `\usage{}` shows `function(listOfFiles, beginTime = 0, …)`. Codoc compares the two and reports a mismatch for all 78.

Because this is a mechanism-level mismatch, the count will track every new `lst_*`/`trk_*` function added. `cran-comments.md` labels it "intentional"; CRAN does not accept an intentional WARNING on a new submission.

Two candidate fixes:

- **(a) Preserve the original formals on the generic (recommended).** Pass a `fun` to `S7::new_generic()` built from `formals(original_fn)` with body `S7::S7_dispatch()`. The installed signature then equals the documented one, all 78 mismatches disappear, and the documentation stays maximally informative. Risk: dispatch calls the method with args matched against the *generic's* formals, so the method must accept the same set — `original_fn` does by construction, but the `class_any` fallback (`function(listOfFiles, ...)`) and the AVAudio method must be re-verified, and R must not report the generic as "not a function" for `args()`, `formals()`, or `utils::help` lookups.
- **(b) Fall back to deleting the `@usage` blocks** in the 75 files and letting roxygen emit `\usage{fn(listOfFiles, ...)}`. Cheap and honest, but drops default values from the rendered usage.

Take (a); keep (b) as the fallback if the dispatch experiment fails.

### F3 — CRAN incoming feasibility NOTE

- `Remotes:` is present in the built tarball (`tar -xzOf superassp_2.9.5.tar.gz superassp/DESCRIPTION` → line 37). CRAN does not recognise the field (`Unknown, possibly misspelled, fields in DESCRIPTION: 'Remotes'`). `.Rbuildignore` does not exclude it.
- `Suggests or Enhances not in mainstream repositories: pladdrr` and `Package suggested but not available for checking: 'pladdrr'` — this NOTE is unavoidable while `pladdrr` is a GitHub-only optional dependency. `pladdrr` **must stay in `Suggests`**: R CMD check's `checking dependencies in R code` reports `'loadNamespace' or 'requireNamespace' call not declared from: 'pladdrr'` if the field is dropped while 20 guarded call sites remain.
- 11 possibly-misspelled words (`AsspDataObj`, `JSTF`, `Praat`, `SSFF`, `Wrassp`, `cepstral`, `formant`, `lst`, `openSMILE`, `toolkits`, `trk`). No `inst/WORDLIST` exists — every one of these is a correct domain term that a WORDLIST would silence.
- `New submission` is expected and permanent.

### F4 — HTML manual NOTE: unescaped angle bracket in an Rd title

`R/sptk_helpers.R` has the roxygen title `Emit a consistent "Applying \if{html}{\out{<fun>}}()" progress message`, which renders a literal `<fun>` element into `superassp.html` at line 4412. Fix the title text; do not touch the `\description{}` (whose `\code{fun()}` is fine).

### F5 — Install WARNING is compiler warnings, not a timeout

`cran-comments.md` attributes W1 to "exceeds the check's default time budget". The log shows it is R's "significant warnings" list. Six distinct sites:

| Site | Warning |
|---|---|
| `src/opensmile/src/lldcore/melspec.cpp:408` | `-Walloc-size-larger-than=` (arg range folded to ~2^64) |
| `src/opensmile/src/include/smileutil/smileUtil.h:731`, `src/opensmile/src/smileutil/smileUtil.c:2922` | `-Wstrict-prototypes` |
| `src/opensmile/src/include/rapidjson/document.h:110,118` | `std::iterator` deprecated (removed in C++17) |
| `src/opensmile/src/include/core/configManager.hpp:486-490`, `smileComponent.hpp:47-67`, `progsrc/smileapi/SMILEapi.cpp:37-46` | `-Wreorder` (10 lines, 3 headers/1 source) |
| `vat_srh_pitch.cpp:94`, `vat_gci_se_vq.cpp:157`, `vat_voice_quality.cpp:119`, `vat_creak_detect.cpp:244`, `vat_creak_detect.cpp:388` | `arma::Mat::max(uword&)` deprecated — **superassp's own code** |
| `src/SPTK/third_party/REAPER/core/float_matrix.cc:53` | `-Wuninitialized` (real UB) |

The five arma sites are ours and trivially fixable. The rest is vendored — but the usual escape hatch of adding `-Wno-…` is unavailable: any warning-suppressing flag trips W6/W3 (see F6).

### F6 — Three WARNINGs share the "flags and Makefiles" cause

- `PKG_CXXFLAGS = $(SHLIB_OPENMP_CXXFLAGS) -Wno-register -Wno-deprecated-register` in both branches of `src/Makevars`. The comment justifies the flag by Apple Clang rejecting the C++-removed `register` storage-class specifier in `SPTK/third_party/Snack/jkGetF0.cc` and `sigproc.cc`. On GCC 14.3 (win-builder) the flag is unnecessary — GCC only warns — and it is what produces W3 *and* W6.
- GNU make extensions: `src/Makevars` (`ifeq`, `$(shell uname -s)`, `$(shell xcrun …)`, `$(shell pwd)`), plus vendored `src/SPTK/Makefile` and `src/SPTK/tools/Makefile` (`+=`, `$(wildcard)`). The two SPTK Makefiles are SPTK's own standalone build files and are unused by R — `src/Makevars` lists every SPTK source explicitly.
- Pragma WARNING: `#pragma GCC diagnostic ignored` in `src/opensmile/src/include/iocore/dataPrintSink.hpp:33`, `src/opensmile/src/iocore/csvSink.cpp:178` ("suppressing important diagnostics"), and `src/opensmile/progsrc/include/tests/catch2/catch.hpp:42-46`. `catch.hpp` is a test-framework header never compiled into the package.

The `register` problem, not the flag, is the real defect: patching the two Snack files removes the need for both flags, which closes W3 and W6 together.

### F7 — Compiled-code WARNING

| Symbol | Object |
|---|---|
| `std::cerr` | `SPTK/src/utils/sptk_utils.o` (`sptk_utils.cc:510`) |
| `std::cout` | `opensmile/build_r/libopensmile.a` (×4) |
| `exit` | (aggregated) |
| `putchar`, `puts` | `tandem/tandem_64/feature.o`, `tandem/tandem_64/mScaleInten.o` |
| `rand`, `srand` | `tandem/tandem_64/*.o`, `snack_formant.o`, `libopensmile.a` |

Vendoring split matters for how these get fixed: `src/SPTK` and `src/tandem` are `humlab-speech` forks (submodules, currently at `fc92384` and `e98f436`) and can be patched upstream + pointer-bumped; `src/opensmile` is vendored in-tree and patchable directly; `src/snack_formant.cc` is ours (`src/` root).

`rand`/`srand` are the one group that must **not** be swapped for R's RNG — `snack_formant.cc:406,612-613` uses them in a real DSP path, and changing the sequence changes numeric output. That constraint is absolute in this repo.

### F8 — `NEWS.md` version-section NOTE

`# superassp (development) — Consistency Refactor` (line 248) and `# Earlier Versions` (line 2188) are not `# <pkg> <version>` headings, so R's news parser cannot extract a version. Either version the development section or remove the headings.

### F9 — `unlockBinding` NOTE

Intrinsic to the `.onLoad` generic-conversion design. Not worth restructuring for a NOTE; explain it in `cran-comments.md`.

---

## Constraints

- **No change to DSP numeric output.** Any change under `src/` must be semantics-preserving; PRNG substitution is out (F7).
- Vendored-fork edits land in `humlab-speech/SPTK` and `humlab-speech/tandem` upstream, then the submodule pointer is bumped here. In-tree `src/opensmile` is edited directly.
- `devtools::document()` after every roxygen edit.
- Do not re-pin or delete tests to make the check pass; the `pladdrr` failures must become *skips*, not deletions.
- Preserve `pladdrr` in `Suggests` (F3).

## Order

Phase 1 closes both ERRORs and is the gate for any further submission. Phases 2–3 are required for 0 WARNING. Phase 4 is NOTES.

---

### Phase 1 — Reproduce the win-builder condition locally

**Why:** the maintainer's library has `pladdrr`, so a local `--as-cran` run cannot see E1/E2. Every task below is unverifiable until the local check runs without it.

- [ ] **Task 1.1 — check as CRAN does.** Run `R CMD build .` then `R CMD check --as-cran superassp_2.9.5.tar.gz` with `pladdrr` unavailable (`_R_CHECK_FORCE_SUGGESTS_=false`, and a library where `pladdrr` is not installed). Record the baseline; it must show 2 ERRORs.
- [ ] **Task 1.2 — pin the baseline to this log.** Copy the fresh `00check.log` result into `planning/` alongside this file so later phases diff against a known state. Add `00install.out.txt` to `.gitignore` (untracked; `00check.log` is already handled).

**Verification:** local check reproduces `2 ERRORs`, and the failing test list matches the 15 rows in F1 exactly.

---

### Phase 2 — Remove the `pladdrr` hard dependency from tests and vignettes (closes E1, E2)

- [ ] **Task 2.1 — one guard helper.** Add `tests/testthat/helper-pladdrr.R` exposing a single `skip_without_pladdrr()` wrapper over the existing house convention (`skip_if_not(pladdrr_available(), "pladdrr not available")`, as used in `R/pladdrr_helpers.R:4` and `tests/testthat/test-lst_voice_report.R`, `test-dysprosody.R:4`). Four test files currently lack any guard; do not invent a second convention.
- [ ] **Task 2.2 — guard the 15 failing tests.** `test-smoke-lst.R` (71, 83, 99), `test-smoke-trk.R` (124, 138, 143, 149), `test_praat.R` (`praat_funs` loop, line 14), `test_praat_slicefunctions.R` (`slicefunctions` loop, lines 15 and 27). The loops must skip per-iteration so the non-`pladdrr` functions in the same loop (`trk_pitch_cc` etc. included) still run where they can.
- [ ] **Task 2.3 — gate the two vignette chunks.** `vignettes/ssff-jstf-io.Rmd` chunk 4 (57-61) and `vignettes/voice-quality.Rmd` (61-71) → `eval = requireNamespace("pladdrr", quietly = TRUE)`, matching `getting_started.Rmd`'s existing conditional style, with one sentence of prose naming the dependency and pointing at the install route from Task 3.3.
- [ ] **Task 2.4 — collapse the guard messages.** Replace the 20 `cli_abort()` sites with one internal `pladdrr_unavailable(caller)` constructor in `R/pladdrr_helpers.R` producing a single message. Derive the version floor from one constant so F1's `4.8.16` vs `4.8.34` divergence cannot recur, and reconcile those two values against `DESCRIPTION:55` and `R/list_pladdrr_pharyngeal.R:72,145`.

**Verification:** `R CMD check --as-cran` with `pladdrr` absent gives `checking tests ... OK` and `re-building of vignette outputs ... OK`; the test summary shows the 15 as `SKIP`. Re-run with `pladdrr` present (local `NOT_CRAN=true` run) and confirm the same 15 tests execute and pass, i.e. the guard is a skip, not a mask.

---

### Phase 3 — DESCRIPTION and submission metadata (closes N1)

- [ ] **Task 3.1 — remove `Remotes:` from `DESCRIPTION`.** Verified present in the submitted tarball at line 37.
- [ ] **Task 3.2 — add `inst/WORDLIST`** with the 11 domain terms from F3, one per line, so the spell NOTE resolves.
- [ ] **Task 3.3 — give users a working install route.** `pladdrr` cannot be installed via `install.packages()`; the `Remotes` field was the only machine-readable pointer. Whatever Task 2.4's message says must be a command that actually works (`remotes::install_github("humlab-speech/pladdrr")` or `pak::pkg_install("humlab-speech/pladdrr")`), and `README.md` plus the vignette prose must agree with it.
- [ ] **Task 3.4 — accept the residual `pladdrr` NOTE deliberately.** `Suggests or Enhances not in mainstream repositories` stays. Document it, with the reason (`pladdrr` is GitHub-only; every call site guarded; the package installs and checks without it).
- [ ] **Task 3.5 — rewrite `cran-comments.md` against the real log.** It currently reports `0 errors | 7 warnings | 4 notes`, misattributes the install WARNING to a timeout, and asserts a `Remotes` removal that did not happen. Regenerate it only after Phase 4, from the final `00check.log`.

**Verification:** `tar -xzOf superassp_*.tar.gz superassp/DESCRIPTION | grep Remotes` returns nothing; `checking CRAN incoming feasibility` NOTE no longer lists `Remotes` or the misspellings.

---

### Phase 4 — `code/documentation mismatches` (closes W2)

- [ ] **Task 4.1 — preserve original formals on the S7 generic.** In `R/s7_methods.R`'s `.convert_to_s7_generic()`, construct `S7::new_generic(name = fn_name, dispatch_args = "listOfFiles", fun = <fn with formals(original_fn) and body S7::S7_dispatch()>)`. Keep the `inherits(original_fn, "S7_generic")` early return, the method registrations, the attribute re-attachment, and the `unlockBinding`/`lockBinding` pair untouched.
- [ ] **Task 4.2 — verify dispatch before believing it.** For a representative sample spanning both families and every input class — a plain-signature `trk_*` (`trk_acf`), an S7-generic-family `lst_*` (`lst_vq`), a function with a `...`-forwarding wrapper (`processMediaFiles_LoadAndProcess`), and one with alias/attribute metadata — confirm: character vector input takes the original method; `AVAudio` input takes the tempfile method; an unsupported input takes the `class_any` fallback; defaulted arguments are actually applied when omitted; `formals()` and `args()` report the full signature.
- [ ] **Task 4.3 — fallback only if 4.2 fails.** Delete the 82 `@usage` blocks across the 75 R files and re-document, accepting `\usage{fn(listOfFiles, ...)}`. Do not do this speculatively; it loses information.

**Verification:** `checking for code/documentation mismatches ... OK` in a fresh check, **and** the Task 4.2 dispatch matrix passes, **and** `devtools::test()` is unchanged from its pre-edit pass count. A green codoc with broken AVAudio dispatch is a worse outcome than the WARNING.

---

### Phase 5 — Build-system WARNINGs (closes W1, W3, W4, W5, W6)

Ordered so that each step removes a flag rather than adding one.

- [ ] **Task 5.1 — delete the `register` storage-class specifier** in `SPTK/third_party/Snack/jkGetF0.cc` and `sigproc.cc` (fork `humlab-speech/SPTK`), bump the submodule pointer. Mechanical and semantics-preserving.
- [ ] **Task 5.2 — drop `-Wno-register -Wno-deprecated-register` from `src/Makevars`** once 5.1 lands, closing W3 and W6. Also verify macOS still builds (the original justification for the flag).
- [ ] **Task 5.3 — portabilise `src/Makevars`.** Move the Darwin `uname`/`xcrun` branch into `configure` (which already exists and already gates on CMake) and ship `src/Makevars.in`, so the shipped `Makevars` contains no `ifeq`/`$(shell)`. Replace `$(shell pwd)` in `PKG_LIBS` with a path R's make resolves relative to `src/`. Closes the `src/Makevars` third of W4.
- [ ] **Task 5.4 — `src/SPTK/Makefile`, `src/SPTK/tools/Makefile`.** They are SPTK's standalone build files, unused by `src/Makevars`. Either portabilise them upstream in the fork or add both to `.Rbuildignore`. Prefer excluding: they are not part of the R build. Closes W4.
- [ ] **Task 5.5 — remove the diagnostic-suppressing pragmas.** `dataPrintSink.hpp:32-37` and `csvSink.cpp:177-182` — fix the underlying `-Wformat-security`/`-Wformat-overflow` trigger instead of suppressing it, or drop the pragma and let the warning surface for a proper fix. Closes two thirds of W5.
- [ ] **Task 5.6 — stop shipping `catch.hpp`** (`.Rbuildignore`: `^src/opensmile/progsrc/include/tests`). A test-framework header that is never compiled. Closes the rest of W5.
- [ ] **Task 5.7 — fix the five arma deprecations (our code).** `vat_srh_pitch.cpp:94`, `vat_gci_se_vq.cpp:157`, `vat_voice_quality.cpp:119`, `vat_creak_detect.cpp:244`, `vat_creak_detect.cpp:388` use the `max(uword&)` overload, which returns the value *and* writes the index. Replace with `index_max()`/`max()` preserving both outputs. Same fix required in `src/vat_*.cpp` only — the `.max()` no-arg overloads (`vat_gci_se_vq.cpp:97,222`, `vat_creak_detect.cpp:329,336,340,396,551`) are not deprecated and must not be touched.
- [ ] **Task 5.8 — vendored compiler warnings.** `smileUtil.h:731`/`smileUtil.c:2922` → `(void)` prototypes. `configManager.hpp:486-490`, `smileComponent.hpp:47-67`, `SMILEapi.cpp:37-46` → reorder member-initialiser lists (semantics-preserving). `rapidjson/document.h:110,118` → update the bundled rapidjson to a C++17-clean version or replace the `std::iterator` base. `melspec.cpp:408` → hoist the inlined call into a local so GCC cannot fold the argument range. REAPER `float_matrix.cc:53` → initialise `data_`.
- [ ] **Task 5.9 — re-check; if W1 survives, isolate the remainder** with a per-file bisect rather than a blanket `-Wno-`.

**Verification:** `checking whether package 'superassp' can be installed ... OK` with no "significant warnings" list; W3–W6 gone; `R CMD INSTALL` still succeeds on macOS (the flag removal is the risky step) and on win-builder's GCC 14.3.

---

### Phase 6 — Compiled-code and remaining NOTEs (W7, N2, N3, N4)

- [ ] **Task 6.1 — `std::cerr` → `Rprintf`/`REprintf`** in `SPTK/src/utils/sptk_utils.cc:510` (fork + pointer bump).
- [ ] **Task 6.2 — `std::cout` ×4 in the opensmile sources** compiled into `libopensmile.a`. Route through openSMILE's logging (`SMILE_MSG`) or R's console. Requires locating all four call sites, since the check reports only the archive.
- [ ] **Task 6.3 — `printf`/`putchar` in `tandem/tandem_64/feature.cpp`, `mScaleInten.cpp`** — these are the `echo`-gated progress prints. Route to `Rprintf` (fork + pointer bump), or pass `echo = 0` permanently and confirm the symbols drop.
- [ ] **Task 6.4 — `exit`** in the archive. Replace with an exception or a return code; confirm the error path still terminates the analysis cleanly.
- [ ] **Task 6.5 — `rand`/`srand` are deliberately retained.** `snack_formant.cc:406,612-613` and the tandem RNG feed DSP output. Substituting R's RNG changes results. If W7 survives Tasks 6.1–6.4 on these symbols alone, accept and document in `cran-comments.md` with the numeric-stability constraint as the reason. **Decision point — do not silently swap the PRNG.**
- [ ] **Task 6.6 — `NEWS.md` (N2).** Version the development section (line 248) or fold it into a released heading; remove `# Earlier Versions` (line 2188), keeping the pointer text.
- [ ] **Task 6.7 — `format_apply_msg` title (N4).** Remove the literal `<fun>` from the roxygen title in `R/sptk_helpers.R`; re-document and confirm `superassp.html` no longer trips the validator.
- [ ] **Task 6.8 — `unlockBinding` (N3).** Accept and document. Intrinsic to the load-time generic conversion; restructuring for a NOTE is not worth the risk to Task 4.1.

**Verification:** `Status:` line shows 0 ERRORs; W7 either gone or reduced to the `rand`/`srand` symbols; N2 and N4 gone; N1 and N3 documented in `cran-comments.md`.

---

### Phase 7 — Resubmit

- [ ] **Task 7.1 — `R CMD build .`, verify the tarball DESCRIPTION** has no `Remotes` and carries `inst/WORDLIST`, then check on win-builder R-devel and R-release.
- [ ] **Task 7.2 — rewrite `cran-comments.md` from the final logs** (supersedes Task 3.5's draft): real status line, the two accepted NOTEs with reasons, the `pladdrr` situation, and the `rand`/`srand` constraint if it survives 6.5.
- [ ] **Task 7.3 — update `NEWS.md`** with the amendment summary and bump the version.

**Verification:** win-builder returns `0 ERRORs` on both R-devel and R-release; every remaining WARNING/NOTE is one that `cran-comments.md` names and justifies.

---

## Residual, accepted after the plan

| Item | Reason |
|---|---|
| `Suggests or Enhances not in mainstream repositories: pladdrr` | GitHub-only optional dependency; every call site guarded; package checks without it (F3) |
| `New submission` | permanent for a first submission |
| `unlockBinding` NOTE | intrinsic to load-time S7 conversion (F9) |
| `rand`/`srand` compiled-code entries | swapping the PRNG changes DSP output — prohibited (F7, Task 6.5) |

## What this plan does not do

- It does not touch DSP numerics — Task 5.7 and 5.8 are restricted to initialiser order, prototype syntax, and argument hoisting.
- It does not delete or re-pin any test; the 15 `pladdrr` failures become skips under one shared helper.
- It does not remove `pladdrr` from `Suggests`.

---

# Implementation record — 2026-09-12

## Result

`R CMD check --as-cran` on macOS (R 4.6.1, Apple clang 21), run under the
win-builder condition, moved from the reviewed baseline to **0 ERRORs,
0 WARNINGs, 4 NOTEs**:

| | baseline (win-builder) | after |
|---|---|---|
| ERRORs | 2 | 0 |
| WARNINGs | 7 | 0 |
| NOTEs | 4 | 4 |

`checking tests`, `re-building of vignette outputs`, `checking examples`,
`checking examples with --run-donttest` and `whether package can be installed`
all report OK.

Final artifact: `superassp_3.0.0.tar.gz`. The release was bumped from 2.9.5 to
3.0.0 by the API migration described below, which is the only user-visible
behaviour change in the submission.

## How pladdrr-absence was reproduced locally

The maintainer's library has `pladdrr` installed, which is why every earlier
local check hid both ERRORs. Faithfully reproducing CRAN/win-builder needs
`pladdrr` to be *findable but unloadable*: absent entirely makes R CMD check
report `Package suggested but not available`, which is a different failure.

A stub package is installed into a private library and put on `R_LIBS`:

```
/tmp/nopladdrr/pladdrr     # DESCRIPTION Version 4.8.34 + .onLoad() that stops
R CMD INSTALL --no-docs --no-help --no-test-load -l /tmp/nopladdrr <stub>
R_LIBS=/tmp/nopladdrr R CMD check --as-cran superassp_3.0.0.tar.gz
```

`find.package("pladdrr")` succeeds (so the dependency check is satisfied) while
`requireNamespace("pladdrr")` returns FALSE (so tests and vignettes take their
skip paths). The real install is never touched. This recipe is worth keeping.

## Per-task outcome

**Phase 1** — done. Baseline reproduced at the testthat level before any edit
(`FAIL 15`, same messages as win-builder) and at the check level via the stub.

**Phase 2 (E1, E2)** — done. `tests/testthat/helper-pladdrr.R` adds
`skip_without_pladdrr()`; the 15 in-file conditions were migrated to it and the
4 unguarded test files now use it. Both vignette chunks are gated on
`has_pladdrr`. The 17 abort sites collapse into `pladdrr_unavailable()`, and
`install_pladdrr()` — which never existed — is gone from every message.

**Phase 3 (N1)** — done. `Remotes:` removed from DESCRIPTION, `inst/WORDLIST`
added. The `pladdrr`-not-in-mainstream NOTE is retained deliberately: it cannot
be removed while `pladdrr` is an optional GitHub-only dependency, and the field
must stay in `Suggests` or `requireNamespace` calls become undeclared.

**Phase 4 (W2)** — done, via option (a). The generic now carries the original
function's formals, and the helper methods reuse the generic's signature (S7
requires an exact match when the generic has no `...`, which the first attempt
missed — it dropped conversion to 3/86 functions until fixed). Conversion
coverage is back to 77/86, identical to before. Two genuine doc/code
divergences were fixed separately: `processMediaFiles_LoadAndProcess`'s stale
`@usage` and four `trk_*` blocks documenting a `listOfFiles = NULL` default.

**Phase 5 (W1, W3, W4, W5, W6)** — done. `register` removed from the Snack
sources (fork), which let both `-Wno-*` flags go; `src/Makevars` is now
generated by `configure` from `src/Makevars.in` and contains no GNU make
constructs; SPTK's standalone Makefiles and `catch.hpp` are excluded. The
arma, REAPER, rapidjson, reorder, prototype and pragma diagnostics were fixed
at the source.

Two findings changed the plan during execution:

- The GNU-extensions check **recreates** what `configure` produces: it scans
  the post-install source tree, and `.Rbuildignore` cannot exclude a directory
  that `configure` regenerates. The CMake tree was therefore moved to a hidden
  directory (`src/opensmile/.build_r`), since the scan provably does not
  descend into hidden directories. Verified with a purpose-built fixture.
- The flags check reads the *effective* Makevars values, not the file text,
  while the GNU-extensions check reads the text. Verified with fixtures.

**Phase 6 (W7, N2, N4)** — N2 and N4 done. W7 is downgraded to a NOTE and
accepted: `rand`/`srand` sit in the bundled Tandem/openSMILE/Snack DSP paths and
replacing them with R's RNG would change numeric output. `sprintf` in Tandem
was converted to `snprintf`, which removed the `[v]sprintf` finding.

**Phase 7** — `cran-comments.md` rewritten from the final log.

## New findings raised during implementation

1. **`test-trk-attributes.R` was passing vacuously.** It asserts every exported
   `trk_*` defaults `toFile = FALSE`. Because the installed functions were
   generics with `(listOfFiles, ...)`, `formals(fn)$toFile` was always absent,
   so the check never examined anything. Preserving the real formals exposed
   the real gap: **42 of the 64** `trk_*` wrappers still defaulted to `TRUE`.
   (An early read of the check log suggested 10 — that was a truncated excerpt
   of the printed list.) Raised with the maintainer, who chose to complete the
   sweep; see "API migration" below.
2. **Real bugs in vendored code.** REAPER's `FloatMatrix` copy constructor read
   `data_` before initialising it, so `clear()` could `delete[]` a garbage
   pointer. `tandem` wrote paths with `sprintf` into 255-byte buffers.
3. **`install_pladdrr()` did not exist** — 13 error messages pointed users at an
   undefined, unexported function.
4. **Two doc claims contradicted the code** — `list_pladdrr_pharyngeal`'s
   version floor and `list_pladdrr_dysprosody`'s "pladdrr is in Imports".

## Vendored fixes — published

Both fixes are now durable. The fork convention turned out to be a dedicated
pin branch rather than the default branch, so the SPTK work went to
`superassp-pin` (which sat exactly on the previously pinned `fc92384`) and the
Tandem work to `master` (which sat exactly on the previously pinned `e98f436`).
Both were fast-forwards:

| Submodule | Pushed to | From → to |
|---|---|---|
| `src/SPTK` | `superassp-pin` | `fc92384` → `1ddf910` |
| `src/tandem` | `master` | `e98f436` → `817652e` |

The superproject pointers were then bumped and staged. This mattered more than
it first appeared: every CI workflow checks out with `submodules: recursive`,
so a pristine clone would otherwise have hit Apple Clang's hard error on
`register` once the `-Wno-register` flags were removed — the macOS CI leg would
have failed to compile despite the tarball being correct.

**Latent inconsistency worth fixing separately:** this checkout's
`.git/config` carries `submodule.src/SPTK.url = sp-nitech/SPTK.git`, an
upstream URL that contradicts `.gitmodules` (which names the humlab fork, and
is what CI uses). It is local-only, but it is what made the first read of the
submodule layout wrong. It was left alone rather than changed.

## API migration — `trk_*` `toFile` default (decision: migrate)

`test-trk-attributes.R` has always asserted that every exported `trk_*`
defaults to `toFile = FALSE`, per the project spec. It was passing vacuously:
the installed functions were S7 generics, so `formals(fn)$toFile` was always
absent and the check examined nothing. Preserving the real formals exposed the
gap — and the real number was **42 of 64**, not the 10 a truncated log excerpt
suggested.

Resolved by completing the sweep: all 42 defaults flipped, `enforce_toFile`
set to `TRUE`, and the release bumped to **3.0.0** because it changes what a
bare `trk_acf("f.wav")` does.

The edit was applied by script with a deliberately narrow scope — the
signature and the `@usage` block only:

* An earlier, broader attempt also rewrote prose such as
  "If \code{toFile = TRUE}: integer count of files written", inverting six
  sentences that describe behaviour *when* `toFile` is TRUE, and added a
  trailing newline to a file that had none. Both were caught by diffing
  against a pre-edit copy and the whole edit was redone with the narrower
  scope. Final diff: 84 changed lines, all `toFile` defaults, zero collateral.
* The `@param toFile` prose still said "Default `TRUE`" in 48 places. That
  phrase also appears on other parameters (`verbose`, `validate`,
  `auto_unbox`, `parallel`), so a second script rewrote it only inside
  `@param toFile` blocks, using literal rather than regex matching after a
  first attempt silently matched nothing.

`inst/QUICK_REFERENCE.md`'s new-function skeleton was the one shipped document
still showing `toFile = TRUE` as a default; it was aligned. Explicit
`toFile = TRUE` calls in `@examples` were left alone — they demonstrate the
file-writing mode deliberately.

## Verification after the migration

* `devtools::test()` — `FAIL 0 | WARN 8 | SKIP 13 | PASS 3282`, including the
  newly-enforced `toFile` contract check.
* `R CMD check --as-cran`, version 3.0.0, under the `pladdrr` stub condition
  recorded above.

Remaining: re-run win-builder on R-devel and R-release to confirm parity with
the local result.

## Accepted residuals

| Item | Reason |
|---|---|
| `Suggests ... not in mainstream repositories: pladdrr` | Optional GitHub-only dependency; every call site guarded; checks without it |
| `New submission` | Permanent for a first submission |
| `unlockBinding` in `R/s7_methods.R` | Intrinsic to the load-time S7 conversion |
| Compiled code: `rand`, `srand`, `printf`, `stderr`, `stdout`, `std::cerr` | Bundled DSP sources; the RNG calls cannot change without altering numeric output |
| `HTML version of manual` NOTE | Locally it only reports that `tidy` is too old to run the validation; the one real finding (`<fun>` in `format_apply_msg`) is fixed |

