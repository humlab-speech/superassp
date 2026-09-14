# GitHub Actions and win-builder remediation — 2026-09-14

Supersedes nothing; this is the follow-up to `2026-09-13-compiled-code-amendment.md`,
which landed the compiled-code fixes and left two residuals: the fork commits unpushed and
the win-builder run outstanding.

## Artifacts

| Artifact | Provenance |
|---|---|
| `00check.log` (121 lines) | win-builder `R CMD check`, R-release, `superassp_3.0.0.tar.gz`, submitted 2026-09-14 08:10, `gcc/g++ 14.3.0` |
| `00install.out.txt` (221 KB) | same run's install log |
| `gh run list --repo humlab-speech/superassp` | whole-workflow failures on every push since 2026-09-12 |
| `gh run view 34812076499 --log-failed` | current head `a53ab43` — fails in `actions/checkout` |
| `gh run view 34707528328 --log-failed` | previous head `4f65d32` — fails in `setup-r-dependencies` |

## Where the compiled-code amendment landed

Win-builder now reports **`Status: 1 WARNING, 2 NOTEs`**, and `checking compiled code ...
OK` — the NOTE this work was aimed at is gone, on the Windows toolchain, in the real
check. Everything else passes there: `checking tests ... [11m] OK`,
`re-building of vignette outputs ... OK`, `examples ... OK`, `PDF/HTML manual ... OK`,
`line endings in C/C++/Fortran sources ... OK`, `pragmas ... OK`,
`compilation flags ... OK`. The two NOTEs are the known, documented ones (CRAN incoming
feasibility; `unlockBinding()` in `R/s7_methods.R`).

## Findings

### F1 — every workflow dies in `actions/checkout` (current head, all four workflows)

`R-CMD-check`, `test-coverage`, `lintr` and `pkgdown` all fail in 22–31 s. The failing
step is `actions/checkout@v4`, whose `submodule update --init --force --depth=1 --recursive`
fetches `src/SPTK` from `https://github.com/humlab-speech/SPTK.git` and `src/tandem` from
`https://github.com/humlab-speech/tandem.git`.

The pins committed in `43da7a5` do not exist on those forks:

```
git ls-remote https://github.com/humlab-speech/SPTK.git   1b3a503  ->  no refs
git ls-remote https://github.com/humlab-speech/tandem.git a1fd952  ->  no refs
```

`1ddf910` (previous SPTK pin) is on the fork's `superassp-pin` branch and `817652e`
(previous tandem pin) is tandem's `master` tip, so the previous state was fetchable. The
two new commits were made in the local submodule checkouts and never pushed — exactly the
residual flagged in the 2026-09-13 record. This is a hard blocker: nothing else in CI can
run until the pins resolve.

### F2 — after F1, dependency resolution fails on `pladdrr` (since 2026-09-12)

`gh run view 34707528328` (macOS and Windows jobs) dies in `setup-r-dependencies`:

```
! Could not solve package dependencies:
* deps::.: Can't install dependency pladdrr (>= 4.8.34)
* pladdrr: Can't find package called pladdrr.
```

`Remotes: github::humlab-speech/pladdrr` was removed from `DESCRIPTION` in the 3.0.0
remediation because CRAN does not accept the field. `setup-r-dependencies` has no
`dependencies:` input in `R-CMD-check.yaml`, so it solves `Suggests` against the
configured repositories, where `pladdrr` does not exist. The library solves nothing else
wrongly: the `any::sessioninfo` / `any::rcmdcheck` / `any::Rcpp*` "dependency conflict"
lines are downstream of that one failure.

This means CI has been red for reasons unrelated to the compiled code since
`4f65d32` (2026-09-12), i.e. the check that should have caught yesterday's work has not
run at all.

### F3 — one install WARNING on Windows, introduced by the amendment

`00check.log:41-58`:

```
* checking whether package 'superassp' can be installed ... WARNING
Found the following significant warnings:
  smileUtil.c:2396|2411|2424|2438|2455|2466: warning: unknown conversion type character 'z' in format [-Wformat=]
  smileUtil.c:2396:99: warning: format '%s' expects argument of type 'char *', but argument 2 has type 'long long unsigned int'
  ...: warning: too many arguments for format [-Wformat-extra-args]
```

All six are the `smilePcm: ... %zu ... '%s' ...` wave-file error messages that the
conversion rewrote from `fprintf(stderr, ...)` to `smile_console_error(...)`. The call
sites are correct — `%zu` with `sizeof(...)`, `%s` with `filename` — and they compiled
warning-free as `fprintf` calls because GCC already knows MinGW's `gnu_printf` semantics
for the `fprintf` family. My `SMILE_CONSOLE_FORMAT` used
`__attribute__((format(printf, M, N)))` unconditionally, which tells GCC to apply MSVCRT
`printf` semantics, where `%zu` does not exist; the `%s`/extra-arguments warnings are
cascades from that. Runtime behaviour is unaffected: the shim formats with the C library's
`vsnprintf`, and Rtools 4.5 is UCRT, whose `vsnprintf` implements `%zu`.

Fix shape (mirrors R's own `R_PRINTF_FORMAT` in `R_ext/Print.h`, validated by
preprocessing a stub under `-D_WIN32 -D_UCRT` and `-D__MSVCRT_VERSION__=0xE10`):

| Build | Attribute produced |
|---|---|
| Linux/macOS | `format(printf, 1, 2)` — unchanged, `%zu` still checked |
| Windows + GCC + UCRT | `format(gnu_printf, 1, 2)` — `%zu` accepted |
| Windows + clang | `format(printf, 1, 2)` — R's own choice for clang |
| Windows + pre-UCRT msvcrt | none — nothing to validate against |

### F4 — the CI warning allowlist can no longer be trusted

`R-CMD-check.yaml`'s `Fail on unexpected R CMD check WARNINGs` step allowlists

```
code/documentation mismatches | GNU extensions in Makefiles | line endings in Makefiles |
pragmas in C/C++ headers and code | compiled code | whether package .* can be installed |
compilation flags (in Makevars|used)
```

Four of those no longer describe reality: `compiled code` is now `OK` on both platforms,
`pragmas` was fixed at the source in the 3.0.0 remediation, `-Wno-register` is gone, and
`code/documentation mismatches` was closed by reusing the original formals on the S7
generic. Keeping them allowlisted means CI silently tolerates exactly the class of
regression this submission worked to remove — F3 would have passed CI unnoticed, because
the check name it lands under (`whether package .* can be installed`) is allowlisted
*for build time*.

### F5 — the fork push path is not configured in this checkout

`src/SPTK` has only `remote.origin = https://github.com/sp-nitech/SPTK.git` (upstream);
the `humlabfork/*` refs in the checkout are stale, and the fork remote is absent. The
superproject's `.gitmodules` points at the fork, so CI and the local checkout disagree
about where SPTK comes from. `gh` is authenticated as `FredrikKarlssonSpeech` with `repo`
and `workflow` scopes, so pushing is possible once a fork remote exists.

### F6 — `pkgcheck.yaml` cannot ever run

It is triggered by `workflow_dispatch` and `push` to `main`, but the default branch is
`master`:

```
on:
  workflow_dispatch:
  push:
    branches:
      - main
```

`gh run list --workflow pkgcheck.yaml` is empty — it has never executed. Either the
branch filter is a leftover from a `main`-era rename and should read `master`, or the
workflow is abandoned and should be deleted. `lint-changed-files.yaml` is
`pull_request`-only, which is by design and not a fault.

## Constraints

- `Remotes:` must stay out of `DESCRIPTION`; the `pladdrr` install has to be expressed in
  the workflows, not in package metadata.
- The pins must name commits that exist on the forks; the package cannot be built from a
  fresh clone otherwise (this is what CRAN's tarball does not care about, but CI and
  `R CMD build` from a clone do).
- The allowlist may only shrink to entries that are still justified; the CI gate is worth
  more than a green tick.
- The `pladdrr` tests must remain skippable — win-builder and CRAN run without it, and
  they already do (`testthat` reports 0 failures with and without it).

## Phases

### Phase 1 — publish the fork commits (unblocks everything)

- [ ] **Task 1.1 — SPTK.** Add the fork remote and publish the branch:
  ```sh
  git -C src/SPTK remote add humlabfork https://github.com/humlab-speech/SPTK.git   # or set-remote if stale
  git -C src/SPTK push humlabfork cran-console-output:superassp-pin
  ```
  `1b3a503` descends from `1ddf910`, which is the current `superassp-pin` tip, so this is a
  fast-forward and keeps the existing branch convention (the superproject pins a SHA, and
  the SHA lives on a named branch). Do **not** force-push.
- [ ] **Task 1.2 — tandem.** `a1fd952` is one commit on top of the fork's `master` tip
  (`817652e`), so a plain fast-forward push:
  ```sh
  git -C src/tandem push origin master
  ```
- [ ] **Task 1.3 — verify reachability**, which is the actual acceptance criterion:
  ```sh
  git ls-remote https://github.com/humlab-speech/SPTK.git   1b3a5031d54b611df3390270522f86c4d85f987e
  git ls-remote https://github.com/humlab-speech/tandem.git a1fd952ff4823465fbf99e01a70042c8fc809c68
  ```
  Both must print the SHA. Then re-run the workflows (`gh run rerun 34812076499` for the
  head, or push a follow-up commit).
- [ ] **Task 1.4 — make the checkout self-describing.** Add the fork remote to the
  submodule config used locally (or at least document it in the planning record), so the
  next person does not have to rediscover that `origin` is upstream.

**Verification:** `actions/checkout` completes on all four workflows; the jobs reach
`setup-r-dependencies`.

### Phase 2 — fix the `pladdrr` dependency solve

- [ ] **Task 2.1 — install the GitHub-only `Suggests` explicitly.** In the workflows that
  solve dependencies (`R-CMD-check.yaml`, `test-coverage.yaml`, `pkgdown.yaml`, and
  `lintr.yml` if it solves them), add to `setup-r-dependencies`:
  ```yaml
          extra-packages: |
            github::humlab-speech/pladdrr
  ```
  pak understands `github::` specs, so the package comes from the fork while `DESCRIPTION`
  stays CRAN-clean. Keep the existing `any::` entries.
- [ ] **Task 2.2 — decide the fallback now, not later.** If the GitHub install proves
  flaky or too slow on some runner (the 09-11 green run took 48 min), the fallback is
  `dependencies: '"hard"'` plus the `any::` list, which installs no `Suggests` and lets
  the `pladdrr` tests skip through the existing `skip_without_pladdrr()` helper. Loss:
  the `pladdrr` code paths stop being exercised in CI. Prefer 2.1; record the fallback in
  the workflow comment.
- [ ] **Task 2.3 — confirm the skip path is the only `pladdrr` dependency.** The vignette
  chunks and tests already gate on `requireNamespace("pladdrr")`; verify with a
  `_R_CHECK_FORCE_SUGGESTS_=false` run that a missing `pladdrr` still yields
  `0 failures`, so Phase 2 cannot make CI dependent on the fork's availability.

**Verification:** `setup-r-dependencies` succeeds on ubuntu/macos/windows; the run
progresses to `check-r-package`.

### Phase 3 — remove the Windows format warnings (CRAN-facing)

- [ ] **Task 3.1 — replace `SMILE_CONSOLE_FORMAT`** in
  `src/opensmile/src/include/smileutil/smileConsole.h` with R's `R_PRINTF_FORMAT` logic
  verbatim (the table in F3), citing `R_ext/Print.h` in the comment. Keep the attribute
  on `smile_console_printf` and `smile_console_error` so Linux/macOS keep full
  `-Wformat` checking.
- [ ] **Task 3.2 — rebuild and re-check locally** (`R CMD check --as-cran`, expect
  `0 errors | 0 warnings | 3 notes`, `compiled code ... OK`), and confirm the gate still
  passes: `Rscript tools/check_cran_symbols.R <installed lib>`.
- [ ] **Task 3.3 — re-submit to win-builder** and require `Status: 2 NOTEs` with
  `whether package ... can be installed ... OK`. This is the only place MinGW's
  `gnu_printf` semantics can be observed.

**Verification:** win-builder R-release and R-devel both report no WARNING.

### Phase 4 — re-tighten the CI warning policy

- [ ] **Task 4.1 — prune the allowlist** in `R-CMD-check.yaml`: drop
  `whether package .* can be installed` (Phase 3 makes it real again),
  `compiled code` (now `OK`), `pragmas in C/C++ headers and code` and
  `compilation flags (in Makevars|used)` (both fixed at the source in 3.0.0), and
  `code/documentation mismatches` (closed by the S7 formals change). Keep
  `GNU extensions in Makefiles` and `line endings in Makefiles` only if the CI log still
  shows them on Windows, and shorten the comment block to the two live justifications —
  the stale paragraphs about `-Wno-register` and vendored `printf` symbols actively
  mislead.
- [ ] **Task 4.2 — prove the gate bites.** After 4.1, the next Windows run must be green
  *with* the install check enforced; if any pruned entry reappears, that is a genuine
  regression and belongs in a new allowance with a fresh justification, not a revert of
  the pruning.

**Verification:** `Fail on unexpected R CMD check WARNINGs` passes with the reduced list on
all matrix cells.

### Phase 5 — ship

- [ ] **Task 5.1 — green CI** on every workflow the push triggers: `gh run list --repo
  humlab-speech/superassp --limit 10` (expect `R-CMD-check`, `test-coverage`, `lintr`,
  `pkgdown`). Separately, decide F6: point `pkgcheck.yaml` at `master` and run it once via
  `gh workflow run pkgcheck.yaml`, or delete the workflow if `pkgcheck` is not wanted —
  a workflow that can never fire is worse than none, because it looks like coverage.
- [ ] **Task 5.2 — update the docs** (`cran-comments.md`: win-builder now checked, status
  `0 errors | 0 warnings | 2 notes`; `NEWS.md` if the format-attribute fix is
  user-visible — it is not; planning record with the push SHAs and CI run URLs).
- [ ] **Task 5.3 — re-submit to CRAN** once win-builder is clean and CI is green.

## Risks

| Risk | Mitigation |
|---|---|
| Pushing to the forks is a shared-state write the maintainer may want to time themselves | Phase 1 is two fast-forward pushes of commits that were already made; nothing is rewritten, and Task 1.3 verifies from the server side |
| `pladdrr` GitHub install is slow or flaky, turning CI into a coin flip | Task 2.2's `dependencies: '"hard"'` fallback keeps CI deterministic; the tests already skip |
| Windows CI's historical silent test crash (the reason `OPENSMILE_TRACE` and the flushing reporter exist) resurfaces once checkout works | Win-builder's Windows test run passed in 11 min, so the package is sound on Windows; if CI still dies, the in-tree diagnostics are already in place to localise it |
| Pruning the allowlist turns a previously-tolerated warning into a red build | That is the intent; do it in the same change as Phase 3 so the install check is genuinely clean first |
| Submodule pins drift again on the next fork commit | The 2026-09-13 record already lists this as a residual; a CI preflight step (`git submodule status --recursive` + `git ls-remote` reachability) would catch it before a push, if it is worth the workflow complexity |
