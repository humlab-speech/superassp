# Vendored Upstream Sources

`superassp` vendors six upstream C/C++ libraries under `src/` as ordinary
tracked files. They are **not** git submodules: there is no
`git submodule update --init` step, a plain `git clone` (or an unpacked source
tarball) already contains every file the build needs, and `R CMD build` sees the
same tree as a checkout.

Every vendored tree is byte-identical to the upstream revision listed below and
keeps its own `.gitignore`, so local build output (`*.o`, `*.so`,
`src/SPTK/build/`, `src/SPTK/bin/`, …) stays untracked. Together the six trees
are 2031 files / ~18 MB.

## Pinned revisions

| Path | Source repository | Pinned revision | Commit date | Tracked files |
|------|-------------------|-----------------|-------------|---------------|
| `src/SPTK` | https://github.com/humlab-speech/SPTK.git (branch `superassp-pin`) | `1b3a5031d54b611df3390270522f86c4d85f987e` | 2026-09-13 | 735 |
| `src/ESTK` | https://github.com/festvox/speech_tools | `9208f35ac1a8cd32c949b069c3ef073364f18a8b` | 2021-10-20 | 971 |
| `src/tcl-snack` | https://github.com/humlab-speech/snack.git | `2625509db5b271c76401eeb402237966f0605556` | 2026-04-16 | 213 |
| `src/Yin-Pitch-Tracking` | https://github.com/ashokfernandez/Yin-Pitch-Tracking.git | `69483b048bea0faac73a49e209577aebeb5e9680` | 2014-01-14 | 11 |
| `src/pyin` | https://github.com/aguai/pyin.git | `ae1df29ea5948002726cdf59e50b94d896084273` | 2018-04-16 | 47 |
| `src/tandem` | https://github.com/humlab-speech/tandem.git | `a1fd952ff4823465fbf99e01a70042c8fc809c68` | 2026-09-13 | 54 |

Notes:

* `src/SPTK`, `src/tcl-snack` and `src/tandem` are `humlab-speech` forks that
  carry the R-integration patches this package depends on (for example
  diagnostics routed through `Rprintf` instead of `stderr`/`stdout`). `src/SPTK`
  is pinned on the `superassp-pin` branch, **not** on `master` — do not
  "update" it to the fork's `master` head, which is unfixed upstream `v4.0`.
* The `src/ESTK` revision `9208f35` exists in `festvox/speech_tools` but is no
  longer reachable from that repository's current `master` head; keep the SHA.
* Only `src/SPTK` and `src/tandem` are compiled into the package (source lists
  and `-I` flags in `src/Makevars.in`). `src/ESTK`, `src/tcl-snack`,
  `src/Yin-Pitch-Tracking` and `src/pyin` are kept as in-tree upstream sources;
  `src/ESTK` and `src/tcl-snack` are additionally excluded from the built
  package by `.Rbuildignore`.
* Each tree ships its own upstream `LICENSE`/`COPYING`; those files are retained
  verbatim.

## Updating a vendored tree

Patch the change in the fork (or upstream) first — editing the vendored copy in
place means the next re-vendor silently drops the patch.

```sh
# /tmp/upstream is a clone of the repository in the table above
git -C /tmp/upstream fetch origin && git -C /tmp/upstream log --oneline -1 <new-revision>
git -C /tmp/upstream archive <new-revision> | tar -x -C src/SPTK
git status src/SPTK            # untracked lines are files the revision added
git add -A src/SPTK
# Ignore rules come from this repository *and* from each tree's own .gitignore,
# so confirm that nothing but build output was skipped. Two tracked files
# (src/tcl-snack/unix/pkgIndex.tcl.dll, src/tandem/tandem_128/.tandem.cpp.swp)
# are matched by those rules; they are tracked, so they only need `git add -f`
# if a revision rewrites them as new files.
git status --ignored --short src/SPTK | grep '^!!'
```

Then update the revision, date and file count in the table above and test.
Vendored DSP code is compiled into `superassp`, so a revision bump can change
the numeric output of `trk_*` / `lst_*` functions.

## Checking a vendored tree

The file counts in the table double as a integrity check — a mismatch means the
tree drifted from the recorded revision:

```sh
git ls-files src/SPTK | wc -l   # 735
```
