# SSFF read/write cross-validation tests against wrassp — 2026-09-19

**Status:** plan. Nothing here is implemented yet.

## 1. Mandate

Extend the test suite so that the SSFF read/write facilities added in `fa289e5`
(`planning/2026-09-19-ssff-read-performance.md`) are pinned to an **external**
reference rather than to our own past output. `wrassp` is that reference: it is
the package the ASSP reader/writer came from, it ships the same file format
through a different audio-loading path and different function names, and it is on
CRAN and already in `Suggests`.

The clause to satisfy is "still behave the same as before", which has three
distinct meanings and therefore three kinds of test:

| meaning | how it is tested |
|---|---|
| our reader agrees with an independent implementation of SSFF | frozen gold files produced by `wrassp`, compared against expectations produced by `wrassp::read.AsspDataObj` |
| our writer produces the bytes the format (and `wrassp`) expects | `write_ssff(read_ssff(gold))` byte-identical to `gold`; `write_ssff(obj)` byte-identical to `wrassp::write.AsspDataObj(obj)` |
| the deliberate 3.1.0 changes are the *only* behavioural changes | explicit assertions that pin `read_track()`'s `0 -> NA` and the writer's `NA/NaN -> 0`, and that `read_ssff()`/the legacy alias stay verbatim like `wrassp` |

## 2. What already exists (coverage inventory)

| file | covers |
|---|---|
| `tests/testthat/test-ssff-read-engine.R` | 63 assertions against an in-test `readBin()` decode: all storage types, interleaved layout, big-endian swap, windowed reads, `tracks=`, `threads=`, zero/NA rules, 300-generic header |
| `test-read_ssff.R`, `test-write_ssff.R`, `test-io-roundtrip.R`, `test-edge-cases.R` | smoke-level read/write, JSTF dispatch, malformed-file rejection |
| `test_wrassp.R` (403 assertions) | *superassp's own* `trk_*` wrappers, not the file format |
| `testthat/golden/voiceanalysis/a1.rds` | the established golden-fixture convention (`testthat::test_path("golden", ...)`) |
| `tests/signalfiles/` | third-party signal files (ESPS, wav) |

Gaps this plan closes: no cross-package check, no externally produced SSFF
fixture, writer bytes not pinned, window/selection/thread semantics not compared
against another implementation, no recorded fixture provenance.

## 3. Scoping probe — what `wrassp` gives us

Measured 2026-09-19: `wrassp` **1.0.6**, R 4.6.1, superassp work tree at `fa289e5`,
input `inst/samples/sustained/a1.wav` (4.035 s, md5 `c2afca3ef1496fa052874caa54f0b2bb`).

### F1 — real analysis files read identically by both packages

| producer | file | shape | `read.AsspDataObj` vs `read_ssff` |
|---|---|---|---|
| `ksvF0` | `a1.f0w` | 805 × 1 FLOAT | equal |
| `mhsF0` | `a1.f0m` | 805 × 1 FLOAT | equal |
| `rmsana` | `a1.rms` | 805 × 1 FLOAT | equal |
| `zcrana` | `a1.zcr` | 805 × 1 FLOAT | equal |
| `acfana` | `a1.acf` | 805 × 48 DOUBLE | equal |
| `cepstrum` | `a1.cep` | 805 × 1025 FLOAT | equal |
| `dftSpectrum` | `a1.dft` | 805 × 1025 FLOAT | equal |
| `forest` | `a1.fms` | 2 tracks × SHORT × 4 | equal |

"equal" = values, dims, types, all five carried attributes (`sampleRate`,
`startTime`, `startRecord`, `endRecord`, `trackFormats`) identical without
tolerance. Headers are identical too (`Machine IBM-PC`, `Record_Freq`,
`Start_Time`, `Column`, `Original_Freq`).

### F2 — the two writers emit identical bytes

For synthetic objects (REAL32 1-field, REAL64 4-field, INT16 1-field, mixed
REAL32+REAL64 multi-track), `write_ssff()` and `wrassp::write.AsspDataObj()`
produce **byte-identical** files (274/410/213/530 bytes).

### F3 — read → write on a wrassp file is byte-identical

| gold file | shape | `write_ssff(read_ssff(gold))` == `gold` |
|---|---|---|
| `a1.f0w` | 805 × 1 FLOAT (33.8 % zeros) | yes (3373 B) |
| windowed `dftSpectrum` | 15 × 1025 FLOAT | yes (61 635 B) |
| windowed `acfana` | 60 × 48 DOUBLE | yes (23 196 B) |
| windowed `forest` | 2 tracks × (60 × 4) SHORT | yes (1131 B) |

This is the strongest writer assertion available and directly tests the
"all information retained in order to safely rewrite without data loss" design
claim in `src/dataobj.c`.

### F4 — windowed reads agree with `wrassp`

On `a1.f0w` (`startTime 0.0025`, 199.55 Hz): `begin=0.5,end=1.5` → both return
records 101–300; `begin=2.0,end=3.0` → both return 401–599; values and *all*
attributes equal.

### F5 — the two deliberate deviations, measured

| behaviour | `wrassp` | superassp 3.1.0 |
|---|---|---|
| `NA`/`NaN` in a track on write | stored as NaN bit pattern | stored as `0` |
| reading that file | NaN comes back as NaN | NaN comes back as NaN (not NA) |
| plain reader (`read_ssff`, legacy alias) | verbatim | verbatim (same) |
| `read_track()` | n/a (`wrassp` has no equivalent) | stored `0` → `NA` (272 of 805 records in `a1.f0w` are `NA`) |
| audio tracks | verbatim | verbatim (never masked) |

Both packages read each other's NA files without error.

### F6 — loading `wrassp` in a superassp session conflicts with our S3 methods

`loadNamespace("wrassp")` prints *"Registered S3 methods overwritten by 'wrassp'"*
for `as_tibble.AsspDataObj` and `print.AsspDataObj`. Under `pkgload::load_all()`
dispatch still resolved to superassp's methods, but the registry conflict is real
and its outcome may differ in an installed-package session. Any test that loads
`wrassp` must therefore not share a session with tests that depend on those
methods.

## 4. Test design

### T1 — golden reader parity (no `wrassp` needed at test time)

New `tests/testthat/golden/ssff/` with, per fixture: the file (`<name>.ssff`) and
the expectation captured from `wrassp::read.AsspDataObj()` at generation time
(`<name>.expected.rds`), plus `PROVENANCE.md` (commands, versions, md5s).

Assertions per fixture: `read_ssff()` equals the expectation on values (tolerance
0), dims, `typeof`, names, and the five carried attributes; `trackFormats` equal;
for the mixed/shape fixtures also assert the storage class per track (REALSXP vs
INTSXP) so a future kernel change cannot silently switch destination types.

### T2 — live cross-package parity (skipped without `wrassp`)

A subprocess (`Rscript tests/testthat/scripts/wrassp_interop.R`, §6) regenerates
the analysis files with `wrassp` from `inst/samples/sustained/a1.wav`, reads them
with `wrassp::read.AsspDataObj()` and `read_ssff()`, and returns a comparison
tibble. The test asserts every row is "equal" and fails with the offending shape
on mismatch. Covers F1 and F4 (windows: full, mid-file, `begin == end`,
`begin > end` error parity).

### T3 — writer parity and byte-frozen round trips

* Live: `write_ssff(obj)` bytes == `wrassp::write.AsspDataObj(obj)` bytes for the
  F2 shapes (skipped without `wrassp`).
* Always: `write_ssff(read_ssff(gold))` == `gold` bytes for every NA-free fixture
  (F3). This is the regression fence that will catch any writer edit that changes
  formatting, endianness, record order or the NA policy.

### T4 — window and selection invariants

* `read_ssff(f, begin, end)` == `read.AsspDataObj(f, begin, end)` (live, T2).
* `read_ssff(f, tracks = t)` == the `t` columns of the full read, for single- and
  multi-track fixtures; `tracks = c()` and `tracks = NULL` read everything;
  unknown names abort with the file's track list (already covered synthetically,
  extended to a gold file).
* `threads = 4` identical to `threads = 1` on the 1025-field fixture (bounds the
  block split with a realistic record size).

### T5 — pinned deviations

* On a fixture written with `NA`/`NaN` (committed, produced by `write_ssff`):
  on-disk values are exactly `0`; `read_ssff()` returns `0`; `read_track()`
  returns `NA`; `wrassp::read.AsspDataObj()` also returns `0` (live, optional) —
  i.e. the file is wrassp-readable, which is the point of the change.
* On a fixture written by `wrassp` with `NA` (NaN on disk): our `read_ssff()`
  returns NaN, `read_track()` returns NaN (not `NA`), and neither errors.
* On the `ksvF0` gold fixture: `sum(is.na(read_track(f)$F0)) == 272` and
  `sum(read_ssff(f)$F0 == 0) == 272`, i.e. the mask hits exactly the stored zeros.
* Audio: `read_track(<wav>)` has no `NA` in the audio track, equal to `read_ssff`.

### T6 — legacy alias parity

`superassp::read.AsspDataObj()` / `getAsspDataObj()` (which deliberately keep the
verbatim semantics) must return what `wrassp::read.AsspDataObj()` returns, for
full and windowed reads (live). This keeps the compatibility surface honest now
that `read_track()` diverges by design.

### T7 — track manipulation round trip

`addTrack()`/`delTrack()` then `write_ssff()`/`read_ssff()`: names,
`trackFormats` order and values survive, and a re-written gold file stays
byte-identical when no track was touched. (`wrassp` exposes the same helpers;
the live variant can compare `delTrack` results between packages.)

### T8 — keep the synthetic edge suite

The cases `wrassp` cannot produce stay in `test-ssff-read-engine.R`: BYTE/LONG
storage types, byte-swapped big-endian fixtures, 300-generic headers, truncated
files, unknown `tracks=`, 0-record files, `samples = TRUE` windowing. Gold files
cover the *realistic* shapes; the synthetic suite covers the *format space*.

## 5. Fixture inventory

All generated from `inst/samples/sustained/a1.wav` unless stated. Sizes measured.

| fixture | produced by | shape | size |
|---|---|---|---|
| `f0_ksv.ssff` | `wrassp::ksvF0(wav, toFile = TRUE)` | 805 × 1 FLOAT, 272 zeros | 3.3 KB |
| `f0_mhs.ssff` | `wrassp::mhsF0` | 805 × 1 FLOAT | 3.3 KB |
| `rms.ssff` | `wrassp::rmsana` | 805 × 1 FLOAT | 3.3 KB |
| `spectrum_dft.ssff` | `wrassp::dftSpectrum(beginTime = 0, endTime = 0.3, windowShift = 20)` | 15 × 1025 FLOAT | 60.2 KB |
| `acf48.ssff` | `wrassp::acfana(beginTime = 0, endTime = 0.3)` | 60 × 48 DOUBLE | 23 KB |
| `forest.fms` | `wrassp::forest(beginTime = 0, endTime = 0.3)` | 2 tracks × (60 × 4) SHORT | 1.1 KB |
| `f0_ksv_be.ssff` | derived: `f0_ksv.ssff` with each 4-byte group reversed and `Machine SPARC  ` | 805 × 1 FLOAT big-endian | 3.3 KB |
| `na_written_by_wrassp.ssff` | `wrassp::write.AsspDataObj` with an `NA` track | 4 × 1 FLOAT (NaN) | < 1 KB |
| `na_written_by_superassp.ssff` | `write_ssff` with `NA`/`NaN` | 4 × 1 FLOAT (0) | < 1 KB |

Plus one `.expected.rds` per fixture (the `wrassp` read result). Total budget
≈ 110 KB, in line with the existing `golden/voiceanalysis/a1.rds`.

`PROVENANCE.md` records: wrassp version, R version, superassp commit, md5 of the
input wav, the exact generation calls, and the md5 of every fixture, so the set
can be audited and regenerated.

## 6. Mechanics

* **Generation script** `tests/testthat/golden/ssff/generate.R`, run manually
  (`Rscript tests/testthat/golden/ssff/generate.R`), refuses to overwrite unless
  `--force`, and writes `PROVENANCE.md`. It requires `wrassp`, so it is not part
  of `R CMD check`; the committed fixtures are.
* **Subprocess isolation** for every live `wrassp` comparison: a helper in
  `tests/testthat/helper-wrassp-subprocess.R` locates
  `file.path(R.home("bin"), "Rscript")`, runs `scripts/wrassp_interop.R` with
  `system2()`, passing input paths and an output RDS path, and the test reads the
  RDS back. One subprocess serves all live checks (start-up ~1 s). This keeps
  `wrassp`'s namespace out of the test session (F6), keeps `load_all()` and
  installed-package runs equivalent, and needs no new dependency (`callr` is
  avoidable).
* **Skip rules**: live tests `skip_if_not_installed("wrassp")`; golden tests run
  everywhere. The interop script itself must not fail the suite if a wrassp
  analysis function is unavailable — it returns a `status` column the test
  inspects.
* **No timing assertions** anywhere (CI variance); performance is tracked by the
  manual benchmark corpus in `planning/2026-09-19-ssff-read-performance.md`, not
  by the suite.
* **Portability**: paths through `normalizePath()`/`shQuote()` (Windows), no
  `system()` shell features beyond `system2()`, no `file.copy` of the gold files
  (read-only use).

## 7. Implementation order (each step independently green)

| # | Step | Gate |
|---|---|---|
| 1 | `generate.R` + fixtures + `PROVENANCE.md` | re-running the generator reproduces the same md5s; total size ≤ 150 KB |
| 2 | T1 golden parity tests | pass; deliberately corrupting one byte of a fixture makes them fail |
| 3 | T3 byte-frozen round trips | pass; the NA-free fixtures survive read→write byte-identically |
| 4 | T5 deviation pins | pass; `read_track` NA count on the F0 fixture is exactly 272 |
| 5 | T4 selection/threads invariants on gold files | pass |
| 6 | T2/T6/T7 live interop in a subprocess | pass with wrassp 1.0.6; skip cleanly when wrassp is absent (verify by unsetting the library path) |
| 7 | Update `CLAUDE.md` testing notes + `NEWS.md` (test-only entry) | `devtools::test()` green |

## 8. Risks and open decisions

* **Don't let goldens become self-fulfilling.** Expectations must come from
  `wrassp::read.AsspDataObj()` at generation time; regenerating them from
  superassp would defeat the purpose. `PROVENANCE.md` states the producer for
  each file.
* **`wrassp` is an oracle for the format, not for our intended deviations** (F5):
  the `0 -> NA` and `NA/NaN -> 0` rules are asserted as *documented divergences*,
  with the wrassp-written NaN fixture pinning that NaN is not masked.
* **Fixture staleness**: if wrassp ever changes its writer, the gold files stay
  valid (they are frozen); only the live T3/T2 checks would flag a difference,
  and the fix is to regenerate deliberately, never to relax the assertion.
* **Subprocess cost**: one `Rscript` start plus eight wrassp analyses
  (≈ 1–2 s total measured) — acceptable; keep it in one file so it can be
  skipped as a unit.
* **Open decision 1**: whether to also freeze a golden for `read_track()`'s
  NA-masked view (duplicative of T5 as written — currently no).
* **Open decision 2**: whether the synthetic BE fixture should be replaced by the
  `f0_ksv_be.ssff` derivative (currently both; the derivative is more realistic,
  the synthetic one exercises REAL64/BYTE swaps as well).

## 9. Implementation record — 2026-09-19

### 9.1 What landed

| item | implementation |
|---|---|
| Fixtures | `tests/testthat/golden/ssff/` (8 fixtures + 8 `.expected.rds`, `PROVENANCE.md`, `generate.R`), 228 KB total. |
| Generator | `generate.R` writes the set, captures the `wrassp::read.AsspDataObj()` expectations (xz-compressed, `filePath` normalised to the file name so they are reproducible), and supports `--force` / `--check`. `--check` was verified: re-running reproduces all 16 files. |
| Golden tests | `test-ssff-wrassp-golden.R`, 180 assertions: provenance/md5 integrity, read parity against every expectation, fixture shape/format inventory, byte-identical read→write round trips, big-endian twin equivalence, the zero→NA rule (272/805 records), wrassp's NaN encoding, `tracks=`/`threads=`, `delTrack()`/`addTrack()` round trips. Requires no wrassp. |
| Live interop | `test-ssff-wrassp-interop.R` (96 assertions) + `helper-wrassp-interop.R` + `scripts/wrassp_interop.R`. |
| Isolation | The subprocess produces the **wrassp side only** (files, `read.AsspDataObj()` results, windowed reads, `delTrack()` result, four writer objects + their bytes) and returns them as an RDS; the parent compares with the code under test. This is a refinement of §6: it also removes the "which superassp do we compare against" problem, since the subprocess never loads superassp at all. |

### 9.2 Gate results

| gate (§7) | result |
|---|---|
| 1 — regeneration reproduces the fixture set | `generate.R --check` → "16 fixtures match" |
| 2 — corrupting a fixture fails the tests | flipping one data byte in `rms.ssff` fails the provenance check *and* the value parity check (2 failures), fixture restored afterwards |
| 3 — NA-free fixtures survive read→write byte-identically | passing (5 fixtures: the two 3.3 KB F0 files, `rms.ssff`, `spectrum_dft.ssff`, `acf48.ssff`) |
| 4 — F0 fixture masks exactly 272 records | passing, `identical(is.na(masked), expected == 0)` |
| 5 — selection/threads invariants on gold files | passing |
| 6 — live interop passes; skips without wrassp | passing with wrassp 1.0.6; with the wrassp installation temporarily hidden the interop file reports 5 skips / 0 failures while the golden file still runs all 180 assertions |
| 7 — docs updated | `CLAUDE.md` testing section, `NEWS.md` (3.1.0 → "Tests") |
| 8 — whole suite green | `test_dir()`: **3235 passed, 0 failed, 0 errors, 40 skipped** (previous green baseline 2958 + 277 new assertions: 180 golden + 97 interop) |

### 9.3 Deviations and findings

* **Fixture budget**: 228 KB rather than the estimated ≤ 150 KB. The expectations
  dominate (`spectrum_dft.expected.rds` 53 KB, `acf48.expected.rds` 22 KB even
  xz-compressed); the fixtures themselves are 94 KB. Still small next to the
  16 MB tarball, so no trimming.
* **`addTrack()` for a *new* track is broken** (pre-existing, inherited from
  wrassp): `attr(dobj, 'trackFormats') <- append(...)` is discarded by R's
  copy-on-modify semantics, so the new format never registers and `write_ssff()`
  then fails with "Not enough format specifiers for the data tracks." Call sites
  work around it (`R/ssff_pladdrr_intensity.R:189` and three in
  `R/ssff_pladdrr_formant.R` append the format by hand). T7 therefore covers
  `delTrack()` and the `addTrack(..., deleteExisting = TRUE)` *replace* path,
  which do maintain the bookkeeping. Fixing the bug is a separate change: it
  would obsolete those workarounds in four files, three of them pladdrr-only and
  therefore not verifiable in this environment.
* **The big-endian fixture carries a wrassp expectation too** — wrassp reads the
  byte-swapped derivative correctly (verified while implementing), so the
  swapped path is compared against the reference implementation rather than
  against our own little-endian result alone.
* **`read.AsspDataObj()`/`getAsspDataObj()` are internal** (not exported, per the
  export policy), so the alias-parity test calls them unqualified from the test
  environment rather than via `::`.
* **F6 is real, and the subprocess requirement covers the availability check
  too.** The first version of `helper-wrassp-interop.R` guarded itself with
  `requireNamespace("wrassp", quietly = TRUE)`. Even without attaching the
  package, that loads wrassp's namespace and registers its
  `as_tibble.AsspDataObj()`/`print.AsspDataObj()` over superassp's, which broke
  four assertions in `test-track-naming-phase2.R` (the tibble view lost
  `track_labels`/`track_descriptions` and the `*_Hz` column names) — reproduced
  in isolation as: `tibble::as_tibble(trk_formant_forest(...))` before and after
  `requireNamespace("wrassp")` are not identical. The helper now checks
  `nzchar(system.file(package = "wrassp"))`, which loads nothing; the test file
  asserts the interop run leaves `"wrassp" %in% loadedNamespaces()` FALSE.
