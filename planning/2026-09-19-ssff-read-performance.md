# SSFF read performance and the 0-is-missing convention — 2026-09-19

**Status:** implemented (2026-09-19) — see §8 for what landed, the measured
before/after numbers and the deviations from this plan. The design sections below
are kept as written so the decisions can be reviewed against the implementation.

Two asks, one plan because they touch the same code path:

1. **Read SSFF tracks faster**, using the ideas that make `simdjson`/`RcppSimdJson` fast
   (already the reader for JSTF in `R/jstf_io.R`) as a design reference.
2. **`read_track()` must return `NA` for stored values of exactly `0` in SSFF tracks.**
   SSFF has no NULL/NA encoding; `0` is this package's substitute (documented in
   `as_tibble.AsspDataObj(na.zeros = TRUE)`, `R/ssff_c_assp_formant_forest.R:54`
   "0 = missing", and the praatsauce unvoiced-F0 convention recorded in `CLAUDE.md`).

## 0. Method and artifacts

| Artifact | Provenance |
|---|---|
| `read_ssff` timings (tables 2.1–2.3) | installed `superassp` 3.0.0, R 4.6.1, Apple M1 Pro (10 cores), macOS 25.6, best/median of 3–7 reps, files in page cache |
| Synthetic corpus | 7 SSFF + 2 WAV files generated in a scratch dir; shapes listed in table 2.1 |
| `sample` profile (2.2) | `sample <pid> 12` on an R loop calling `read_ssff("multispec.ssff")` continuously for 90 s |
| Kernel microbenchmarks (2.4) | standalone C++ (`clang++ -O3 -std=c++17 -fno-fast-math`, Apple clang 21), warm buffers, best of 5 |
| Fused-fill experiment (2.5) | ad-hoc `R CMD SHLIB` kernel (`mmap` + single pass + fused NA mask), driven from R, median of 7, correctness-checked against `read_ssff` incl. the NA mask |
| xsimd assessment (§3.1) | installed `RcppXsimd` 7.1.6-2 (ships xsimd 7.1.3; `XSIMD_VERSION_*` printed at runtime), compile probes for `batch_cast`, `clang++ -S` instruction-mix comparison, upstream `xtensor-stack/xsimd` `master` (14.3) header check |
| R transpose reference (§3.1) | `t(matrix(rnorm(20000*1025), 20000, 1025))` median of 5 on the same machine |
| `RcppSimdJson` reference | `RcppSimdJson` 7.1.6.2 (Suggests; used today via `RcppSimdJson::fload()` in `R/jstf_io.R:139`) |

Harness: [Appendix A](#appendix-a--reproduction). Numbers are single-machine and should be
treated as ratios, not absolutes.

## 1. Where SSFF reading happens today

```
read_track(file)                        R/jstf_io.R:207   dispatch on extension
  └── read_jstf()                        R/jstf_io.R:85   JSTF (JSON) — RcppSimdJson:131 / jsonlite:155
  └── read_ssff(file, begin, end, samples) R/read_ssff.R:34
        └── .External("getDObj2", …)     R/assp_dataobj_methods.R:22 is the legacy alias
              └── getDObj2()             src/dataobj.c:48
                    ├── asspFOpen()      src/assp/asspfio.c:131
                    │     └── getSSFFhdr()      src/assp/headers.c:2568  (text header, fgetl)
                    ├── allocDataBuf()   src/assp/dataobj.c:828         (calloc)
                    ├── asspFFill()      src/assp/asspfio.c:777         (one fread + swapDataBuf)
                    │     ├── asspFSeek()/asspFRead()  src/assp/asspfio.c:321,408
                    │     └── swapDataBuf()            src/assp/dataobj.c:943 (per-record swap pass)
                    └── dobj2AsspDataObj()   src/dataobj.c:154
                          └── getDObjTrackData() per descriptor   src/dataobj.c:391
                                 memcpy(record → tempBuffer)      src/dataobj.c:451
                                 switch(desc->format) inside the record loop
                                 scalar per-field loops           src/dataobj.c:474-536
```

Write path (relevant to §4.9): `write_ssff()` → `writeDObj_` → `sexp2dobj()` (`src/dataobj.c:594`)
→ `addTrackData()` (`src/dataobj.c:934`, `(float) numPtr[...]` / `f64Ptr[n] = numPtr[...]`)
→ `asspFFlush()` (`src/assp/asspfio.c:843`), which swaps on the way out when the file
endianness differs from the host.

Facts that shape the design:

* **The data section is a fixed-layout record array.** Record `r` of descriptor `d` lives at
  `headerSize + r*recordSize + d->offset`, `d->numFields` values of `dd->format` each
  (`setRecordSize()`, `src/assp/dataobj.c:607`). There is **no structure to discover in the
  data** — the "parse" is offset arithmetic.
* **Binary only** (`dop->fileData = FDF_BIN`, `src/assp/headers.c:2592`), so the ASCII branch
  of `asspFPrint` is irrelevant to SSFF.
* **Byte order comes from the header's `Machine` line** (`src/assp/headers.c:2621`):
  `IBM-PC`/`VAX` = little-endian (what this package writes), anything else (`SPARC`, …) =
  big-endian and triggers `swapDataBuf()`.
* **Windowed reads already seek** (`getDObj2`, `src/dataobj.c:121-123`), so `begin`/`end`
  only read the window. That part is fine and must stay.
* Every reader call materialises **all** descriptors; there is no track selection anywhere in
  the API.
* `getDObjTrackData()` copies each record into a `malloc`'d `tempBuffer` before touching it
  (`src/dataobj.c:451`) — pure overhead on the read path (the temp buffer exists for the
  legacy in-place-swap workflow, not for reading).

## 2. Baseline measurements

### 2.1 Throughput of `read_ssff` (installed 3.0.0)

| workload | shape | size | time | throughput |
|---|---|---|---|---|
| 1-track spectrum, REAL32 | 20 000 rec × 1025 fields | 82.0 MB | 0.036–0.047 s | 1.8–2.3 GB/s |
| 4-track spectrum, REAL32 | 5 000 rec × 1025 × 4 | 82.0 MB | 0.037–0.059 s | 1.4–2.2 GB/s |
| 1-track spectrum, REAL64 | 10 000 rec × 513 | 41.0 MB | 0.010–0.011 s | 3.7–4.1 GB/s |
| 1-track spectrum, REAL32, big-endian | as row 1, `Machine SPARC` | 82.0 MB | 0.055–0.061 s | 1.4–1.5 GB/s |
| 10 tracks × 1 field, REAL32 | 100 000 rec | 4.0 MB | 0.006–0.007 s | 0.6–0.7 GB/s |
| audio 600 s mono, INT16 (WAV via `read_ssff`) | 26.5 M samples | 52.9 MB | 0.083–0.085 s | 0.62 GB/s |
| audio 60 s stereo, INT16 | 2.65 M × 2 | 10.6 MB | 0.015 s | 0.70 GB/s |
| tiny f0 track | 400 rec × 1 | 1.7 KB | 38–47 µs/call | — |

Reference points on the same machine: `readBin(<same file>, "raw")` = 2.9–4.6 GB/s;
pure-R push of the data through `readBin` + `matrix()` = 230 MB/s (0.36 s for row 1).

### 2.2 Profile of a repeated read (4-track spectrum)

`sample` on a tight `read_ssff` loop, top-of-stack collapsed:

| symbol | samples | share |
|---|---|---|
| `RunGenCollect` (R GC, triggered by allocating 4 × 41 MB matrices/iteration) | 4436 | 54 % |
| `getDObjTrackData` (scalar convert loops, leaf) | 2241 | 27 % |
| `_platform_memmove` (the per-record temp `memcpy`, 4 tracks × 82 MB) | 1141 | 14 % |
| `__read_nocancel` (the actual `fread`) | 192 | 2.3 % |
| `__bzero` (`calloc` of the data buffer) | 122 | 1.5 % |

**The bytes are not the problem; the copies and the allocations are.** The `fread` of 82 MB
from page cache is ~1.5 ms; everything else is ≥ 25 ms.

### 2.3 Per-call overhead and the header path

`read_ssff` on a 400-record f0 file costs **38.5 µs**, and the cost scales with the number of
header lines, not with their length:

| filler header lines | file size | µs/read |
|---|---|---|
| 0 | 1 731 B | 38.5 |
| 50 | 2 971 B | 76.0 |
| 500 | 14 621 B | 509 |
| 2000 (13-byte lines) | 36 621 B | 3145 |
| 2000 (27-byte lines) | 54 621 B | 3447 |

≈ **1.6–1.7 µs per header line**. Cause: `fgetl()` reads one character per `fgetc()`
(`src/assp/fgetl.c:80`), and each generic variable appends by walking the whole linked list
(`src/assp/headers.c:2727` → `addTSSFF_Generic()`, `src/assp/dataobj.c:461-476`, "search last
element in list"), i.e. O(n²) in the number of header lines — which is why the per-line cost
rises from 1.0 µs at 500 lines to 1.7 µs at 2000 lines. Typical headers (5–30 lines) therefore cost 10–50 µs — the
dominant part of the per-call overhead. `getFileSize()` itself is a seek, not a scan
(`src/assp/headers.c:771`).

### 2.4 Kernel microbenchmarks (warm buffers, standalone C++)

| kernel | time | rate |
|---|---|---|
| f32→f64 contiguous, + zero→NaN mask (scalar C, auto-vectorised) | 5.2 ms / 26.5 M | 20.5 GB/s |
| same, hand-written NEON (`vcvt_f64_f32` + `vbslq`) | 6.2 ms | 16.9 GB/s |
| i16→i32 contiguous, scalar (auto-vectorised) | 2.6 ms | 20.4 GB/s |
| i16→i32, hand NEON (`vmovl_s16`) | 2.7 ms | 19.7 GB/s |
| i16→i32, scalar on 10 threads | 1.3 ms | 39.5 GB/s |
| f32 big-endian → f64 native, scalar `bswap` | 5.2 ms | 20.5 GB/s |
| 1025-field record-major convert, current shape (incl. temp `memcpy`) | 5.4 ms/track | 15.3 GB/s |
| same, without the `memcpy` | 4.2 ms/track | 19.4 GB/s |
| same, NEON 4-field convert + mask | 4.3 ms/track | 18.9 GB/s |

Conclusions: (a) LLVM already auto-vectorises the contiguous loops to within noise of hand
NEON, so **hand-written SIMD buys nothing on Apple silicon for these kernels** and adds ISA
maintenance; (b) the multi-field (spectra) case is **store-pattern bound** — converting four
fields at a time does not help because the destination stride is `nRecord`; (c) the current
per-track scalar loop is within ~1.25× of a stripped-down loop on warm data.

### 2.5 Fused fill experiment (what a rewrite can reach)

Ad-hoc kernel: `mmap` the file, one pass over records, hoist the format switch, write directly
into pre-allocated R matrices, fuse byte swap and the zero→NA mask, optional record
blocking (`b`) and worker threads. Median of 7; correctness-checked against `read_ssff`.

| case | `read_ssff` | fill, 1 thread | fill, 8 threads |
|---|---|---|---|
| 1-trk spectrum f32 | 0.047 s | **0.029 s** (b=64) | **0.008 s** |
| 4-trk spectrum f32 | 0.037 s | **0.028 s** (b=512) | **0.009 s** |
| 1-trk spectrum f64 | 0.011 s | **0.009 s** (b=64) | **0.003 s** |
| 1-trk spectrum f32 BE | 0.061 s | **0.029 s** (b=64) | **0.008 s** |
| 10-trk 1-field f32 | 0.007 s | 0.001 s | 0.001 s |

Block size is load-bearing (same kernel, 1 thread): 1-track spectrum — no blocking 0.050 s,
b=64 0.029 s, b=512 0.035 s; 4-track spectrum — no blocking 0.074 s, b=64 0.036 s,
b=512 0.028 s. Hence the `RB ≈ 2 MB / recordSize` rule in §4.1.

→ **3.7–7.6× on the C-side work.** Caveat, measured: with a *fresh* R allocation per call,
the end-to-end gain shrinks to **1.2–2.3×** because page faults on the (unavoidable) result
matrices and R's GC dominate; only the 10-track case keeps its 2.3×. The rewrite removes work
the reader owns; it cannot remove the cost of the matrix R must hand back.

## 3. What simdjson actually buys, and which parts transfer

| `simdjson` technique (RcppSimdJson exposes it via `fload`/`fparse`) | SSFF analogue | verdict |
|---|---|---|
| **Stage 1** structural indexing: classify 64-byte blocks with `pshufb`/`cmpeq` + `movemask`, branch on masks instead of bytes | The data section has no structure to find — offsets are closed-form (`recordSize`, `dd->offset`) | **N/A** — go straight to stage 2. Do not build a tokenizer for binary records |
| **Stage 2** tape building off the index; per-value work is a table-driven load | Per-value work is a typed load from a computed address → same idea, implemented as per-format kernels | **Transfers** (D1) |
| **On Demand / `skip_child`**: don't materialise what wasn't asked for; `raw_json_token` passes bytes through | Skip descriptors the caller did not request — never touch their bytes, never allocate their matrices | **Transfers, highest value** (D4) |
| Whole-buffer read into one padded buffer; no per-token I/O | One `mmap` (or one `fread`) per file; drop the intermediate DOBJ buffer on the read path | **Transfers** (D3) |
| One allocation, linear writes, no per-value copies; no per-value `strdup`/temp buffers | Allocate the result once; fold the byte swap and the NA mask into the same pass; delete the per-record `memcpy` | **Transfers** (D2, D5) |
| Fused validation on the same blocks that are being parsed (UTF-8 during stage 1) | Fuse endian swap + widening + 0→NA masking into the conversion pass | **Transfers** (D5, D9) |
| Cache-sized blocks (keeps working set in L1/L2) | Record-blocked loops; block size chosen so source block + destination columns stay in L2 | **Transfers** (D1) |
| `parse_many`, thread-parallel stages for big inputs | Block-parallel fill of disjoint record ranges for large inputs | **Transfers, gated** (D7) |
| Branch-free token iteration; format specialised at compile time | Hoist `switch(dd->format)` out of the record loop → one specialised kernel per `DF_*` | **Transfers** (D1) |
| SIMD string scanning (`fgetl`-style line reading) | Only the **text header** is text: replace char-by-char `fgetl` with one bulk read + `memchr` line scan | **Transfers, small but cheap** (D6) |

Non-goals copied from this repo's policy: no `-march=native`, xsimd baseline ISA only
(`CLAUDE.md` → SIMD), and faithful numeric results (double precision, same values as today).

### 3.1 Would xsimd / RcppXsimd help the reader? — measured, no

xsimd is already a build dependency (`LinkingTo: RcppXsimd`, `-DRCPPXSIMD_AVAILABLE` in
`src/Makevars.in`), so adopting it in the reader costs no new dependency. It still is not the
right tool here, for three measured/verified reasons:

1. **The available version cannot express the kernels.** RcppXsimd 7.1.6-2 (CRAN,
   2026-01-22) ships **xsimd 7.1.3**; upstream is 14.3. In 7.1.3 the only batch conversions
   are same-width float↔int32 and double↔int64 (`xsimd/types/xsimd_neon_conversion.hpp:27-32`;
   a `batch_cast` compile probe fails with "no member named 'batch_cast'"). The reader's two
   primary kernels — `float`→`double` and `int16`→`int` **widening** — are therefore
   inexpressible, as is a byte-shuffle for the big-endian swap (only `bitwise_cast` exists).
   Upstream 14.3 has `batch_cast` (`include/xsimd/types/xsimd_api.hpp:262`); getting it means
   vendoring a newer header-only xsimd or waiting on RcppXsimd.
2. **Where it *is* expressible, it is slower than the compiler.** Same machine, same flags,
   26.5 M elements: `f64` load → compare-to-zero → select NA → store as plain C++ that LLVM
   auto-vectorizes = **5.85 ms**; the xsimd version of the same loop = **8.36 ms** (43 %
   slower); `memcpy` = 4.66 ms. `clang++ -S` shows both bodies contain the same NEON
   compare/select sequence (`fcmeq.2d` + `bit.16b`), i.e. identical ISA work — the gap is loop
   structure and memory behaviour, not vector width. xsimd 7's NEON backend has no advantage
   over LLVM's vectorizer, and the loop is memory-bound anyway.
3. **The remaining cost is not arithmetic.** With the destination pre-allocated, the plan's
   fused fill takes 0.029 s on the 20 000 × 1025 spectrum shape; R's own cache-blocked `t()`
   on the identical shape (pure data movement, zero arithmetic) takes 0.039 s. The multi-field
   case is bound by the transpose/store pattern, so no ISA-level change can pay for itself
   there.

**Decision:** the new reader module stays plain portable C++ (per-format kernels, unit-stride
inner loops) and relies on the compiler's auto-vectorizer; xsimd stays where it already is —
the compute-bound arithmetic primitives in `src/simd_utils.hpp` (`simd_dot`, `simd_energy`,
`simd_fir`), which are the kernels where SIMD actually pays. Revisit only if (a) the rewritten
reader profiles the multi-field convert as the top cost, and (b) a newer xsimd (≥ 8) is
available; even then the realistic ceiling is the memory traffic, not the ISA.

Caveat on generality: the timings above are NEON (Apple silicon) [INFERENCE for x86-64, not
measured here]. The two decisive facts are platform-independent: the 7.1.3 conversion API is
the same across its SSE/AVX/NEON backends (`to_int`/`to_float` are same-width only,
`xsimd_sse_conversion.hpp` and `xsimd_neon_conversion.hpp`), and the kernels are bound by
memory traffic and the store pattern rather than by vector width.

## 4. Design

### 4.1 D1 — per-format conversion kernels, one per shape

New module `src/ssff_convert.hpp` / `src/ssff_convert.cpp` (registered in
`src/Makevars.in`'s `C_SOURCES`/`CXX_SOURCES`; there is no separate `Makevars.win` — both
`./configure` and `configure.win` regenerate `src/Makevars` from `Makevars.in`). Public entry:

```c
/* Fill R matrices (REALSXP for float tracks, INTSXP for integer tracks) for the
   descriptor list [first, last) from the record buffer. */
int ssff_fill_tracks(const void *records, long nrec, size_t recordSize,
                     DDESC *first, DDESC *last, SEXP *dest, int zeroToNa);
```

One kernel per (`dd->format` × shape), selected once per track — never per record:

* **contiguous** (`dd->numFields == 1 && recordSize == elemsize`, the case for every
  `trk_*` single-value track and all audio but the last channel): straight widening loop
  (`float`→`double`, `int16`→`int`, …) that LLVM vectorises to memory bandwidth (20.4 GB/s
  measured, §2.4).
* **strided single field** (`numFields == 1`, `recordSize > elemsize`): same loop with a
  runtime stride.
* **multi-field** (spectra): loop nest *record block → descriptor → field → record in block*
  (unit-stride destination stores, strided source loads). The order is not cosmetic —
  measured on the 4-track spectrum file (§2.5): no blocking 0.074 s → block 64 recs 0.036 s →
  block 512 recs 0.028 s; on the 1-track file block 64 recs 0.029 s beat block 512 recs
  0.035 s. Size the block for ~1–2 MB of source (`RB = clamp(2 MB / recordSize, 16, 4096)`),
  which puts both files near their measured optimum. The current reader's per-record
  `memcpy` + field-inner nest is what this replaces.
* **REAL64 multi-field, native endian, no mask** → `memcpy` the whole track when
  `recordSize == 8*numFields` (single-descriptor file, the common `f64` spectrum case).
* integer destinations stay `INTSXP` exactly as today (`getDObjTrackData` semantics), float
  destinations `REALSXP`.

Each kernel takes `bool swap, bool zeroToNa` as compile-time (template) parameters so the
inner loop has no branches.

Numeric contract: identical results to today (widening casts are exact; `float`→`double` of
the same bit pattern). Verified by a SIMD-vs-scalar test in the style of
`tests/testthat/test-simd.R`.

Vectorisation policy (§3.1): write the kernels as plain C++ so the compiler's auto-vectorizer
emits the NEON/SSE/AVX sequence (measured to match or beat hand-written xsimd 7.1.3 and raw
NEON); do **not** add xsimd to this module. If a target compiler does not vectorize at its
default optimisation level, the loops still produce correct scalar results.

### 4.2 D2 — single pass, no temp buffer, no zeroing

* Delete the per-record `memcpy` into `tempBuffer` (`src/dataobj.c:451`) and the `tempBuffer`
  allocation; read descriptors straight out of the record buffer. (Measured: 14 % of read
  time, 8 ms per 82 MB per track.)
* `allocDataBuf()` uses `calloc` (`src/assp/dataobj.c:838`) → zero-then-overwrite. Use
  `malloc` on the read path (the whole buffer is filled by `fread`; `asspFFill` already calls
  `swapDataBuf` when endianness differs, which reads every byte).
* Hoist the `switch (desc->format)` out of both the record loop in `getDObjTrackData()` and
  `addTrackData()`.

### 4.3 D3 — one buffer, `mmap` where available

Read path only (`getDObj2`), not the DOBJ semantics used by DSP C code:

* POSIX: `open`+`mmap(PROT_READ, MAP_PRIVATE)` the **window** (`headerSize + begin*recordSize`,
  `numRecs*recordSize` bytes) — page-cache direct, no intermediate buffer, no copy, no zeroing.
  This is the analogue of simdjson's "read the whole document into your own buffer", minus the
  copy. Requires the kernels to handle unaligned record starts (they use `memcpy`-based loads)
  and forbids in-place swapping — which the kernels already do per value (D5).
* Windows / `mmap` failure / non-regular files: fall back to today's
  `allocDataBuf`+`asspFFill` path into a `malloc` buffer. Select with a single
  `#ifdef _WIN32` branch; the kernels take a pointer either way.
* Keep `asspFFill()`'s behaviour untouched for C-level DSP callers.

### 4.4 D4 — on-demand materialisation (skip what was not asked for)

The single biggest win available, and the direct analogue of `On Demand`:

* `read_ssff(fname, begin, end, samples, tracks = NULL)` /
  `read_track(file, begin, end, samples, validate, tracks = NULL)`: `NULL` = all (today's
  behaviour), otherwise a character vector of descriptor identifiers.
* C side: `getDObj2` gains `tracks` (character or NULL) and passes a per-descriptor "wanted"
  flag to `dobj2AsspDataObj()`/`ssff_fill_tracks()`; unwanted descriptors are skipped — no
  allocation, no conversion, no bytes touched.
* Attributes stay honest: names, matrices and `trackFormats` cover the selected tracks only
  (unselected tracks are dropped, never returned as empty matrices); `startRecord`/
  `endRecord`/`sampleRate`/`startTime` are unchanged.
* Unknown names → `cli::cli_abort()` listing the available track names (error-reporting
  standard, `R/error_helpers.R`).
* Multi-track files with a mixed shape (1 f0 field + 1025 spectrum fields) currently pay for
  everything; after D4 a `read_track(f, tracks = "f0")` reads 0.4 % of the work.

### 4.5 D5 — fused byte swap

`swapDataBuf()` (`src/assp/dataobj.c:943`) walks the whole buffer with a per-record `switch`
*before* conversion, doubling the passes over big-endian files (measured 61 ms vs 47 ms on the
same bytes). Fold the swap into the conversion kernel as a `__builtin_bswap32/64` on the loaded
value so a big-endian file costs one pass. No ISA-specific code is warranted: scalar `bswap`
measured 20.5 GB/s (memory-bandwidth bound, §2.4) and xsimd 7.1.3 has no byte-shuffle
primitive (§3.1); leave `swapDataBuf()` in place for in-memory DOBJs and the write path.

### 4.6 D6 — header parse

Rewrite `getSSFFhdr()`'s loop (`src/assp/headers.c:2607`) as:

1. `fread` up to a bounded header window in one call (grow once if the EOH marker is not
   found within it; SSFF headers are line-oriented text, so a 64 KB first read covers any real
   file and the marker scan is `memchr`),
2. split lines with `memchr('\n')`, tokenise in place into `field[]` (no `strdup` for locals),
3. keep the O(n²) generic-variable append out of the hot path (tail pointer or identifier
   index, `src/assp/dataobj.c:461`),
4. parse only the keys SSFF defines (`Machine`, `Record_Freq`, `Start_Time`, `Column`,
   `Original_Freq`, generic variables) with the same grammar and the same error messages.

Contract: every header the old parser accepted parses to an identical `DOBJ` metadata set
(same descriptors, types, units, factors, `dataRate`, `Start_Time`, `sampFreq`, generic
variables — compared with `all.equal` on the resulting `AsspDataObj` attributes), plus the
same `asspMsgNum`/`applMessage` on rejection (`tests/testthat/test-edge-cases.R:69-86` covers
the rejection path). Target: ≤ 5 µs for a 30-line header (today ≈ 38 µs including open).

### 4.7 D7 — optional block-parallel fill

For inputs above a size threshold (~8 MB of returned data), split the record range across
worker threads (each writes disjoint rows of disjoint matrices; no locks, no R API calls, no
allocation off the main thread).

* Guard with `#ifdef _OPENMP` + `#pragma omp parallel for`; `SHLIB_OPENMP_CXXFLAGS` is already
  in `src/Makevars.in`. **Do not `#include <omp.h>`** — the Homebrew libomp header fails to
  compile with Apple clang 21 in this environment (verified), and nothing in the package
  currently includes it.
* Public knob: `read_ssff(..., threads = 1L)` with `1` = serial default for the first release
  (parallelism changes timing only, never values — the worker split is by record block).
* Measured ceiling: 8 threads on the fused fill gave 3.6–4.6× over 1 thread (2.5); page
  faults on the result matrices serialise, so treat parallel as a large-file optimiser.

### 4.8 D8 — API / ABI

* `getDObj2` is registered with 4 arguments (`src/superassp_init.c:269`) and called by name
  from `R/read_ssff.R:38` and the legacy alias `R/assp_dataobj_methods.R:22`. New arguments
  (`zero_to_na`, `tracks`, `threads`) are appended as **optional tagged arguments**: the C
  side treats a missing tag as the old default, so the legacy alias keeps working unchanged.
* No new Rcpp exports are needed; the module is called from C. If any binding is added, it is
  internal (`@keywords internal` + `@noRd`, `Rcpp::compileAttributes()`) per `CLAUDE.md`.
* `src/assp/*` is the in-tree libassp fork (locally patched; it is **not** one of the six
  pinned trees in `inst/VENDORED_SOURCES.md`), so edits there are allowed — keep them
  minimal, in `getSSFFhdr`/`fgetl`/`allocDataBuf` only, and list them in `NEWS.md`.

### 4.9 D9 — `0` means missing

**Read rule (new).** In `read_track()`'s SSFF branch, for every descriptor that is **not**
sampled audio (`dd->type != DT_SMP` — the flag libassp sets for `audio`/`samples` in every
container header, `src/assp/headers.c:123-124,194,980,…`), a stored value of exactly `0`
(including `-0.0`; not NaN) is returned as `NA`:
`NA_real_` for `REAL32`/`REAL64` tracks, `NA_integer_` for integer tracks.

* The mask is fused into the conversion kernel (`v == 0.0 ? NA_REAL : v`), so it is free
  (measured: the fused fill in 2.5 already includes it and stays memory-bound).
* Audio is excluded because `0` is digital silence there — `read_track("x.wav")` and
  `read_ssff()` on a WAV must not turn silent samples into `NA`.
* NaN in the file is **not** converted (it is a computed value, not the missing convention);
  it stays NaN, which is still `is.na()` — the "no valid value" contract holds either way.
* Escape hatch: `read_ssff(file)` keeps returning stored values verbatim (default
  `zero_to_na = FALSE`), so pipelines that need the raw encoding (e.g. to distinguish
  "unvoiced" from "not analysed") are unaffected and the change is confined to the
  documented user-facing reader.
* Documented consequences, all inherent to the format: integer *indicator* tracks
  (`pm` 0/1 in `create_pitchmark_asspobj()`, `R/helpers_av_sptk.R:42`) read as `NA`/1;
  tracks where 0 is a legitimate value (a pitchmark at sample 0, an exact 0 dB frame) read as
  `NA`. `read_ssff()` is the documented way to get the encoded values.

**Write rule (new).** `addTrackData()` (`src/dataobj.c:934`) maps `NA` **and NaN** to `0`
before writing floats/integers, making the round trip stable:
`NA → 0 → NA`, `NaN → 0 → NA` (documented), `0 → NA` in `read_track` only.
Today NA_real_ is cast to a `float` NaN and lands in the file as a NaN bit pattern that no
other SSFF tool interprets and that does not follow the package's own convention.

**Everything else stays:** `as_tibble.AsspDataObj(na.zeros = TRUE)` becomes a no-op on track
data (already NA) and stays as-is for compatibility; JSTF keeps its own `null` encoding
(`R/jstf_io.R`, `na = "null"`); `read_audio()` is untouched (audio path).

## 5. Expected result

Rows marked **measured** are the fused-fill numbers from §2.5 (same files, same machine);
rows marked *estimated* are kernel floor + the measured allocation/fault overhead and carry
[INFERENCE] status.

| case | today | after D1–D7 (1 thread) | after D1–D7 (threads) | evidence |
|---|---|---|---|---|
| 82 MB 1-track spectrum | 0.047 s | 0.029 s (1.6×) | 0.008 s (5.9×) | measured |
| 82 MB 4-track spectrum | 0.037 s | 0.028 s (1.3×) | 0.009 s (4.1×) | measured |
| big-endian spectrum | 0.061 s | 0.029 s (2.1×) | 0.008 s (7.6×) | measured |
| 4 MB, 10 single-field tracks | 0.007 s | 0.001 s (7×) | 0.001 s | measured |
| 53 MB INT16 audio | 0.084 s | ~0.020 s (4×) | ~0.015 s | estimated (kernel floor 2.6 ms + allocation) |
| 400-record f0 file | 38 µs | ~15 µs (2.5×) | — | estimated (header target ≤ 5 µs + open/alloc) |
| one track out of a 4-track file (`tracks=`) | 0.037 s | ~0.010 s | — | estimated (skip 3 of 4 descriptors) |

End-to-end numbers include R's matrix allocation and GC, which the reader cannot remove;
the C-side wins in §2.5 are the reliable upper bound, the 1.2–2.3× seen with cold
allocations is the conservative floor.

## 6. Implementation order (each step lands green on its own)

| # | Change | Verification gate |
|---|---|---|
| 1 | `src/ssff_convert.cpp` + `src/dataobj.c` rewrite of `getDObjTrackData` (no `memcpy`, hoisted switch, per-format kernels, `malloc` not `calloc`, no xsimd — §3.1) | `devtools::test()`; new test comparing kernel output vs. the old implementation over all formats (float/int, 1-field/multi-field, LE/BE); benchmark before/after; `R CMD check`'s compiled-code check stays clean |
| 2 | Fused swap (`swapDataBuf` not called on the read path) | big-endian fixture read must equal the little-endian fixture value-for-value; old path removed |
| 3 | `mmap` window path + Windows fallback | windowed reads (`begin`/`end`, `samples=TRUE`) equal today's for a 100 MB file at both ends and mid-file; ASAN run |
| 4 | Header rewrite (`getSSFFhdr` bulk read, O(n) generic list) | existing header tests + `test-edge-cases.R` rejection paths; all `inst/samples` and `tests/signalfiles` files parse to identical metadata (`all.equal` on attributes) |
| 5 | `zero_to_na` in C + `read_track()` default `TRUE` + writer NA/NaN→0 | new tests (§4.9), `README`/`NEWS`; `write_ssff(read_track(f))` round-trip fixture |
| 6 | `tracks=` selection (R + C) | selection returns identical matrices to full read for the selected names; unknown name aborts with the available names |
| 7 | Optional OpenMP block fill (`threads=`) | serial == parallel results bit-for-bit; CRAN build without OpenMP still compiles |

## 7. Risks, non-goals, open decisions

* **R-side allocation is the floor.** For big matrices the remaining cost is page faults and
  GC (~50 % of the profile). D4 (skip) and windowed reads are the only levers left; a
  copy-on-write return (e.g. `ALTREP`) is explicitly out of scope.
* **Parallelism** must never run under `R CMD check`'s CRAN settings by default, must not
  allocate in workers, and must not regress small inputs — hence threshold + `threads = 1`
  default.
* **NA semantics change** is user-visible for `read_track()` only; flagged in `NEWS.md` as a
  3.1.0 change with the `read_ssff()` escape hatch documented in both man pages.
* **`mmap` on Windows** is not attempted (fallback path is the same code the package has
  today).
* **Not in scope:** the C-DSP input path (`performAssp.c:1085` `asspFOpen` per call, which
  re-reads and re-parses the file for every analysis function), JSTF/`RcppSimdJson` tuning,
  and any writer-side performance work beyond NA/NaN→0.
* **xsimd is not adopted in the reader** (§3.1, measured): the available RcppXsimd pins xsimd
  7.1.3, which cannot express the widening kernels the reader needs, and where it can express
  them it is 43 % slower than the compiler's own vectorization. xsimd remains in use for the
  compute-bound arithmetic kernels in `src/simd_utils.hpp`.
* **Open decision 1:** whether `read_ssff()` should also default to `zero_to_na = TRUE` in a
  later release (kept `FALSE` here to preserve internal callers and the escape hatch).
* **Open decision 2:** whether the writer should emit a once-per-file `cli` message when it
  rewrites NA/NaN as 0 (silent today, consistent with `as_tibble(na.zeros = TRUE)`).

## Appendix A — reproduction

```r
# corpus (shapes as in table 2.1), then baseline
library(superassp)
# 1-track spectrum: matrix(rnorm(20000*1025), 20000, 1025), trackFormats "REAL32"
# written with write_ssff() after attr(obj, "fileInfo") <- c(20L, 2L)
system.time(read_ssff("spectrum.ssff"))                      # 2.1
system.time(for (i in 1:2000) read_ssff("tiny.f0")) / 2000   # 2.3 (µs/call)
```

```bash
# profile (2.2)
Rscript -e 'for (;;) invisible(superassp::read_ssff("multispec.ssff"))' &
sample <pid> 12 -file prof.txt
```

The fused-fill kernel of §2.5 is ~90 lines of C++ (`mmap` + `nrec × (track, field)` loop with
`NA_REAL` masking) built with `R CMD SHLIB`, driven from R with `system.time()` medians over
7 reps and validated with `all.equal()` against `read_ssff()`.

xsimd comparison (§3.1), ad-hoc, outside the package build:

```bash
XI=$(Rscript -e 'cat(system.file("include", package = "RcppXsimd"))')
clang++ -arch arm64 -isystem "$(xcrun --show-sdk-path)/usr/include/c++/v1" -I"$XI" \
        -O3 -std=c++17 -fno-fast-math -o xsimd_kernels xsimd_kernels.cpp   # scalar vs xsimd loop
clang++ ... -fsyntax-only xsimd_cast_probe.cpp    # xsimd::batch_cast<double>(batch<float>) -> error on 7.1.3
clang++ ... -S vec_probe.cpp                      # compare emitted NEON (fcmeq.2d / bit.16b)
```

R-side reference for the spectra shape (§3.1): `system.time(t(matrix(rnorm(20000*1025), 20000, 1025)))`.

## Appendix B — invariants the rewrite must not break

1. `read_ssff()` returns exactly the same `AsspDataObj` today's reader returns for every file
   in `inst/samples`, `tests/signalfiles` and `tests/testthat/golden` (values, dims, names,
   `sampleRate`, `startTime`, `startRecord`, `endRecord`, `trackFormats`, `fileInfo`,
   `filePath`).
2. Windowed reads are record-exact and unchanged (`getDObj2` timing arithmetic,
   `src/dataobj.c:99-127`).
3. `read.AsspDataObj`/`getAsspDataObj` (legacy alias) keep their 4-argument call.
4. Integer tracks stay `INTSXP`; float tracks stay `REALSXP`; no change of in-memory type.
5. C-level DSP callers see unchanged `asspFOpen`/`asspFFill`/`swapDataBuf` semantics.
6. No `-march=native`, no new hard dependency; `RcppSimdJson` stays in `Suggests` with the
   jsonlite fallback.

## 8. Implementation record — 2026-09-19

### 8.1 What landed

| Plan item | Implementation |
|---|---|
| D1 kernels | `src/ssff_convert.{hpp,cpp}`: one templated loop per (storage format × destination type) with the swap/NA flags as compile-time parameters; the format switch sits outside the loops. Multi-field descriptors use the **record-major** nest (contiguous loads inside a record, one 64-byte-line-fitting run per destination column), single-field descriptors the contiguous/strided column loop. Both were measured; field-major lost by 8–15 % on the 1025-field shapes, so `numFields > 1` now takes the record-major path. |
| D2 single pass | `dobj2AsspDataObjEx()` allocates every requested matrix up front and fills them in one call to `ssff_convert_records()`; the per-record `memcpy` into `tempBuffer` and `getDObjTrackData()` are gone. |
| D3 mmap | `ssffMapWindow()` maps the requested records when the window is ≥ 256 kB (below that the syscall overhead is not worth it) and calls `madvise(MADV_WILLNEED)` to keep faults clustered; the previous `allocDataBuf`+`asspFFill` path is the fallback (Windows, non-regular files, small windows) and already yields host-ordered data, so conversion then runs with swapping off. |
| D4 selection | `tracks=` on `read_ssff()`/`read_track()`; unselected descriptors are skipped before allocation. Unknown names abort from C with the file's available track list (a C-level `error()` rather than the planned `cli_abort()` — the R side never knows the file's track names without reading the file; the message carries the same information). |
| D5 fused swap | Conversion loads are endian-aware (`__builtin_bswap32/64` on the loaded word); `swapDataBuf()` is no longer called on the file-read path. |
| D6 header | `SSFFLineReader` in `src/assp/headers.c`: one buffered read of up to 64 kB (clamped to the file size), line splitting with the same fgetl semantics (LF/CR/CRLF, over-long-line truncation) and an automatic fgetl fallback when the complete header is not in the buffer; generic variables append through a tail pointer (`addTSSFF_GenericAfter()`), removing the O(n²) walk. |
| D7 threads | `threads =` argument; `#pragma omp parallel for schedule(static) num_threads(n) if (n > 1)` over record blocks, guarded by `#ifdef _OPENMP`, no `<omp.h>` include. Serial by default. |
| D8 API | `getDObj2` is registered with a variable arity (`-1`) so the legacy `read.AsspDataObj` four-argument call keeps working; the new tags (`zero_to_na`, `tracks`, `threads`) are optional. |
| D9 semantics | `zero_to_na` is fused into the kernels per descriptor (`desc->type != DT_SMP`); `read_track()` defaults it to `TRUE`, `read_ssff()` to `FALSE`. `addTrackData()` maps `NA` and `NaN` to `0` on write. |

Tests: `tests/testthat/test-ssff-read-engine.R` (63 assertions) checks the reader against
an independent `readBin()` decode (all storage types, multi-track/multi-field layout,
windowed reads down to individual records, big-endian fixtures, thread-count invariance,
track selection, the zero/NA rules and an NA/NaN round trip), plus a 300-generic-variable
header. The full suite passes on the final build: **2958 passed, 0 failed, 0 errors,
40 skipped**, and `test-ssff-read-engine.R` passes 63/63 both at `-O0` and at `-O2`.

### 8.2 Measured before/after (same machine, fresh process per measurement)

| case | old reader | new reader | change |
|---|---|---|---|
| 82 MB 1-track spectrum REAL32 | 36–39 ms | 36–40 ms | parity |
| 82 MB 4-track spectrum REAL32 | 31–34 ms | 28 ms | 1.1–1.2× |
| 82 MB big-endian spectrum | 56–60 ms | 56–58 ms | parity |
| 41 MB 1-track spectrum REAL64 | 10–12 ms | 6–9 ms | 1.4–2× |
| 4 MB, 10 single-field tracks | 7 ms | 1 ms | 7× |
| 53 MB INT16 audio (WAV) | 82 ms | 11–13 ms | 6.5× |
| 10.6 MB INT16 audio (stereo) | 15 ms | 3–4 ms | 4× |
| 400-record f0 track | 40 µs | 40–50 µs | parity |
| 2000-line header | 3.5 ms | 1.2 ms | 3× |

The audio case is the one that matters most in daily use: every `trk_*` wrapper loads the
signal through `read_ssff()`/`av_to_asspDataObj()` before doing any DSP, so a 6.5× faster
audio read is a per-call win, not a micro-benchmark.

### 8.3 Deviations from the plan and what they cost

* **The f32 spectrum shapes are at parity, not 3–8×.** The plan's §2.5 prototype measured a
  single-thread fill of 29 ms vs 47 ms for `read_ssff`, but that excluded the 164 MB
  destination allocation and R's GC, which are identical for both readers and dominate this
  shape (a freshly allocated 20 000 × 1025 double matrix alone costs ~20–30 ms of page
  faults). The conversion itself did get 1.2–2× faster (record `memcpy` removed, dispatch
  hoisted); the end-to-end number is therefore flat for this extreme shape. Every other
  shape is faster, up to 7×.
* **`read_ssff(..., tracks=)` unknown names** report from C (see 8.1).
* **A diagnostic worth keeping in mind:** `pkgbuild::compile_dll()` (and therefore
  `devtools::load_all()`) compiles with `-O0` debug flags by default. All benchmarks and
  disassembly checks above were made with `debug = FALSE`; at `-O0` the templated kernels are
  not inlined and the new reader looks ~2× slower than the old one, which is an artefact of
  the development build, not of the code. CRAN and `R CMD INSTALL` use `-O2`.

### 8.4 Not done (still open)

* The C-DSP input path (`performAssp.c` opening the file per analysis) is untouched, as
  planned.
* `read_ssff()` still defaults to `zero_to_na = FALSE` (open decision 1 in §7); the writer
  does not emit a message when it maps `NA`/`NaN` to `0` (open decision 2).
* `threads =` is serial by default and was verified for identical results, not benchmarked
  end to end (the record split only pays off for files far larger than the test corpus;
  the prototype measured 3.6–4.6× on the block-parallel fill).
