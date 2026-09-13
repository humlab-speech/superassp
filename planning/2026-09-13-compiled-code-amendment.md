# Compiled-code amendment plan — clearing the CRAN `checking compiled code` NOTE

Supersedes Phase 6 (Tasks 6.1–6.5) of `planning/2026-09-12-winbuilder-amendment.md`.
That phase assumed the only residue would be `rand`/`srand` in `snack_formant.cc`. This
plan removes every banned symbol family, including `rand`/`srand`, without changing DSP
output on the maintainer's platform.

## Artifacts and provenance

| Artifact | Provenance |
|---|---|
| maintainer-reported NOTE | `R CMD check` on a UNIX build, `superassp 3.0.0`, symbols: `std::cerr`, `stderr`, `stdout`, `printf`, `putchar`, `puts`, `rand`, `srand` |
| `00check.log` (6943 lines) | win-builder `R CMD check --as-cran`, `superassp_3.0.0.tar.gz`, `gcc/g++ 14.3.0`, same symbol families plus `std::cout`, `exit` |
| `src/superassp.so` (local) | macOS build present in the tree; re-analysed here with `tools:::check_so_symbols()` |

Instrument used for every claim below: `tools:::.check_compiled_code` scans
`<pkg>/libs/*.{so,dylib,dll}` only (never `inst/`), maps undefined symbols through
`tools:::so_symbol_names_table`, and reports each hit with its object file when
`libs/symbols.rds` exists.

Banned set relevant here (from `tools:::so_symbol_names_table` on R 4.6.1):

- C: `abort`, `assert`, `exit`, `_exit`, `_Exit`, `printf`, `puts`, `putchar`, `stderr`,
  `stdout`, `sprintf`, `vprintf`, `vsprintf`, `rand`, `random`, `rand_r`, `srand`,
  `srandom`, `srand48`
- C++: `std::cout`, `std::cerr`

**Not banned** — and therefore usable: `fprintf`, `fputs`, `fputc`, `vfprintf`,
`snprintf`, `vsnprintf`, `fflush`, `fwrite`. This is what makes the fix below mostly a
one-token-per-site change instead of a rewrite of every call.

Verified empirically in this session: `R CMD SHLIB` of a TU calling `abort()` produces
`Found '_abort', possibly from 'abort' (C)` from `tools:::check_so_symbols()`; a
`fprintf(f, ...)`-only TU produces nothing.

## Baseline: authoritative per-object inventory

`nm -u` over every object that ends up in the package shared object (macOS build;
`/tmp` scratch build maps the same sources). `_printf`/`_fprintf`/`_snprintf` shown as
raw nm names.

### SPTK (`src/SPTK`, fork `humlab-speech/SPTK` @ `1ddf910`)

| Object | Banned |
|---|---|
| `SPTK/src/utils/sptk_utils.o` | `std::cerr` |
| `SPTK/third_party/REAPER/core/track.o` | `stderr` |
| `SPTK/third_party/REAPER/epoch_tracker/epoch_tracker.o` | `stderr`, `stdout` |
| `SPTK/third_party/REAPER/epoch_tracker/fd_filter.o` | `stderr` |
| `SPTK/third_party/REAPER/wave/codec_riff.o` | `stderr` |
| `SPTK/third_party/REAPER/wave/wave.o` | `stderr` |
| `SPTK/third_party/REAPER/wave/wave_io.o` | `stderr` |
| `SPTK/third_party/SWIPE/swipe.o` | `stderr` |
| `SPTK/third_party/SWIPE/vector.o` | `printf` |
| `SPTK/third_party/Snack/jkGetF0.o` | `stderr`, `printf` |
| `SPTK/third_party/Snack/sigproc.o` | `stderr` |

Header sites reachable from those TUs (same fix applies): `REAPER/core/track.h`,
`REAPER/core/float_matrix-inl.h`, `REAPER/epoch_tracker/epoch_tracker.h`,
`REAPER/wave/codec_api-inl.h`, `REAPER/wave/wave.h`, `REAPER/wave/wave_io-inl.h`,
`Snack/jkGetF0.h`.

Not compiled today and therefore out of scope: `SPTK/src/main/*`, `SPTK/tools/*`,
`src/ESTK/*`, `src/tcl-snack/*`, `src/pyin/*`, `src/Yin-Pitch-Tracking/*`.

### tandem (`src/tandem`, fork `humlab-speech/tandem` @ `817652e`)

| Object | Banned |
|---|---|
| `tandem/tandem_64/feature.o` | `printf` (on GCC also `putchar` — `printf("\n")` folds to `putchar`) |
| `tandem/tandem_64/mScaleInten.o` | `printf` (idem) |

Only these two of the eight `TANDEM_SOURCES` reference the printf family.

### ours

| Object | Banned |
|---|---|
| `snack_formant.o` | `rand` (`src/snack_formant.cc:406`, used at `:430` dither; `:612-613` root-finder restart) |

### openSMILE (in-tree `src/opensmile`, built by `configure` → `libopensmile.a`)

| Object | Banned |
|---|---|
| `src/dsp/vadV1.cpp.o` | `printf` |
| `src/rnn/rnnVad2.cpp.o` | `printf` |
| `src/lldcore/mfcc.cpp.o` | `printf` |
| `src/functionals/functionalSegments.cpp.o` | `printf`, `puts` |
| `src/functionals/functionalPeaks2.cpp.o` | `printf`, `puts` |
| `src/examples/exampleSink.cpp.o` | `printf` |
| `src/iocore/dataPrintSink.cpp.o` | `printf`, `stdout` |
| `src/lld/pitchDirection.cpp.o` | `puts` |
| `src/smileutil/smileUtilSpline.c.o` | `puts` |
| `src/core/componentManager.cpp.o` | `puts` |
| `src/smileutil/smileUtil.c.o` | `stderr` |
| `src/core/smileLogger.cpp.o` | `stderr` (+ `isatty(fileno(stderr))`, `fflush(stderr)`) |
| `src/other/vectorOperation.cpp.o` | `rand` |
| `src/other/maxIndex.cpp.o` | `rand` |
| `src/dsp/signalGenerator.cpp.o` | `rand`, `srand` |
| `src/newmat/newmatnl.cpp.o` | `std::cout` |
| `src/newmat/myexcept.cpp.o` | `std::cout`, `exit` |

`progsrc/smilextract/SMILExtract.cpp.o` also has `puts`, but the executable is built into
`inst/opensmile/bin/` (`.Rbuildignore`d from the tarball) and is **not** scanned by
`checking compiled code`. It must keep working, which constrains how openSMILE is fixed
(see F2).

### locally built `src/superassp.so`, as the checker sees it

```
Found '___stderrp'  (stderr)   Found '___stdoutp' (stdout)   Found '__ZNSt3__14cerrE' (std::cerr)
Found '_printf'     (printf)   Found '_puts'      (printf/puts)  Found '_rand' (rand)  Found '_srand' (srand)
Found '_abort'      (abort)
```

`_abort` does **not** come from the package: it is `Rcpp/r_cast.h:74`, compiled under the
library's `#ifndef NDEBUG` branch (upstream Rcpp guard is inverted: without `-DNDEBUG`
the "impossible cast" path calls `abort()`, with `-DNDEBUG` it throws
`Rcpp::not_compatible`). CRAN and win-builder builds pass `-DNDEBUG` (see
`00install.out.txt:1389`), a plain local `R CMD INSTALL` here does not
(`R CMD config CXX17FLAGS` → `-falign-functions=64 -Wall -g -O2`). Hence it is absent from
the maintainer's pasted NOTE and from the win-builder log, and present in the tree's
`src/superassp.so`. See Task 1.3.

## Findings

### F1 — the offending sites are message sinks, not data paths

Every `stderr`/`stdout`/`printf`/`puts`/`putchar`/`cout`/`cerr` site in the compiled set is
a diagnostic or progress message. No site writes to a *data* file through those entry
points — REAPER/tandem/openSMILE data output goes through `fprintf(fr->fp(), …)`,
`fwrite`, or C++ streams on a caller-supplied file. `fprintf` is not banned, so those
calls are untouched. This is what makes a mechanical, one-token substitution safe.

Site counts (raw token matches, comments/`#if 0` regions not excluded — the executing
script filters both, and the `nm` inventory above is the contract):

| File | sites |
|---|---|
| `SPTK/third_party/REAPER/core/track.cc` (+ `track.h`, `float_matrix-inl.h`) | 20 + 4 + 2 |
| `SPTK/third_party/REAPER/epoch_tracker/epoch_tracker.cc` | 18 |
| `SPTK/third_party/REAPER/wave/codec_riff.cc` | 16 |
| `SPTK/third_party/REAPER/epoch_tracker/fd_filter.cc` | 9 |
| `SPTK/third_party/REAPER/wave/wave_io.cc` / `wave.cc` (+ `-inl.h`, `wave.h`, `codec_api-inl.h`) | 5 / 3 (+4) |
| `SPTK/third_party/Snack/jkGetF0.cc` | 22 (all through `Fprintf`, `jkGetF0.h:86`; one direct `printf` at `:593`) |
| `SPTK/third_party/Snack/sigproc.cc` | 6 (all `Fprintf`) |
| `SPTK/third_party/SWIPE/swipe.cc` | 9 enabled (`fprintf(stderr, …)`); the CLI demo (`main`, 30 prints, 11 `exit`) is inside `#if 0` at `:489-713` and is not compiled |
| `SPTK/third_party/SWIPE/vector.cc` | 6 enabled (`printv:106`, `printiv:225`, `printm:360-361`, `printim:442-443`); demo `main` is `#ifdef DEBUG`, not compiled |
| `SPTK/src/utils/sptk_utils.cc` | 1 (`std::cerr << stream.str()` at `:510`) |
| `tandem/tandem_64/feature.cpp` / `mScaleInten.cpp` | 10 / 5, all `if (echo)`-gated progress dots |
| `snack_formant.cc` | 3 `rand()` |
| openSMILE | 17 objects, ≤ a handful of sites each; `smileLogger.cpp` is the concentrated one |

Multi-line `fprintf(stderr,` sites that a line-oriented rewrite cannot handle:
`epoch_tracker.cc:1163`, `fd_filter.cc:165,182,185,229`, plus the corresponding
`Fprintf(stderr,` continuations in Snack. Six hand edits.

### F2 — `libopensmile.a` must not reference the R API

`configure` builds `SMILExtract` from the same CMake tree and copies it to
`inst/opensmile/bin/SMILExtract`; `R/list_cpp_opensmile_emobase_helper.R:16` runs it
through `system()`. The executable links `libopensmile.a` and nothing else R-side, in a
process where no R runtime is initialised. Any `Rprintf`/`REprintf`/`R_Consolefile`
reference inserted into library objects therefore either breaks the link or calls into an
uninitialised R at runtime. `R_Consolefile`/`R_Outputfile` are additionally declared only
in `Rinterface.h`, which is not a public package header.

Consequence: openSMILE's console output must go through a host-installed sink, not
through R directly. The R side supplies an `Rvprintf`/`REvprintf` writer from *our own*
TU (which is compiled with R headers), the executable supplies a `vfprintf` writer from
*its* TU.

### F3 — `rand()` in `snack_formant.cc` can be replaced with bit-identical output on macOS

`frand_val()` feeds a dither term, `sig[i] = data[i] + 0.016*frand_val() - 0.008`
(`snack_formant.cc:430`) in the stabilised-covariance LPC, and `:612-613` perturb `p`/`q`
when the Newton iteration fails to converge. There is no `srand()` call, so the sequence
is libc's default-seeded one.

macOS's `rand()` is MINSTD (Park–Miller, `a=16807`, `m=2^31-1`) — verified in this
session: `16807, 282475249, 1622650073, 984943658, 1144108930, 470211272` is exactly
macOS `rand()`'s output and exactly the MINSTD recurrence. glibc's `rand()` is the
TYPE_3 additive-feedback generator and produces a *different* sequence, as does MSVC's.

So a package-local MINSTD generator is **bit-identical to today's macOS output**
(the platform the regression suite and the maintainer's checks run on), while making the
formant pipeline platform-invariant. Nothing else changes: the amplitude stays
`/RAND_MAX`-scaled and the dither term is ~1e-7 of full scale on an int16-scaled signal
(`0.016` counts out of `32768`).

### F4 — openSMILE's `rand`/`srand` sites are unreachable from the bundled configs

`grep -rl "signalGenerator\|vectorOperation\|maxIndex\|dataPrintSink" inst/opensmile/config`
over the 12 shipped `.conf` files returns nothing. The three sites are: `signalGenerator`
(a synthetic-signal *source*, seeded from `time(0)` at `signalGenerator.cpp:256` and
already non-reproducible), `vectorOperation` (optional polar-method noise), and `maxIndex`
(`randNoise` defaults to `0.0`, so its `rand()` call is multiplied by zero). Replacing the
generator cannot change any documented/validated analysis path; only a user-supplied config
that instantiates those components sees different random draws.

### F5 — `exit` is a single site, in dead code for this build

`exit` appears only in `src/newmat/myexcept.cpp:231` (`Terminate()`), which
`myexcept.h:208-209` reaches from `Throw()`. `newmat` macro-dispatches `THROW` to `throw`
when exceptions are enabled and to `Throw()` otherwise; no `-fno-exceptions` appears in
any build log. Removing the call is therefore a symbol fix with no runtime path in this
package, and `throw std::runtime_error` is the correct replacement in either case.

## Constraints

- **No DSP numeric change on macOS.** F3 is the only place where a banned API sits on a
  numeric path, and the MINSTD replacement is sequence-identical there.
- **`SMILExtract` keeps working**, and stays free of R references (F2). It is not scanned
  by the check, so its own `printf`-family use needs no change.
- **No warning-suppressing flags.** The prior amendment removed `-Wno-register`
  specifically to close W3/W6; do not reintroduce `-Wno-…` to hide format diagnostics
  (it would also be a policy violation in its own right).
- Vendored-fork edits land in `humlab-speech/SPTK` and `humlab-speech/tandem` and are
  pointer-bumped here; `src/opensmile` is in-tree and edited directly.
- `Rcpp/r_cast.h` is an installed dependency: nothing in this package can change it.
  The only lever is `-DNDEBUG` (Task 1.3).

## Phases

### Phase 1 — Acceptance tooling and baseline (blocks every later phase)

- [ ] **Task 1.1 — symbol gate.** Add `tools/check_cran_symbols.R`: takes an installed
  package dir (default: `find.package("superassp")`), runs
  `tools:::check_so_symbols()` on each `libs/**/*.<dynlib.ext>`, prints the symbol/object
  table, and `quit(status = 1)` when anything is found. This is the pass/fail contract for
  the whole plan and is reusable on win-builder logs (`nm` alone is not enough — the
  checker's own table is the ground truth).
- [ ] **Task 1.2 — record the baseline** by running Task 1.1 against an
  `R CMD INSTALL`-built library (never a `devtools::load_all()` build — see 1.3), and
  paste the output into this file's implementation record.
- [ ] **Task 1.3 — decide `NDEBUG`.** Add `-DNDEBUG` to `PKG_CPPFLAGS` in
  `src/Makevars.in` (it is already present in every CRAN/win-builder compile command, see
  `00install.out.txt:1389`) so local builds stop emitting `Rcpp/r_cast.h`'s `abort()`. If
  the flag is deliberately *not* added, `_abort` must be listed as accepted residue with
  the `r_cast.h` provenance — but then the gate in Task 1.1 fails on machines whose R
  default flags omit `-DNDEBUG`, which is a bad gate to hand a co-maintainer.

**Verification:** Task 1.1 fails on the current build with exactly the symbol families in
the Baseline table, and the object attribution matches that table.

### Phase 2 — SPTK and tandem forks (submodule edits + pointer bumps)

One scripted pass, executed in each fork checkout, then reviewed as a diff. Substitution
rules, applied only to files listed as offenders in the Baseline inventory (plus the six
hand-edited multi-line sites):

| From | To |
|---|---|
| `fprintf(stderr,` | `REprintf(` |
| `Fprintf(stderr,` | `REprintf(` |
| `fprintf(stdout,` | `Rprintf(` |
| `fflush(stderr);` / `fflush(stdout);` | *(deleted)* |
| `printf(` | `Rprintf(` |
| `puts(x)` | `Rprintf("%s\n", x)` |
| `putchar(c)` | `Rprintf("%c", c)` |
| `std::cerr << stream.str();` | `REprintf("%s", stream.str().c_str());` |

Plus `#include <R_ext/Print.h>` in each touched file and a trailing `\n` where a message
lacked one.

- [ ] **Task 2.1 — REAPER** (`core/`, `epoch_tracker/`, `wave/`, including the `-inl.h`
  and plain headers), message-only sites per F1.
- [ ] **Task 2.2 — Snack** (`jkGetF0.cc`, `sigproc.cc`). After the rewrite the
  `#define Fprintf (void)fprintf` at `jkGetF0.h:86` has no users; delete it so it cannot
  be reintroduced. `jkGetF0.cc:593`'s direct `printf` → `Rprintf`.
- [ ] **Task 2.3 — SWIPE** (`swipe.cc` enabled sites, `vector.cc` print helpers). Leave
  the `#if 0` CLI demo and the `#ifdef DEBUG` demo alone: they are not compiled and
  rewriting them would put dead code in the diff.
- [ ] **Task 2.4 — SPTK utils** (`sptk_utils.cc:510`) → `REprintf("%s", …)`.
- [ ] **Task 2.5 — tandem** (`feature.cpp`, `mScaleInten.cpp`) → `Rprintf`. Keep the
  `echo` gate exactly as is; only the sink changes.
- [ ] **Task 2.6 — bump both submodule pointers** and record the new SHAs here.

**Verification:** rebuild; `nm -u` on the eleven SPTK and two tandem objects (or the
`symbols.rds` attribution from Task 1.1) shows none of the banned names; the package still
builds on macOS and Windows; `devtools::test()` count unchanged.

### Phase 3 — `snack_formant.cc` PRNG (closes the `rand` family)

- [ ] **Task 3.1 — add the generator.** Anonymous-namespace MINSTD in
  `src/snack_formant.cc`: `unsigned int` state seeded to 1, `state = state*16807 % 2147483647`,
  returning `int`, with a matching `rand_max()` constant of `2147483647` used in place of
  `RAND_MAX` at the three sites (`:406`, `:612`, `:613`). Document at the definition why
  libc `rand` is unusable (CRAN policy) and why MINSTD was chosen (byte-for-byte the macOS
  sequence, per F3).
- [ ] **Task 3.2 — verify numerics.** Run the formant-related tests, and diff
  `trk_formant`/Snack-formant output on a fixed corpus (e.g. the bundled example audio)
  against the pre-change binary on macOS. Expected: identical values. If a restart path
  (`:612-613`) fires on the corpus, report the observed delta rather than assuming zero.
- [ ] **Task 3.3 — report the platform note.** State in `NEWS.md`/`cran-comments.md` that
  on glibc/Windows the dither sequence changes from libc's to the fixed MINSTD one, i.e.
  formant output goes from platform-dependent to platform-invariant above the 1e-7 level.

**Verification:** Task 1.1 no longer names `snack_formant.o`; the corpus diff from 3.2 is
empty on macOS; `devtools::test()` unchanged.

### Phase 4 — openSMILE console sink (the largest surface, and the one with a build constraint)

- [ ] **Task 4.1 — add the sink.** New
  `src/opensmile/src/include/smileutil/smileConsole.h` + `src/opensmile/src/smileutil/smileConsole.cpp`,
  compiled into `libopensmile.a`:
  ```c
  typedef void (*smile_console_writer)(const char *fmt, va_list ap);
  void smile_console_set_writer(smile_console_writer out, smile_console_writer err);
  int  smile_console_printf(const char *fmt, ...);   /* printf replacement, same signature/return */
  int  smile_console_puts(const char *s);
  int  smile_console_putchar(int c);
  void smile_console_error(const char *fmt, ...);    /* stderr replacement */
  int  smile_console_is_tty(void);
  #ifdef __cplusplus
  std::ostream &smile_console_out();                 /* cout replacement, formatting preserved */
  std::ostream &smile_console_err();
  #endif
  ```
  Implementation writes to the installed writers; with none installed it is silent. It
  uses `vfprintf`/`fwrite` on writers supplied by the host — no banned symbol inside the
  library, and no R dependency (F2).
- [ ] **Task 4.2 — install the writers from the executable.**
  `progsrc/smilextract/SMILExtract.cpp` installs `vfprintf(stdout, …)` /
  `vfprintf(stderr, …)` writers at the top of `main`, preserving today's CLI behaviour.
- [ ] **Task 4.3 — install the writers from R.** In `src/opensmile_wrapper.cpp`
  (compiled with R headers), install writers that call `Rvprintf` / `REvprintf`
  (`R_ext/Print.h` — public API, no `R_Consolefile`). Install once, before the first
  `smile_new()` (`opensmile_wrapper.cpp:127`), via the existing init/entry path.
- [ ] **Task 4.4 — convert the C/C++ sites in the 17 library objects** using
  identifier-level substitutions that preserve signatures and formatting:
  `printf(` → `smile_console_printf(`, `puts(` → `smile_console_puts(`,
  `putchar(` → `smile_console_putchar(`, `cout` → `smile_console_out()`,
  `cerr` → `smile_console_err()`, `fprintf(stderr, …)` → `smile_console_error(…)`,
  `fprintf(stdout, …)` → `smile_console_printf(…)`.
  Site-specific work: `smileLogger.cpp:71-72` (`isatty(fileno(stderr))` →
  `smile_console_is_tty()`, `:245-246` → `smile_console_error` and drop the `fflush`),
  `dataPrintSink.cpp` (`stdout` stream), `newmat/newmatnl.cpp` (4 `cout` sites),
  `newmat/myexcept.cpp` (`cout` sites plus `exit(1)` at `:231` → `smile_console_error(…)`
  followed by `throw std::runtime_error(…)`, per F5).
- [ ] **Task 4.5 — RNG sites.** Same MINSTD as Task 3.1, exposed as
  `smileRandomUniform()` / `smileRandomSeed()`; `signalGenerator.cpp:205,243,244,272,296`,
  `vectorOperation.cpp:561`, `maxIndex.cpp:149`. Justify with F4 (no bundled config
  reaches these components; `maxIndex` multiplies by a zero default).
- [ ] **Task 4.6 — toolchain sweep.** Run the same site pattern over *all* opensmile
  sources that go into `libopensmile.a` (not just the 17 objects) so that a different
  optimiser's `printf`→`puts`/`putchar` folding (the reason the win-builder log shows
  symbols the macOS build does not) cannot surface new entries. Exclude
  `src/video/`, `include/android|ios`, `progsrc/include/tests/`, and the `smilextract`
  target.
- [ ] **Task 4.7 — build both consumers.** Confirm `configure` still produces both
  `libopensmile.a` and `SMILExtract`, that the executable links without R, and that
  `R/list_cpp_opensmile_emobase_helper.R` still returns its expected features.

**Verification:** Task 1.1 no longer names any `libopensmile.a` object; `SMILExtract -h`
and a full config run behave as before; `lst_cpp_opensmile_emobase*` tests pass; the
openSMILE diagnostics that used to appear on stderr now appear on R's console via
`REvprintf`.

### Phase 5 — Re-check on both toolchains

- [ ] **Task 5.1 — local `--as-cran` check** with the `pladdrr` stub recipe from the
  2026-09-12 plan; the `checking compiled code` section must be absent.
- [ ] **Task 5.2 — win-builder** (R-devel and R-release). The Windows/Linux symbol set is
  a superset of macOS's (`puts`, `putchar`, `std::cout` as `_ZSt4cout`), so this is the
  real test of Task 4.6.
- [ ] **Task 5.3 — rewrite the compiled-code paragraph of `cran-comments.md`** against
  the final logs: symbol families fixed, mechanisms used, and any residue with its reason
  (none expected; if Task 1.3 is declined, `abort` with its `Rcpp/r_cast.h` provenance).

**Verification:** `R CMD check --as-cran` shows no `checking compiled code` NOTE/WARNING
on macOS and on win-builder; `Status:` otherwise unchanged from the 3.0.0 baseline.

## Risks

| Risk | Mitigation |
|---|---|
| `Rprintf`'s 8 kB buffer truncating a long vendored message | Every site is a short diagnostic; Task 4.1's sink is a thin wrapper, and the R writer can chunk if ever needed |
| `REprintf` without a trailing newline buffering oddly in Rgui | Append `\n` at sites whose messages lacked one (F1 tail) |
| Fork patches lost on submodule update | Patches live upstream in `humlab-speech/SPTK` and `humlab-speech/tandem`; record the SHAs in Task 2.6 |
| openSMILE build breakage on Windows (`configure.win` uses the same CMake tree) | Task 4.2 keeps writers in the executable TU; Task 4.7 builds both consumers before the win-builder run |
| Numeric drift from the PRNG change | MINSTD is sequence-identical to macOS libc `rand()` (F3), verified by Task 3.2's corpus diff |
| A future edit reintroduces a banned symbol | Task 1.1's gate is a repo script, runnable in CI alongside the existing checks |

## What this plan does not do

- It does not touch any data path, any `fprintf` on a real `FILE *`, or any DSP
  mathematics. `snack_formant.cc`'s three `rand` sites are the only numeric-adjacent
  changes and they are sequence-preserving on macOS.
- It does not silence compiler warnings or add `-Wno-*`; the vendored warnings are
  separate items in the 2026-09-12 plan (Phase 5).
- It does not patch `Rcpp`, `RcppArmadillo`, or any installed dependency; `-DNDEBUG`
  (Task 1.3) is the only dependency-side interaction.
- It does not change `SMILExtract`'s user-visible behaviour, and does not remove it from
  `inst/` (a separate question from this NOTE).


# Implementation record — 2026-09-13

## Result

| | before | after |
|---|---|---|
| macOS `src/superassp.so`, `tools/check_cran_symbols.R` | `stderr`, `stdout`, `std::cerr`, `printf`, `puts`, `rand`, `srand`, `abort` → **exit 1** | no entries → **exit 0** |
| win-builder object set (`00check.log:1460-1476`) | `std::cerr`, `std::cout`, `exit`, `putchar`, `puts`, `rand`, `srand`, `stderr`, `stdout`, `printf` | none left in the compiled sources (symbol and source sweeps below) |
| `testthat` | — | `FAILED 0 | WARN 8 | SKIP 40 | PASS 2724` (763 blocks) |
| Snack formant output, libc `rand()` vs MINSTD | — | 0 diff lines over 9 tracker configurations (1791 lines each) |

Verification commands:

```sh
R CMD INSTALL --library=/tmp/superassplib --no-docs --no-help --install-tests .
Rscript tools/check_cran_symbols.R /tmp/superassplib/superassp
nm -u <each .o> | awk '{print $NF}' | grep -xE '_printf|_puts|_putchar|___stderrp|___stdoutp|_rand|_srand|_exit|__ZNSt3__1(4cout|4cerr)E'
```

## Phase 1 — gate and baseline

- `tools/check_cran_symbols.R` added: runs `tools:::check_compiled_code()` over an
  installed package and exits non-zero on any hit. Baseline against the pre-change
  `src/superassp.so` reproduced the reported list exactly, plus `_abort`.
- `-DNDEBUG` added to `PKG_CPPFLAGS` in `src/Makevars.in`. Confirmed safe against
  `tools:::.check_make_vars`' bad-flag regexp (`-O*|-W*|-w|-ansi|-pedantic|-traditional|-f*|-m*|-std*|-isystem|-x|-pipe|-cpp|-g|-q`).
  The `_abort` reference is `Rcpp/r_cast.h:74` under `#ifndef NDEBUG`; CRAN/win-builder
  compile lines already carry `-DNDEBUG`, a local `R CMD INSTALL` here did not.

## Phase 2 — SPTK and tandem

19 files rewritten, all message-only sites (verification: zero banned tokens remain in the
compiled SPTK/tandem/package sources, and `assp/`'s remaining `printf`/`stdout` sites are
inside `#ifndef WRASSP`, which `-DWRASSP` excludes):

- REAPER: `core/track.{cc,h}`, `core/float_matrix-inl.h`, `epoch_tracker/epoch_tracker.cc`,
  `epoch_tracker/fd_filter.cc`, `wave/codec_riff.cc`, `wave/codec_api-inl.h`,
  `wave/wave.{cc,h}`, `wave/wave_io.cc`, `wave/wave_io-inl.h`
- Snack: `jkGetF0.cc`, `sigproc.cc`, and `#define Fprintf (void)fprintf` deleted from
  `jkGetF0.h:86` after its last user was rewritten
- SWIPE: `swipe.cc` (9 enabled sites; the CLI demo is inside `#if 0`), `vector.cc`
  (`printv`/`printiv`/`printm`/`printim`; the `#ifdef DEBUG` demo was rewritten too)
- `SPTK/src/utils/sptk_utils.cc:510` (`std::cerr << stream.str()` → `REprintf("%s", …)`)
- tandem: `tandem_64/feature.cpp`, `tandem_64/mScaleInten.cpp` (echo-gated progress prints)

Seven lines that the scripted pass had changed inside comments or `#if 0` bodies were
reverted (`componentManager.cpp`, `vadV1.cpp`, `pitchDirection.cpp`), so the diffs stay
message-only. `R_ext/Print.h` was added to each touched file; in `jkGetF0.cc` the include
had to move out of the file's leading `#if 0` block.

The two fork submodules are still at `1ddf910` (SPTK) and `817652e` (tandem) with these
changes in the working tree: **upstreaming them and bumping the pointers is the one open
task.**

## Phase 3 — formant PRNG

`src/snack_formant.cc` gained an anonymous-namespace MINSTD generator
(`state = state * 16807 mod 2^31-1`, seeded to 1) used by `frand_val()` (`:406`) and the
Newton restart (`:612-613`), replacing `rand()`/`RAND_MAX`.

Equivalence proof: the HEAD and working-tree copies of `snack_formant.cc` were compiled
into two standalone harnesses (`/tmp/fab/main.cpp`) driving `snack_formant::compute_formants()`
over a 2 s synthetic vowel with a moving f0, for `lpc_type` × `window_type` = 3 × 3, and
their outputs compared with `%.17g` precision: **identical** (0 diff lines). MINSTD with
state 1 reproduces macOS `rand()`'s sequence (`16807, 282475249, 1622650073, …`), which is
why. glibc and MSVC use different sequences, so there the sub-1e-6 dither changes.

## Phase 4 — openSMILE

- New `src/smileutil/smileConsole.{h,cpp}` and `src/smileutil/smileRandom.{h,cpp}` under
  `src/opensmile/`, registered in `opensmile_SOURCES`.
  `smileConsole` formats messages itself (stack buffer, heap fallback) and hands finished
  text to writers installed by the host; with none installed, output is discarded, so
  `libopensmile.a` keeps no host dependency and no banned symbol. C++ code uses
  `smile_console_out()`/`smile_console_err()`, stream formatting preserved via a small
  `streambuf`.
- Hosts: `src/opensmile_wrapper.cpp` installs `Rprintf`/`REvprintf` writers (chunked at
  4 kB, and gated to the installing thread — R's console API is not thread-safe, and
  openSMILE runs worker threads); `progsrc/smilextract/SMILExtract.cpp` installs
  `stdout`/`stderr` writers plus `isatty()` for colour, preserving CLI behaviour.
- 18 files converted: `dataPrintSink.{cpp,hpp}` (including the `fputs(text, stdout)`
  helper), `smileUtil.c` (23 sites), `smileUtilSpline.c`, `smileLogger.cpp`
  (`isatty(fileno(stderr))` → `smile_console_is_tty()`), `componentManager.cpp`,
  `vadV1.cpp`, `rnnVad2.cpp`, `mfcc.cpp`, `functionalSegments.cpp`, `functionalPeaks2.cpp`,
  `pitchDirection.cpp`, `exampleSink.cpp`, `vectorOperation.cpp`, `maxIndex.cpp`,
  `signalGenerator.cpp`, `newmat/newmatnl.cpp`, `newmat/myexcept.cpp`.
- `myexcept.cpp`'s four `exit(1)` calls (one in `Terminate()`, three in `DO_FREE_CHECK`
  routines) now raise exceptions and report through `smile_console_err()`.
- Left untouched, with reason: `fprintf`/`fputs` to real files (not banned), print sites
  inside `#ifdef DM_DEBUG_LOGGER` / `BUILD_RNN` / `BUILD_MODELCRYPT` / `DO_REPORT` blocks
  that no build configuration defines, and `fftsg.c`'s `exit(1)` inside the unbuilt
  `cdft_thread_create` macro. `SMILExtract.cpp`'s own `printf`/`puts` stay: the executable
  is not scanned by the check (it globs `libs/**` only) and needs no R runtime.
- `smileRandom` backs `signalGenerator`, `vectorOperation` and `maxIndex`; no bundled
  config in `inst/opensmile/config` instantiates those components, and `maxIndex`'s
  `rand()` is multiplied by a `0.0` default.

## Phase 5 — checks

- `tools/check_cran_symbols.R` on the installed library: **exit 0**, both on the
  `/tmp/superassplib` install and on the library `R CMD check` built.
- `testthat` suite: **0 failures** (`PASS 2724`, 8 warnings and 40 skips as before),
  both standalone and inside the check (`checking tests ... [361s/342s] OK`).
- `SMILExtract -h` prints its banner through the new writers; `configure` builds both the
  static library and the executable; `lst_cpp_opensmile_*` examples and tests pass.
- `R CMD check --as-cran superassp_3.0.0.tar.gz` (Apple clang 21, R 4.6.1):

  ```
  * checking whether package ‘superassp’ can be installed ... [198s/170s] OK
  * checking compiled code ... OK
  * checking examples ... OK
  * checking examples with --run-donttest ... OK
  * checking tests ... [361s/342s] OK
  * checking re-building of vignette outputs ... OK
  Status: 3 NOTEs
  ```

  The three remaining NOTEs are the two carried over from 3.0.0 (CRAN incoming
  feasibility; the `unlockBinding()` call) plus the HTML-manual/tidy note — i.e. the
  compiled-code NOTE is gone and nothing else regressed.
- win-builder (R-devel and R-release) still to be run: the Windows/Linux symbol set is a
  superset of macOS's (`puts`, `putchar` from `printf` folding, `std::cout` as
  `_ZSt4cout`), which the preprocessor-exact source sweep covers.

### One regression found and fixed during the check

The first `--as-cran` run reported `checking whether package ... can be installed ...
WARNING` with exactly two significant warnings: `smile_console_out`/`smile_console_err`
"has C-linkage specified, but returns user-defined type `std::ostream &`". Cause: in
`opensmile_wrapper.cpp` the new include had been placed *inside* the existing
`extern "C" { #include <smileapi/SMILEapi.h> }` block, giving the C++ stream declarations C
linkage. Moving the include below the block cleared it; the re-run is the log above. The
same run also showed unrelated warnings in `00install.out` (rapidjson's
`-Wnontrivial-memcall`, assert-orphaned variables in `assp/*.c`, `dataobj.c`,
`smileUtil.c`) which R does not count as "significant" and which come from files this
amendment does not touch.

## Residual, accepted

| Item | Reason |
|---|---|
| Fork patches live in the working trees, not upstream | Requires a push to `humlab-speech/SPTK` and `humlab-speech/tandem` plus a pointer bump in this repo. SPTK (17 files): `src/utils/sptk_utils.cc`, `third_party/REAPER/{core/float_matrix-inl.h,core/track.cc,core/track.h,epoch_tracker/epoch_tracker.cc,epoch_tracker/fd_filter.cc,wave/codec_api-inl.h,wave/codec_riff.cc,wave/wave.cc,wave/wave.h,wave/wave_io-inl.h,wave/wave_io.cc}`, `third_party/SWIPE/{swipe.cc,vector.cc}`, `third_party/Snack/{jkGetF0.cc,jkGetF0.h,sigproc.cc}`. tandem (2 files): `tandem_64/feature.cpp`, `tandem_64/mScaleInten.cpp`. |
| Disabled-configuration print sites (`DM_DEBUG_LOGGER`, `BUILD_RNN`, `BUILD_MODELCRYPT`, `DO_REPORT`, `#if 0` demos) | Never compiled by any configuration this package builds |
| `fftsg.c` `exit(1)` in the unbuilt `cdft_thread_create` macro | Threaded CDFT path is not enabled |

