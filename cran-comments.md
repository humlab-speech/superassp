# CRAN submission comments — superassp 3.0.0

## Breaking change in this version

Every exported `trk_*` wrapper now defaults to `toFile = FALSE`, matching the
`lst_*` functions and the documented function contract. 42 of the 64 `trk_*`
wrappers previously defaulted to `toFile = TRUE`, so a bare `trk_acf("f.wav")`
wrote an SSFF file and returned the number of files written; it now returns an
`AsspDataObj` and writes nothing. Callers who want the files add
`toFile = TRUE`.

This is the only user-visible behaviour change in the release, and it is the
reason for the major version bump. It is documented at the top of `NEWS.md`.

## Test environments

* local: macOS (aarch64-apple-darwin), R 4.6.1, Apple clang 21 — `R CMD check
  --as-cran`, with `pladdrr` installed. The `pladdrr`-absent condition is also
  covered: the 3.0.0 remediation ran the full check against a stubbed `pladdrr`
  (see "Optional dependency").
* win-builder: R-release and R-devel (Windows), resubmitted 2026-09-14 —
  `0 errors | 0 warnings | 2 notes`; the compiled-code note is gone and the one
  WARNING from the previous run (six MinGW `-Wformat` warnings from an openSMILE
  message the compiled-code work rewrote) is fixed
* GitHub Actions, `R CMD check --as-cran --no-manual`: ubuntu-latest
  (R-release and R-devel), macos-latest, windows-latest — all four green

## R CMD check results

| environment | result |
|---|---|
| win-builder, R-devel (2026-09-14) | 0 errors \| 0 warnings \| 2 notes |
| win-builder, R-release (2026-09-14) | submitted alongside R-devel, same sources |
| local macOS, R 4.6.1, Apple clang 21 | 0 errors \| 0 warnings \| 3 notes |

The previous win-builder run reported 2 errors and 7 warnings; the submission
before that reported the compiled-code note. All of those are fixed. The two
notes that remain in both environments are the ones below; the third local note
is the `tidy` version note, which does not arise on win-builder.

### Notes

1. **CRAN incoming feasibility** — "New submission", plus
   "Suggests or Enhances not in mainstream repositories: pladdrr". See
   "Optional dependency" below. Win-builder additionally lists three domain
   terms (`AsspDataObj`, `JSTF`, `Praat`) that `inst/WORDLIST` records; the note
   is unavoidable while `pladdrr` is not in a mainstream repository.
2. **Possibly unsafe call** — one `unlockBinding()` in `R/s7_methods.R`. This
   is intrinsic to the design: every exported `lst_*`/`trk_*` function is
   converted into an S7 generic during `.onLoad()`, which requires replacing
   the binding in the package namespace.
3. **HTML version of manual** — locally this only reports that the `tidy` on
   the checking machine is too old to run the validation. The one real finding
   from the previous run (`format_apply_msg.Rd` emitting a literal `<fun>` HTML
   element) was fixed by rewording the title.

## Changes since the previous submission

* **The two errors had a single cause.** `pladdrr` is an optional, GitHub-only
  dependency, but the test suite and two vignettes called its functions
  unconditionally, so the package failed to check wherever `pladdrr` was not
  installed. Every affected test now skips through one shared helper, and the
  two vignette chunks are gated on `requireNamespace("pladdrr")`. The suite
  reports 0 failures with `pladdrr` present and 0 failures with it absent.
* **The 78 "code/documentation mismatch" warnings are gone.** The `.onLoad()`
  conversion to S7 generics produced a `(listOfFiles, ...)` signature while the
  Rd files documented the real parameter list. The generic now reuses the
  original function's formals, so the installed signature matches the
  documentation, and the documented usage stays as informative as before.
* **Build portability.** `-Wno-register`/`-Wno-deprecated-register` were
  removed by deleting the C++-removed `register` keyword from the bundled Snack
  sources instead of suppressing the diagnostic. The macOS SDK include path is
  now supplied by `configure`, so the generated `src/Makevars` contains no GNU
  make constructs, and `src/Makevars.in` contains none either.
* **The compiled-code note is gone.** The bundled Tandem, openSMILE,
  SPTK/REAPER, SPTK/Snack and SPTK/SWIPE sources no longer reference
  `printf`/`puts`/`putchar`, `stdout`/`stderr`, `std::cout`/`std::cerr`,
  `exit` or `rand`/`srand`:
  - diagnostics go through `Rprintf`/`REprintf`;
  - openSMILE's library — whose static archive is also linked into the
    `SMILExtract` executable shipped in `inst/opensmile/bin`, a standalone
    process with no R runtime — formats messages itself and hands them to
    writers that each host installs (R-backed writers in the package,
    `stdout`/`stderr` writers in the executable);
  - newmat's `Terminate()` raises a C++ exception instead of calling `exit()`;
  - the Snack formant tracker's `rand()` calls use a self-contained MINSTD
    generator, which is the generator macOS's `rand()` implements, so formant
    output is bit-identical on macOS (checked over 9 tracker configurations)
    and platform-invariant from here on; on glibc/MSVC, whose `rand()`
    sequences differ, the sub-1e-6 dither changes;
  - `-DNDEBUG` is set explicitly, matching CRAN's builders, which also keeps a
    local build from compiling `Rcpp/r_cast.h`'s `abort()` path.
  `tools/check_cran_symbols.R` runs R's own scan against an installed library
  and fails if any banned entry point reappears.
* `Remotes:` was removed from DESCRIPTION (it is not a CRAN field) and
  `inst/WORDLIST` was added for the domain terms the spell check flagged.
* **The Windows `-Wformat` warning from the previous win-builder run is gone.**
  Rewriting openSMILE's messages to the console shim changed which printf
  archetype GCC validates against, so the six `%zu`-bearing `smilePcm` messages
  were checked with MSVCRT semantics. `smileConsole.h` now mirrors R's own
  `R_PRINTF_FORMAT` (`gnu_printf` for GCC on UCRT, `printf` elsewhere), which is
  what the shim's `vsnprintf` actually accepts. CI installs the GitHub-only
  `pladdrr` `Suggests` explicitly and its R CMD check gate no longer allowlists
  any WARNING, so this class of regression fails the build instead of passing it.
* The `trk_*` `@param toFile` documentation said "Default `TRUE`" and the
  `@usage` blocks showed `toFile = TRUE`; both now match the new default.

## Optional dependency: pladdrr

`pladdrr` (Praat bindings, <https://github.com/humlab-speech/pladdrr>) is in
`Suggests` and is strictly optional. It is not on CRAN, so it cannot be
installed with `install.packages()`; users who want the Praat-backed functions
install it from GitHub:

```r
remotes::install_github("humlab-speech/pladdrr")
```

Every call site is guarded by `pladdrr_available()`, and the functions raise a
single actionable error naming that command when it is absent. The package
installs, loads, builds its vignettes and passes its full test suite with
`pladdrr` not installed, which is how the check was run locally.

There is no separate Praat installation requirement: Praat's C++ sources are
vendored inside `pladdrr` and compiled into it.

## System requirements

`SystemRequirements: C++17, CMake (>= 3.15)`. CMake is needed because the
bundled openSMILE library is built from source during installation; the
configure script fails with an explicit message naming the install command when
CMake is missing.
