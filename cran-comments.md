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
  --as-cran`, run with `pladdrr` unavailable (see "Optional dependency")
* win-builder: R-devel and R-release (Windows) — previous submission
* GitHub Actions, `R CMD check --as-cran --no-manual`: ubuntu-latest
  (R-release and R-devel), macos-latest, windows-latest

## R CMD check results

0 errors | 0 warnings | 4 notes

The previous win-builder run reported 2 errors and 7 warnings. Both errors and
all seven warnings have been fixed; the notes below are the four that remain.

### Notes

1. **CRAN incoming feasibility** — "New submission", plus
   "Suggests or Enhances not in mainstream repositories: pladdrr". See
   "Optional dependency" below.
2. **Possibly unsafe call** — one `unlockBinding()` in `R/s7_methods.R`. This
   is intrinsic to the design: every exported `lst_*`/`trk_*` function is
   converted into an S7 generic during `.onLoad()`, which requires replacing
   the binding in the package namespace.
3. **Compiled code** — the bundled third-party DSP sources (Tandem, openSMILE,
   SPTK/Snack) reference `rand`, `srand`, `printf` and the C++ streams. The
   `rand`/`srand` calls are in the pitch and formant analysis paths; routing
   them through R's RNG would change numeric output, so they are deliberate.
4. **HTML version of manual** — locally this only reports that the `tidy` on
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
* **Vendored compiler diagnostics were fixed at the source, not suppressed:**
  five deprecated `arma::Mat::max(uword&)` calls, an uninitialised member read
  by REAPER's `FloatMatrix` copy constructor, member-initialiser order in three
  openSMILE headers, rapidjson's use of the C++17-deprecated `std::iterator`,
  `sprintf` in the Tandem sources, and two diagnostic-suppressing pragmas.
* `Remotes:` was removed from DESCRIPTION (it is not a CRAN field) and
  `inst/WORDLIST` was added for the domain terms the spell check flagged.
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
