# CRAN submission comments — superassp 2.9.5

## Test environments

* local: macOS (aarch64-apple-darwin23), R 4.6.1 — `R CMD check --as-cran`, PDF manual built
* win-builder: R-release and R-devel (Windows) — see results below
* GitHub Actions, `R CMD check --as-cran --no-manual`: ubuntu-latest (R-release and
  R-devel), macos-latest, windows-latest — all passing
* test suite: ~3,200 assertions, 0 failures

## R CMD check results

0 errors | 7 warnings | 4 notes

The warnings all come from the bundled third-party DSP sources and from one
deliberate documentation choice; none indicate a defect in superassp's own code.

* **"checking whether package 'superassp' can be installed"** — the package
  compiles libassp, ESTK, SPTK (incl. REAPER/WORLD/Snack), tandem, pyin/YIN and
  openSMILE from source, which exceeds the check's default time budget on shared
  machines. It builds in ~3 minutes on a dedicated machine.
* **"checking for code/documentation mismatches"** — intentional. Exported
  `trk_*`/`lst_*` functions are converted into S7 generics at load time, so their
  installed signature is always `(listOfFiles, ...)`. The Rd files deliberately
  document the full, pre-conversion parameter list (via explicit `@usage`) so
  users see the real arguments. Comparing docs to the post-conversion generic
  will always disagree.
* **"checking for GNU extensions in Makefiles"**, **"checking compilation flags
  in Makevars"**, **"checking compilation flags used"**, **"checking pragmas in
  C/C++ headers and code"**, **"checking compiled code"** — all originate in the
  vendored third-party trees (openSMILE, SPTK, tandem, libassp), which ship their
  own build systems, diagnostic-suppressing pragmas and `printf`/`rand` family
  calls. `-Wno-register`/`-Wno-deprecated-register` in `PKG_CXXFLAGS` is required:
  bundled code uses the C++-removed `register` keyword, which Apple Clang treats
  as a hard error under `-std=c++17`.

The four notes are: the usual CRAN incoming feasibility note (new submission),
the absence of a recent HTML Tidy on the check machine, one `unlockBinding()`
call that is intrinsic to the S7 generic conversion, and two historical
`NEWS.md` section titles that predate the versioned format.

## Dependencies

`Remotes` has been removed from the DESCRIPTION submitted here; the repository
keeps it so `pak`/`devtools` can install the two optional GitHub dependencies
during development. The only non-CRAN dependency is `pladdrr` (Praat bindings,
<https://github.com/humlab-speech/pladdrr>), which is in `Suggests` and is
strictly optional: every call site is guarded by `pladdrr_available()` and
raises an actionable error when it is absent.

## Size

The tarball is ~17 MB, essentially all vendored third-party DSP sources. They
are bundled rather than downloaded because CRAN policy forbids network access at
install time, and the package is self-contained by design (no external Praat
installation required).
