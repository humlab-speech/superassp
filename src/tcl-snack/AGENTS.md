# AGENTS.md — Snack Sound Toolkit

This file provides orientation for AI coding agents working in this repository.

## What This Repository Is

**Snack** (v2.2) is a C-based sound processing extension for scripting languages, developed at KTH (Royal Institute of Technology, Sweden). It provides commands to play, record, process, and visualize audio from **Tcl/Tk** and **Python** scripts. The library compiles to a shared library (`libsnack.so` / `libsnack.dll`) that is loaded by the scripting runtime.

License: GPL (due to the MP3 decoding code in `generic/jkFormatMP3.c`; without that file the code is BSD-compatible).

---

## Directory Layout

```
snack/
├── generic/      # Platform-independent C source — the core library
├── unix/         # Unix build system (autoconf/Makefile) + audio drivers + Tcl package files
├── win/          # Windows build files (MSVC .dsp/.vcproj and Mingw configure)
├── mac/          # macOS-specific notes and patches
├── python/       # Python wrapper module (tkSnack.py) + setup.py
├── ext/          # Example of writing a third-party Snack extension (square-wave filter)
├── demos/
│   ├── tcl/      # Tcl/Tk demo scripts
│   └── python/   # Python demo scripts
├── tests/        # Tcl-based test suite (tcltest framework)
└── doc/          # HTML API documentation
```

---

## Key Source Files (`generic/`)

| File | Purpose |
|---|---|
| `snack.c` / `snack.h` | Package entry point and public API surface |
| `snackDecls.h` / `snackStubInit.c` / `snackStubLib.c` | Stubs mechanism for binary compatibility |
| `jkSound.c` / `jkSound.h` | Core **Sound object**: creation, configuration, memory management |
| `jkSoundEdit.c` | Destructive editing: cut, copy, paste, insert, crop |
| `jkSoundFile.c` | File I/O: WAV, AU, AIFF, SND, RAW, and others |
| `jkSoundEngine.c` | Playback and recording engine |
| `jkSoundProc.c` | Signal processing: mix, convert, normalize |
| `jkAudio.c` / `jkAudIO.h` | Audio device abstraction layer |
| `jkCanvWave.c` | Tk canvas **waveform** item |
| `jkCanvSpeg.c` | Tk canvas **spectrogram** item |
| `jkCanvSect.c` | Tk canvas **spectral section** item |
| `jkFilter.c` / `jkFilterIIR.c` | FIR and IIR digital filter framework |
| `jkFormant.c` / `jkFormant.h` | Formant analysis |
| `jkPitchCmd.c` / `jkGetF0.c` / `jkGetF0.h` | Pitch (F0) estimation |
| `jkSynthesis.c` | Speech synthesis primitives |
| `jkMixer.c` | System mixer / volume control |
| `jkFormatMP3.c` / `jkFormatMP3.h` | MP3 decoding (GPL) |
| `SnackOgg.c` | OGG/Vorbis format support (optional) |
| `SphereFile.c` | NIST/Sphere format support (optional) |
| `g711.c` | G.711 μ-law / A-law codec |
| `ffa.c` / `sigproc.c` / `sigproc2.c` | FFT and other signal-processing utilities |
| `shape.c` | Amplitude envelope shaping |
| `sound.c` | Low-level sound data buffer management |

**Platform audio drivers** (in `unix/`):

| File | Platform |
|---|---|
| `jkAudIO_oss.c` | Linux OSS |
| `jkAudIO_alsa.c` | Linux ALSA |
| `jkAudIO_sun.c` | Solaris / SunOS |
| `jkAudIO_hp.c` | HP-UX |
| `jkAudIO_sgi.c` | IRIX / SGI |
| `jkAudIO_osx.c` | macOS CoreAudio |
| `jkAudIO_skel.c` | Skeleton template for new platforms |

---

## Build System

### Unix / Linux / macOS

```sh
cd unix/
./configure              # detects Tcl/Tk and platform audio
make                     # builds libsound.so and libsnack.so
make install             # installs into the Tcl package directory
```

Useful configure flags:
- `--with-tcl=DIR` / `--with-tk=DIR` — explicit Tcl/Tk paths
- `--enable-alsa` — enable native ALSA support (Linux)
- `--with-nist-include=DIR` / `--with-nist-lib=DIR` — NIST/Sphere support
- `--with-ogg-include=DIR` / `--with-ogg-lib=DIR` — OGG/Vorbis support
- `--disable-stubs` — required for Tcl/Tk 8.0.x

### Windows (MSVC)

Open `win/snack.dsw` in MSVC++ 6.0.

### Windows (Mingw)

```sh
cd win/
./configure --with-tcl=$TCL_DIR --with-tk=$TK_DIR --prefix=$INST_DIR
make && make install
```

---

## Running Tests

Tests use the Tcl `tcltest` framework. Run from the Unix build directory after a successful `make`:

```sh
make test
```

Test files live in `tests/` and are named `<feature>.test` (e.g., `fileio.test`, `filter.test`, `play.test`). The entry point is `tests/all.tcl`.

---

## Python Wrapper (`python/`)

- `tkSnack.py` — Pure-Python module; delegates all audio calls to the Snack Tcl library through Tkinter's `tk.eval()`. Install it on the Python path.
- `setup.py` — Distutils installer for the module.
- Requires a matching Snack + Tcl installation already loaded in the Tkinter interpreter.
- Demo scripts are in `demos/python/`.

---

## Writing Extensions (`ext/`)

Third-party commands that operate on Snack `Sound` objects are written in C and registered with `Snack_AddSubCmd()`. The `ext/` directory contains a fully worked example (`square.c`). Link against the stub library (`libsnackstub2.2.a` on Unix, `snackstub22.lib` on Windows) for binary compatibility across Snack versions.

---

## Coding Conventions

- All C source files begin with the standard GPL copyright header (see any file in `generic/`).
- Public API functions use the prefix `Snack_`; internal functions use `jk` (e.g., `jkSound`, `jkAudio`).
- The Tcl/Tk stubs mechanism is used throughout; do not call Tcl or Tk functions directly when a stub variant exists.
- Audio driver files follow the naming convention `jkAudIO_<platform>.c`.
- Memory allocation uses Tcl's `ckalloc`/`ckfree` rather than `malloc`/`free`.

---

## Documentation

HTML documentation is in `doc/`:
- `tcl-man.html` — complete Tcl/Tk API reference
- `python-man.html` — Python API reference
- `SnackLib.html`, `SoundObj.html`, etc. — C library API for extension authors
