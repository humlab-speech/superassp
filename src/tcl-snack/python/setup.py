"""
Build the snack Python C extension.

Searches for a Tcl installation in order:
  1. TCL_INCLUDE_DIR / TCL_LIB_DIR environment variables
  2. pkg-config (tcl)
  3. tclsh / tclsh8.6 introspection
  4. Platform-specific well-known paths
"""

import os
import subprocess
import sys
from setuptools import Extension, setup

GENERIC = os.path.join("..", "generic")
UNIX    = os.path.join("..", "unix")
MAC     = os.path.join("..", "mac")
WIN     = os.path.join("..", "win")


def _run(*cmd, **kw):
    try:
        r = subprocess.run(list(cmd), capture_output=True, text=True,
                           timeout=10, **kw)
        return r.stdout.strip() if r.returncode == 0 else None
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return None


def find_tcl():
    """Return (include_dir, lib_dir, lib_name) for Tcl, or raise."""

    # 1. Explicit env override
    env_inc = os.environ.get("TCL_INCLUDE_DIR")
    env_lib = os.environ.get("TCL_LIB_DIR")
    if env_inc and env_lib:
        # Pick the first lib name whose library file exists in the lib dir.
        if sys.platform == "win32":
            candidates = [
                (name, os.path.join(env_lib, f"{name}.lib"))
                for name in ("tcl86", "tcl85", "tcl")
            ]
        else:
            candidates = [
                (name, os.path.join(env_lib, f"lib{name}.so"))
                for name in ("tcl8.6", "tcl86", "tcl8.5", "tcl")
            ] + [
                (name, os.path.join(env_lib, f"lib{name}.a"))
                for name in ("tcl8.6", "tcl86", "tcl8.5", "tcl")
            ]
        for name, libfile in candidates:
            if os.path.isfile(libfile):
                return env_inc, env_lib, name
        return env_inc, env_lib, "tcl"

    # 2. pkg-config
    inc = _run("pkg-config", "--variable=includedir", "tcl")
    lib = _run("pkg-config", "--variable=libdir", "tcl")
    name = _run("pkg-config", "--libs", "--short-errors", "tcl")
    if inc and lib and name:
        libname = name.lstrip("-l").split()[0]
        return inc, lib, libname

    # 3. tclsh introspection
    for tclsh in ("tclsh", "tclsh8.6", "tclsh8.5"):
        libdir = _run(tclsh, "-c",
                      "puts [file normalize [file join [info library] ..]]")
        incdir = _run(tclsh, "-c",
                      "puts [file normalize "
                      "[file join [info library] .. .. include]]")
        if libdir and incdir and os.path.isfile(
                os.path.join(incdir, "tcl.h")):
            return incdir, libdir, "tcl8.6"

    # 4. Platform-specific fallbacks
    if sys.platform == "darwin":
        for prefix in (
            "/opt/homebrew/opt/tcl-tk",
            "/usr/local/opt/tcl-tk",
        ):
            inc_path = os.path.join(prefix, "include")
            lib_path = os.path.join(prefix, "lib")
            if os.path.isfile(os.path.join(inc_path, "tcl.h")):
                return inc_path, lib_path, "tcl8.6"

    if sys.platform.startswith("linux"):
        for inc_path in (
            "/usr/include/tcl8.6",
            "/usr/include/tcl",
            "/usr/local/include",
        ):
            if os.path.isfile(os.path.join(inc_path, "tcl.h")):
                return inc_path, "/usr/lib/x86_64-linux-gnu", "tcl8.6"

    if sys.platform == "win32":
        for prefix in (r"C:\Tcl", r"C:\ActiveTcl"):
            inc_path = os.path.join(prefix, "include")
            lib_path = os.path.join(prefix, "lib")
            if os.path.isfile(os.path.join(inc_path, "tcl.h")):
                return inc_path, lib_path, "tcl86"

    raise RuntimeError(
        "Cannot locate a Tcl installation.\n"
        "Install tcl-dev (Debian/Ubuntu), tcl-tk (Homebrew), or set\n"
        "TCL_INCLUDE_DIR and TCL_LIB_DIR environment variables."
    )


def platform_audio():
    """Return (source_file, defines, libs, link_args) for the target platform."""
    if sys.platform == "darwin":
        return (
            os.path.join(UNIX, "jkAudIO_osx.c"),
            ["MAC_OSX_TCL"],
            [],
            ["-framework", "CoreAudio", "-framework", "AudioUnit",
             "-framework", "AudioToolbox"],
        )
    if sys.platform.startswith("linux"):
        alsa = os.path.join(UNIX, "jkAudIO_alsa.c")
        oss  = os.path.join(UNIX, "jkAudIO_oss.c")
        src  = alsa if os.path.isfile(alsa) else oss
        define = "ALSA" if "alsa" in src else "OSS"
        libs = ["asound"] if define == "ALSA" else []
        return src, [define], libs, []
    if sys.platform == "win32":
        return (
            os.path.join(WIN, "jkAudIO_win.c"),
            ["WIN32", "_WINDOWS"],
            ["winmm"],
            [],
        )
    raise RuntimeError(f"Unsupported platform: {sys.platform}")


try:
    tcl_inc, tcl_lib, tcl_libname = find_tcl()
    audio_src, audio_defines, audio_libs, audio_link = platform_audio()

    sources = [
        "snack/_snack.c",
        os.path.join(GENERIC, "sound.c"),
        os.path.join(GENERIC, "jkSound.c"),
        os.path.join(GENERIC, "jkSoundEngine.c"),
        os.path.join(GENERIC, "jkSoundEdit.c"),
        os.path.join(GENERIC, "jkSoundFile.c"),
        os.path.join(GENERIC, "jkSoundProc.c"),
        os.path.join(GENERIC, "ffa.c"),
        os.path.join(GENERIC, "jkPitchCmd.c"),
        os.path.join(GENERIC, "jkAudio.c"),
        os.path.join(GENERIC, "jkMixer.c"),
        os.path.join(GENERIC, "shape.c"),
        os.path.join(GENERIC, "jkFilter.c"),
        os.path.join(GENERIC, "jkSynthesis.c"),
        os.path.join(GENERIC, "jkFilterIIR.c"),
        os.path.join(GENERIC, "jkGetF0.c"),
        os.path.join(GENERIC, "sigproc.c"),
        os.path.join(GENERIC, "jkFormant.c"),
        os.path.join(GENERIC, "sigproc2.c"),
        os.path.join(GENERIC, "g711.c"),
        os.path.join(GENERIC, "snackStubInit.c"),
        audio_src,
    ]

    ext = Extension(
        "snack._snack",
        sources=sources,
        include_dirs=[tcl_inc, GENERIC],
        library_dirs=[tcl_lib],
        libraries=[tcl_libname] + audio_libs,
        define_macros=[(d, None) for d in audio_defines] + [("TCL_81_API", None)],
        extra_link_args=audio_link,
    )

    ext_modules = [ext]

except RuntimeError as _e:
    import warnings
    warnings.warn(
        f"Snack C extension will not be built: {_e}\n"
        "Install tcl-dev + libasound2-dev (Linux), tcl-tk (Homebrew), or set\n"
        "TCL_INCLUDE_DIR and TCL_LIB_DIR, then rebuild.",
        stacklevel=2,
    )
    ext_modules = []

setup(ext_modules=ext_modules)
