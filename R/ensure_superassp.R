#' Ensure All superassp Python Dependencies Are Installed
#'
#' Installs all Python module dependencies for superassp in one unified environment.
#' This is a best-effort installer that attempts to set up ~12 independent Python
#' module groups in a single shared environment. Core numeric packages (numpy, scipy,
#' etc.) must install successfully; optional packages (PyTorch, TensorFlow, etc.)
#' are installed on a best-effort basis.
#'
#' @details
#' superassp includes ~12 independent Python module groups:
#' \itemize{
#'   \item \strong{Core} (hard requirement): numpy, scipy, soundfile, pandas, scikit-learn
#'   \item \strong{Audio}: librosa, SPTK, resampy, spectral features
#'   \item \strong{Voice Quality}: brouhaha-VAD, COVAREP, voxit, voice_analysis (VAT)
#'   \item \strong{Pitch/Formants}: CREPE, pYIN, YIN, deepformants
#'   \item \strong{PyTorch stack}: torch, torchaudio, pyannote.audio (brouhaha)
#'   \item \strong{TensorFlow stack}: tensorflow, phonet
#'   \item \strong{Local modules}: bundled editable installations from inst/python/
#' }
#'
#' `ensure_superassp()` attempts to install all groups, aborting only if core packages
#' fail. Optional packages (PyTorch, TensorFlow, local modules) may fail gracefully without
#' stopping the process.
#'
#' @param envname Name of Python environment to create or use. Default: `"r-speech"`
#' @param method Installation method: \describe{
#'   \item{`"auto"`}{Try conda first, then virtualenv, then system Python (default)}
#'   \item{`"conda"`}{Use conda environment manager only}
#'   \item{`"virtualenv"`}{Use Python virtualenv only}
#' }
#' @param install_numba Logical. Install Numba JIT compiler for optimized functions?
#'   Default: TRUE. Provides 10-20x speedups with no compilation required.
#' @param install_cython Logical. Install Cython for compiled extensions?
#'   Default: TRUE. Provides 15-25x additional speedups (requires C compiler).
#' @param force Logical. Reinstall all packages even if already present?
#'   Default: FALSE.
#' @param verbose Logical. Print installation progress and summary?
#'   Default: TRUE.
#'
#' @return Invisibly returns a named logical vector with pass/fail status for each
#'   package group. Names include: `core`, `audio`, `small`, `pytorch`, `tensorflow`,
#'   `phonet`, `numba`, `cython`, and one per local module (e.g., `voice_analysis`,
#'   `brouhaha`, `voxit`).
#'
#' @section Environment Selection:
#'
#' When `method = "auto"`:
#' 1. Check if conda is available → use conda
#' 2. Check if virtualenv is available → use virtualenv
#' 3. Fall back to system Python
#'
#' Explicitly setting `method = "conda"` or `method = "virtualenv"` skips detection.
#'
#' @section What Gets Installed:
#'
#' **Phase A — Core (hard requirement, aborts on failure)**:
#' numpy, scipy, soundfile, pandas, scikit-learn, joblib
#'
#' **Phase B — Audio/Signal (soft failure)**:
#' librosa, SPTK, onnxruntime (swift-f0), pywt, python_speech_features,
#' resampy, tqdm, six
#'
#' **Phase C — Numba JIT (optional accelerator)**:
#' numba (if `install_numba = TRUE`)
#'
#' **Phase D — Cython (optional, requires C compiler)**:
#' cython (if `install_cython = TRUE` and C compiler detected)
#'
#' **Phase E — Small/Specialized**:
#' lempel_ziv_complexity, nolds, EMD-signal, pyworld, matplotlib
#'
#' **Phase F — PyTorch Stack (soft failure)**:
#' torch~=2.8.0, torchaudio~=2.8.0, pyannote.audio>=3.3.2,<4.0.0 (for brouhaha-VAD)
#'
#' **Phase G — TensorFlow Stack (soft failure)**:
#' tensorflow, tf-keras, phonet (GitHub)
#'
#' **Phase H — Local Module Editable Installs**:
#' voice_analysis_python, covarep_python, voxit, brouhaha-vad, ftrack_tvwlp,
#' legacy_STRAIGHT (from inst/python/)
#'
#' @section Conflict Handling:
#'
#' **numpy version**: Pinned to `>=1.24.0,<2.0` to satisfy legacy_STRAIGHT's upper
#' bound while meeting all other modules' lower bounds (>=1.19-1.24).
#'
#' **PyTorch vs TensorFlow**: Both can be installed in the same environment with
#' careful version selection. If you encounter conflicts, consider:
#' - Running `ensure_superassp()` in a fresh environment with only one framework
#' - Using separate environments for GPU-intensive work
#' - Checking NVIDIA CUDA compatibility documentation
#'
#' @section Error Behavior:
#'
#' The function aborts (via `rlang::abort()`) if:
#' - `reticulate` package is not available
#' - Python >= 3.8 cannot be found
#' - Core numeric packages (numpy, scipy, etc.) fail to install
#'
#' The function continues with warnings if:
#' - Optional packages (PyTorch, TensorFlow, etc.) fail
#' - Cython compilation fails (module still works without it)
#' - Local module editable installs fail
#'
#' @examples
#' \dontrun{
#' # Install everything in a new "r-superassp" conda environment
#' ensure_superassp(method = "conda", verbose = TRUE)
#'
#' # Check status
#' result <- ensure_superassp(verbose = FALSE)
#' print(result)
#' #> core audio small pytorch tensorflow numba cython brouhaha voxit
#' #> TRUE  TRUE  TRUE    TRUE       TRUE TRUE   TRUE       TRUE   TRUE
#'
#' # Reinstall all packages (useful after updating Python)
#' ensure_superassp(force = TRUE)
#'
#' # Minimal install: core + audio only (skip optional Numba/Cython)
#' ensure_superassp(install_numba = FALSE, install_cython = FALSE)
#' }
#'
#' @seealso
#' Environment management: \code{\link[reticulate]{conda_create}},
#' \code{\link[reticulate]{virtualenv_create}}.
#'
#' @export
ensure_superassp <- function(
  envname    = "r-speech",
  method     = c("auto", "conda", "virtualenv"),
  install_numba  = TRUE,
  install_cython = TRUE,
  force      = FALSE,
  verbose    = TRUE
) {

  method <- match.arg(method)

  # ========== STEP 1: Guard reticulate ==========
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    rlang::abort(c(
      "Package 'reticulate' is required for superassp Python integration.",
      "i" = "Install with: install.packages('reticulate')"
    ))
  }

  # ========== STEP 2: Resolve installation method ==========
  resolved_method <- .resolve_method(method)

  if (verbose) {
    cli::cli_h1("superassp Python Environment Setup")
    cli::cli_alert_info("Environment: {envname}")
    cli::cli_alert_info("Method: {resolved_method}")
    cli::cli_alert_info("Numba JIT: {if(install_numba) 'enabled' else 'disabled'}")
    cli::cli_alert_info("Cython: {if(install_cython) 'enabled' else 'disabled'}")
    cat("\n")
  }

  # ========== STEP 3: Check Python version ==========
  py_ver <- tryCatch(
    numeric_version(reticulate::py_config()$version),
    error = function(e) NULL
  )

  if (is.null(py_ver) || py_ver < numeric_version("3.8.0")) {
    rlang::abort(c(
      "Python >= 3.8 is required for superassp Python modules.",
      "i" = if(is.null(py_ver)) "Found: none" else paste("Found:", as.character(py_ver)),
      "i" = "Install Python 3.8+ from https://python.org or via conda",
      "i" = "  conda install python=3.11"
    ))
  }

  if (verbose) {
    cli::cli_alert_success("Python {as.character(py_ver)} detected")
    cat("\n")
  }

  # ========== STEP 4: Define package groups ==========
  CORE_PKGS <- c(
    "numpy>=1.24.0,<2.0",    # STRAIGHT requires <2.0; others >=1.19-1.24
    "scipy>=1.7.0",
    "soundfile>=0.11.0",
    "pandas>=1.3.0",
    "scikit-learn>=1.0.0",
    "joblib>=1.0.0"
  )

  AUDIO_PKGS <- c(
    "librosa>=0.9.0",
    "onnxruntime>=1.12.0",   # swift-f0
    "pysptk>=0.1.0",
    "pywt>=1.1.1",            # voice_analysis, vat
    "python_speech_features", # phonet
    "resampy>=0.2.0",         # covarep optional
    "tqdm",
    "six"
  )

  SMALL_PKGS <- c(
    "lempel_ziv_complexity>=0.2.0",  # voxit
    "nolds>=0.5.0",                  # voice_analysis
    "EMD-signal>=1.3.0",             # voice_analysis
    "pyworld>=0.3.0",                # covarep optional
    "matplotlib>=3.3.0",             # ftrack_tvwlp, STRAIGHT
    "praat-parselmouth>=0.4.0"       # DisVoice
  )

  PYTORCH_PKGS <- c("torch~=2.8.0", "torchaudio~=2.8.0", "pyannote.audio>=3.3.2,<4.0.0")
  TENSORFLOW_PKGS <- c("tensorflow", "tf-keras")
  OPT_PKGS_NUMBA  <- c("numba>=0.55.0")
  OPT_PKGS_CYTHON <- c("cython>=0.29.0")

  PHONET_GITHUB <- "git+https://github.com/jcvasquezc/phonet.git"

  # Local modules to install from inst/python/
  LOCAL_MODULES <- list(
    voice_analysis = "voice_analysis_python",
    covarep        = "covarep_python",
    voxit          = "voxit",
    brouhaha       = "brouhaha-vad",
    ftrack         = "ftrack_tvwlp",
    straight       = "legacy_STRAIGHT"
  )

  # ========== STEP 5: Install phases ==========

  # Phase A — Core (hard requirement)
  if (verbose) {
    cli::cli_h2("Phase A: Core Packages")
  }
  core_ok <- .try_install(CORE_PKGS, envname, resolved_method, verbose = verbose)
  if (!core_ok) {
    rlang::abort(c(
      "!" = "Could not install core numeric packages.",
      "i" = paste("Environment:", envname, "(method:", resolved_method, ")"),
      "i" = "Check internet connection and Python environment setup.",
      "i" = "Try manually: pip install numpy scipy soundfile pandas scikit-learn"
    ))
  }
  if (verbose) {
    cli::cli_alert_success("Core packages installed")
    cat("\n")
  }

  # Phase B — Audio/signal (soft failure)
  if (verbose) {
    cli::cli_h2("Phase B: Audio/Signal Packages")
  }
  audio_ok <- .try_install(AUDIO_PKGS, envname, resolved_method, verbose = verbose)
  if (!audio_ok && verbose) {
    cli::cli_alert_warning("Some audio packages failed (non-critical)")
  } else if (verbose) {
    cli::cli_alert_success("Audio packages installed")
  }
  if (verbose) cat("\n")

  # Phase C — Numba (optional)
  numba_ok <- FALSE
  if (install_numba) {
    if (verbose) {
      cli::cli_h2("Phase C: Numba JIT Compiler")
    }
    numba_ok <- .try_install(OPT_PKGS_NUMBA, envname, resolved_method, verbose = verbose)
    if (!numba_ok && verbose) {
      cli::cli_alert_warning("Numba installation failed (non-critical)")
    } else if (verbose) {
      cli::cli_alert_success("Numba installed")
    }
    if (verbose) cat("\n")
  }

  # Phase D — Cython (optional, needs C compiler)
  cython_ok <- FALSE
  if (install_cython) {
    if (verbose) {
      cli::cli_h2("Phase D: Cython Extensions")
    }
    has_compiler <- .check_c_compiler(verbose)
    if (!has_compiler) {
      if (verbose) {
        cli::cli_alert_info("No C compiler detected; skipping Cython")
        cli::cli_alert_info("Install gcc, clang, or Visual Studio Build Tools for full speed")
      }
    } else {
      cython_ok <- .try_install(OPT_PKGS_CYTHON, envname, resolved_method, verbose = verbose)
      if (!cython_ok && verbose) {
        cli::cli_alert_warning("Cython installation failed (non-critical)")
      } else if (verbose) {
        cli::cli_alert_success("Cython installed")
      }
    }
    if (verbose) cat("\n")
  }

  # Phase E — Small/specialized packages
  if (verbose) {
    cli::cli_h2("Phase E: Small/Specialized Packages")
  }
  small_ok <- .try_install(SMALL_PKGS, envname, resolved_method, verbose = verbose)
  if (!small_ok && verbose) {
    cli::cli_alert_warning("Some small packages failed (non-critical)")
  } else if (verbose) {
    cli::cli_alert_success("Small packages installed")
  }
  if (verbose) cat("\n")

  # Phase F — PyTorch stack
  if (verbose) {
    cli::cli_h2("Phase F: PyTorch Stack (Brouhaha-VAD)")
  }
  pytorch_ok <- .try_install(PYTORCH_PKGS, envname, resolved_method, verbose = verbose)
  if (!pytorch_ok && verbose) {
    cli::cli_alert_warning("PyTorch stack failed (needed for brouhaha-VAD)")
  } else if (verbose) {
    cli::cli_alert_success("PyTorch stack installed")
  }
  if (verbose) cat("\n")

  # Phase G — TensorFlow stack
  if (verbose) {
    cli::cli_h2("Phase G: TensorFlow Stack (phonet)")
  }
  tensorflow_ok <- .try_install(TENSORFLOW_PKGS, envname, resolved_method, verbose = verbose)
  phonet_ok <- FALSE
  if (tensorflow_ok) {
    phonet_ok <- .try_install(PHONET_GITHUB, envname, resolved_method, verbose = verbose)
    if (!phonet_ok && verbose) {
      cli::cli_alert_warning("phonet GitHub install failed (non-critical)")
    } else if (verbose) {
      cli::cli_alert_success("TensorFlow + phonet installed")
    }
  } else if (verbose) {
    cli::cli_alert_warning("TensorFlow stack failed (needed for phonet)")
  }
  if (verbose) cat("\n")

  # Phase H — Local module editable installs
  if (verbose) {
    cli::cli_h2("Phase H: Local Python Modules")
  }

  local_ok <- vapply(names(LOCAL_MODULES), function(nm) {
    path <- system.file("python", LOCAL_MODULES[[nm]], package = "superassp")
    if (!dir.exists(path) || !file.exists(file.path(path, "setup.py"))) {
      if (verbose) {
        cli::cli_alert_warning("{nm}: setup.py not found")
      }
      return(FALSE)
    }

    if (verbose) {
      cli::cli_alert_info("Installing {nm} (editable)...")
    }

    py_exe <- reticulate::py_config()$python
    res <- tryCatch(
      system2(py_exe, c("-m", "pip", "install", "-e", shQuote(path)),
              stdout = FALSE, stderr = FALSE),
      error = function(e) 1L
    )

    ok <- (res == 0L)
    if (ok && verbose) {
      cli::cli_alert_success("{nm} installed")
    } else if (!ok && verbose) {
      cli::cli_alert_warning("{nm} installation failed (non-critical)")
    }
    ok
  }, logical(1))

  if (verbose) cat("\n")

  # ========== STEP 6: Summary report ==========
  if (verbose) {
    cli::cli_h1("Installation Summary")

    .report_status("Core packages", core_ok)
    .report_status("Audio/signal", audio_ok)
    .report_status("Numba JIT", numba_ok)
    .report_status("Cython", cython_ok)
    .report_status("Small packages", small_ok)
    .report_status("PyTorch stack", pytorch_ok)
    .report_status("TensorFlow stack", tensorflow_ok)
    .report_status("phonet", phonet_ok)

    for (nm in names(local_ok)) {
      .report_status(paste("Local module:", nm), local_ok[[nm]])
    }

    cat("\n")
    cli::cli_alert_info("Environment: {envname}")
    cli::cli_alert_info("Method: {resolved_method}")
    cat("\n")

    if (core_ok) {
      cli::cli_alert_success("superassp is ready to use!")
      cat("\n")
      cli::cli_alert_info("Next steps:")
      cli::cli_alert_info("1. Load package: library(superassp)")
      cli::cli_alert_info("2. Check availability: ?ensure_superassp")
      cli::cli_alert_info("3. Use functions: trk_rapt(), lst_vat(), etc.")
    } else {
      cli::cli_alert_danger("Installation failed; core packages unavailable")
    }

    cat("\n")
  }

  # ========== STEP 7: Return value ==========
  result <- c(
    core = core_ok,
    audio = audio_ok,
    numba = numba_ok,
    cython = cython_ok,
    small = small_ok,
    pytorch = pytorch_ok,
    tensorflow = tensorflow_ok,
    phonet = phonet_ok,
    local_ok
  )

  invisible(result)
}


# ============================================================================
# Internal Helpers
# ============================================================================

#' Resolve Installation Method
#'
#' Detects available Python environment managers and returns the best one.
#' Priority: conda > virtualenv > auto
#'
#' @param method User-specified method or "auto"
#'
#' @return Character string: "conda", "virtualenv", or "auto"
#'
#' @keywords internal
.resolve_method <- function(method) {
  if (method != "auto") return(method)

  # Try conda
  tryCatch({
    conda_list <- reticulate::conda_list()
    if (nrow(conda_list) > 0) return("conda")
  }, error = function(e) NULL)

  # Try virtualenv
  tryCatch({
    venv_list <- reticulate::virtualenv_list()
    if (length(venv_list) > 0) return("virtualenv")
  }, error = function(e) NULL)

  # Fallback
  "auto"
}


#' Try Installing Packages with Error Capture
#'
#' Wrapper around reticulate::py_install() that captures errors and returns
#' a logical indicating success.
#'
#' @param packages Character vector of package names (pip-style)
#' @param envname Environment name
#' @param method Installation method
#' @param pip Logical. Use pip? Default: TRUE
#' @param verbose Logical. Print details?
#'
#' @return Logical. TRUE if installation succeeded, FALSE otherwise.
#'
#' @keywords internal
.try_install <- function(packages, envname, method, pip = TRUE, verbose = TRUE) {
  tryCatch({
    reticulate::py_install(
      packages = packages,
      envname = envname,
      method = method,
      pip = pip
    )
    TRUE
  }, error = function(e) {
    if (verbose) {
      cli::cli_alert_warning("Installation failed: {e$message}")
    }
    FALSE
  })
}


#' Check for C Compiler
#'
#' Detects if a C compiler (gcc, clang, or MSVC) is available on the system.
#'
#' @param verbose Logical. Print diagnostic messages?
#'
#' @return Logical. TRUE if a compiler is available, FALSE otherwise.
#'
#' @keywords internal
.check_c_compiler <- function(verbose = TRUE) {
  # Try common compilers
  compilers <- c("gcc", "clang", "cl")
  for (compiler in compilers) {
    result <- tryCatch(
      system2(compiler, "--version", stdout = FALSE, stderr = FALSE),
      error = function(e) 1L,
      warning = function(w) 1L
    )
    if (result == 0) return(TRUE)
  }
  FALSE
}


#' Report Installation Status
#'
#' Prints a single status line using cli formatting.
#'
#' @param label Character. Description of the component
#' @param ok Logical. Did it install successfully?
#'
#' @return Invisible NULL
#'
#' @keywords internal
.report_status <- function(label, ok) {
  if (ok) {
    cli::cli_alert_success("{label}")
  } else {
    cli::cli_alert_warning("{label}")
  }
  invisible(NULL)
}
