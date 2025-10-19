#' Install voice_analysis Python module
#'
#' This function installs the voice_analysis Python module from the inst/python/voice_analysis_python
#' directory. The module is a faithful reimplementation of the MATLAB Voice Analysis Toolbox
#' for computing 132 dysphonia measures from sustained vowel recordings.
#'
#' The function can install the module either with or without Cython optimizations.
#' Cython optimizations provide 2-3x performance improvements but require a C compiler.
#' If Cython compilation fails, the module falls back to pure Python/Numba implementations.
#'
#' @param method Installation method, one of:
#'   \describe{
#'     \item{"auto"}{Try Cython first, fallback to pure Python (default)}
#'     \item{"cython"}{Force Cython build (requires C compiler)}
#'     \item{"pure"}{Pure Python/Numba only (no compilation)}
#'   }
#' @param envname Name of Python environment to use (default: "r-reticulate")
#' @param pip_options Additional options to pass to pip install
#' @param force_reinstall If TRUE, reinstall even if already installed
#'
#' @return Invisibly returns TRUE on success, throws error on failure
#' @export
#'
#' @examples
#' \dontrun{
#' # Install with automatic method selection
#' install_voice_analysis()
#'
#' # Force Cython build
#' install_voice_analysis(method = "cython")
#'
#' # Pure Python only (no compilation)
#' install_voice_analysis(method = "pure")
#'
#' # Reinstall
#' install_voice_analysis(force_reinstall = TRUE)
#' }
install_voice_analysis <- function(method = c("auto", "cython", "pure"),
                                   envname = "r-reticulate",
                                   pip_options = NULL,
                                   force_reinstall = FALSE) {

  method <- match.arg(method)

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required. Install with: install.packages('reticulate')")
  }

  # Locate the voice_analysis_python directory
  module_dir <- system.file("python", "voice_analysis_python", package = "superassp")

  if (!dir.exists(module_dir)) {
    stop("voice_analysis_python module not found at: ", module_dir)
  }

  # Check if already installed
  if (!force_reinstall && reticulate::py_module_available("voice_analysis")) {
    message("voice_analysis module is already installed.")

    # Check if r_interface is available
    tryCatch({
      va <- reticulate::import("voice_analysis.r_interface")
      message("  - R interface module: available")

      # Get system info
      sys_info <- va$get_system_info()
      message("  - Cython extensions: ", if(sys_info$cython_available) "available" else "not available")
      message("  - Numba support: ", if(sys_info$numba_available) "available" else "not available")
      message("  - Recommended workers: ", sys_info$recommended_workers)

      message("\nUse force_reinstall=TRUE to reinstall.")
      return(invisible(TRUE))
    }, error = function(e) {
      warning("voice_analysis installed but r_interface not found. Reinstalling...")
    })
  }

  # Determine setup file based on method
  if (method == "pure") {
    setup_file <- file.path(module_dir, "setup.py")
    message("Installing voice_analysis (pure Python/Numba mode)...")
  } else {
    # Try Cython setup
    setup_file <- file.path(module_dir, "setup_cython.py")

    if (!file.exists(setup_file)) {
      warning("setup_cython.py not found, falling back to setup.py")
      setup_file <- file.path(module_dir, "setup.py")
      method <- "pure"
    } else {
      message("Installing voice_analysis with Cython optimizations...")
      message("This may take a few minutes to compile...")
    }
  }

  # Build pip install command
  pip_args <- c("install")

  if (force_reinstall) {
    pip_args <- c(pip_args, "--force-reinstall", "--no-deps")
  }

  # Add upgrade flag
  pip_args <- c(pip_args, "--upgrade")

  # Add user-provided options
  if (!is.null(pip_options)) {
    pip_args <- c(pip_args, pip_options)
  }

  # Install in editable mode from local directory
  if (method == "cython") {
    # For Cython, we need to copy setup_cython.py to setup.py temporarily
    # or use SETUPTOOLS_USE_DISTUTILS environment variable
    pip_args <- c(pip_args, "-e", module_dir)

    # Set environment to use setup_cython.py
    old_wd <- getwd()
    on.exit(setwd(old_wd), add = TRUE)

    setwd(module_dir)

    # Create temporary symbolic link or copy
    if (file.exists("setup.py.bak")) {
      file.remove("setup.py.bak")
    }

    if (file.exists("setup.py")) {
      file.rename("setup.py", "setup.py.bak")
    }

    file.copy("setup_cython.py", "setup.py", overwrite = TRUE)

    on.exit({
      if (file.exists("setup.py")) {
        file.remove("setup.py")
      }
      if (file.exists("setup.py.bak")) {
        file.rename("setup.py.bak", "setup.py")
      }
      setwd(old_wd)
    }, add = TRUE)
  } else {
    pip_args <- c(pip_args, "-e", module_dir)
  }

  # Install dependencies first
  message("\nInstalling dependencies...")

  required_packages <- c(
    "numpy>=1.21.0",
    "scipy>=1.7.0",
    "soundfile>=0.10.0",
    "librosa>=0.9.0",
    "pywt>=1.1.1",
    "pysptk>=0.1.0",
    "nolds>=0.5.0",
    "EMD-signal>=1.3.0",
    "numba>=0.54.0",
    "joblib>=1.0.0",
    "pandas>=1.3.0"
  )

  if (method == "cython") {
    required_packages <- c(required_packages, "cython>=0.29.0")
  }

  tryCatch({
    reticulate::py_install(
      packages = required_packages,
      envname = envname,
      pip = TRUE
    )
  }, error = function(e) {
    warning("Some dependencies may have failed to install: ", e$message)
    message("Continuing with module installation...")
  })

  # Install the voice_analysis module
  message("\nInstalling voice_analysis module...")

  tryCatch({
    # Use system pip for editable installs
    python_path <- reticulate::py_config()$python

    install_cmd <- paste(
      shQuote(python_path),
      "-m pip",
      paste(pip_args, collapse = " ")
    )

    message("Running: ", install_cmd)
    result <- system(install_cmd, intern = TRUE)

    if (length(result) > 0) {
      message(paste(result, collapse = "\n"))
    }

  }, error = function(e) {
    if (method == "cython") {
      warning("Cython installation failed: ", e$message)
      message("\nRetrying with pure Python installation...")
      return(install_voice_analysis(
        method = "pure",
        envname = envname,
        pip_options = pip_options,
        force_reinstall = TRUE
      ))
    } else {
      stop("Installation failed: ", e$message)
    }
  })

  # Verify installation
  message("\nVerifying installation...")

  if (!reticulate::py_module_available("voice_analysis")) {
    stop("Installation appeared to succeed but voice_analysis module is not available")
  }

  # Try to import and get system info
  tryCatch({
    va <- reticulate::import("voice_analysis.r_interface")
    sys_info <- va$get_system_info()

    message("\n", strrep("=", 60))
    message("Installation successful!")
    message(strrep("=", 60))
    message("Platform: ", sys_info$platform, " (", sys_info$machine, ")")
    message("CPU cores: ", sys_info$cpu_count_physical, " physical, ",
            sys_info$cpu_count, " logical")
    message("Cython extensions: ", if(sys_info$cython_available) "available" else "not available")
    message("Numba support: ", if(sys_info$numba_available) "available" else "not available")
    message("Recommended workers: ", sys_info$recommended_workers)
    message(strrep("=", 60))
    message("\nYou can now use lst_vat() to analyze voice recordings.")

  }, error = function(e) {
    warning("Module installed but system info could not be retrieved: ", e$message)
  })

  invisible(TRUE)
}


#' Check if voice_analysis module is available
#'
#' @return Logical indicating if the module is available
#' @export
#'
#' @examples
#' if (voice_analysis_available()) {
#'   message("voice_analysis is ready to use")
#' } else {
#'   message("Install with: install_voice_analysis()")
#' }
voice_analysis_available <- function() {
  reticulate::py_module_available("voice_analysis") &&
    reticulate::py_module_available("voice_analysis.r_interface")
}


#' Get voice_analysis system information
#'
#' @return A list with system information and capabilities
#' @export
#'
#' @examples
#' \dontrun{
#' info <- voice_analysis_info()
#' print(info)
#' }
voice_analysis_info <- function() {
  if (!voice_analysis_available()) {
    stop("voice_analysis module not available. Install with: install_voice_analysis()")
  }

  va <- reticulate::import("voice_analysis.r_interface")
  sys_info <- va$get_system_info()

  as.list(sys_info)
}
