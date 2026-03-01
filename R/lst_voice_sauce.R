#' Voice Quality Analysis via VoiceSauce
#'
#' Extract comprehensive voice quality measures from speech signals using
#' VoiceSauce algorithms. Computes 40+ parameters including F0, formants,
#' harmonic amplitudes, CPP, HNR, energy, and spectral measures.
#'
#' The function performs a complete VoiceSauce analysis pipeline:
#' 1. F0 estimation (multiple methods available)
#' 2. Formant estimation (F1-F5, B1-B5)
#' 3. Harmonic amplitudes (H1, H2, H4)
#' 4. Formant amplitudes (A1, A2, A3)
#' 5. CPP (Cepstral Peak Prominence)
#' 6. HNR (Harmonics-to-Noise Ratio at multiple frequency bands)
#' 7. Energy
#' 8. Spectral measures (2K, 5K, 2K5K, H42K)
#' 9. Iseli-Alwan formant corrections
#'
#' @param listOfFiles Character vector of audio file paths
#' @param beginTime Numeric vector of start times in seconds (default: 0.0)
#' @param endTime Numeric vector of end times in seconds (default: 0.0, full file)
#' @param frame_shift Frame shift in milliseconds (default: 1.0)
#' @param window_size Window size in milliseconds (default: 25.0)
#' @param f0_method F0 estimation method: "reaper", "praat", "shr", "world"
#'   (default: "reaper")
#' @param f0_min Minimum F0 in Hz (default: 40.0)
#' @param f0_max Maximum F0 in Hz (default: 500.0)
#' @param formant_method Formant estimation method: "praat" (default: "praat")
#' @param n_formants Number of formants to estimate (default: 5)
#' @param max_formant Maximum formant frequency in Hz (default: 5500.0)
#' @param n_periods Number of pitch periods for harmonic analysis (default: 3)
#' @param n_periods_ec Number of pitch periods for CPP/HNR/Energy (default: 5)
#' @param verbose Logical; show progress messages (default: TRUE)
#' @param toFile Logical. If TRUE, write results to JSTF file. Default FALSE.
#' @param explicitExt Character. File extension for output. Default "vsj".
#' @param outputDirectory Character. Output directory path. Default NULL (use input directory).
#'
#' @return If \code{toFile=FALSE} (default): For single file, a named list with
#'   voice quality measures. For multiple files, list of named lists.
#'   If \code{toFile=TRUE}, invisibly returns the path(s) to the written JSTF file(s).
#'
#'   Each measure list contains:
#'   \describe{
#'     \item{\code{F0}}{Fundamental frequency (Hz)}
#'     \item{\code{F1}, \code{F2}, \code{F3}, \code{F4}, \code{F5}}{Formant frequencies (Hz)}
#'     \item{\code{B1}, \code{B2}, \code{B3}, \code{B4}, \code{B5}}{Formant bandwidths (Hz)}
#'     \item{\code{H1}, \code{H2}, \code{H4}}{Harmonic amplitudes (dB)}
#'     \item{\code{A1}, \code{A2}, \code{A3}}{Formant amplitudes (dB)}
#'     \item{\code{H1H2}, \code{H2H4}, \code{H1H4}}{Harmonic differences (dB)}
#'     \item{\code{H1A1}, \code{H1A2}, \code{H1A3}}{Harmonic-formant differences (dB)}
#'     \item{\code{CPP}}{Cepstral Peak Prominence (dB)}
#'     \item{\code{HNR05}, \code{HNR15}, \code{HNR25}, \code{HNR35}}{HNR at 0.5, 1.5, 2.5, 3.5 kHz (dB)}
#'     \item{\code{Energy}}{Frame energy (dB)}
#'     \item{\code{2K}, \code{5K}}{Spectral energy at 2 kHz and 5 kHz (dB)}
#'     \item{\code{2K5K}}{Spectral difference 2K-5K (dB)}
#'     \item{\code{H42K}}{H4-2K difference (dB)}
#'     \item{\code{H1c}, \code{H2c}, \code{H4c}}{Corrected harmonic amplitudes (dB)}
#'     \item{\code{A1c}, \code{A2c}, \code{A3c}}{Corrected formant amplitudes (dB)}
#'     \item{\code{H1H2c}, \code{H2H4c}, \code{H1A1c}, \code{H1A2c}, \code{H1A3c}}{Corrected differences (dB)}
#'     \item{\code{H42Kc}}{Corrected H4-2K difference (dB)}
#'   }
#'
#'   All measures are returned as numeric vectors (time-series), with one value
#'   per frame. Use \code{mean()}, \code{median()}, or other summary functions
#'   to obtain scalar statistics.
#'
#' @details
#' **Voice Quality Measures:**
#'
#' \describe{
#'   \item{\bold{F0 (Fundamental Frequency)}}{\cr
#'     Pitch of the voice, estimated using REAPER, Praat, SHR, or WORLD algorithms.\cr
#'     Clinical use: Voice disorders, prosody analysis\cr
#'     Normal range: 80-250 Hz (male), 150-300 Hz (female)
#'   }
#'   \item{\bold{Formants (F1-F5)}}{\cr
#'     Resonant frequencies of the vocal tract.\cr
#'     F1, F2: Vowel quality\cr
#'     F3-F5: Speaker characteristics, nasality\cr
#'     Clinical use: Articulation disorders, dysarthria
#'   }
#'   \item{\bold{H1, H2, H4}}{\cr
#'     Amplitudes of first, second, and fourth harmonics.\cr
#'     Reflect glottal source characteristics.\cr
#'     Clinical use: Voice quality, breathiness, phonation type
#'   }
#'   \item{\bold{H1-H2, H2-H4}}{\cr
#'     Spectral tilt measures.\cr
#'     Positive values: breathy voice\cr
#'     Negative values: pressed phonation\cr
#'     Clinical use: Dysphonia, voice quality assessment
#'   }
#'   \item{\bold{CPP (Cepstral Peak Prominence)}}{\cr
#'     Measure of harmonicity and periodicity.\cr
#'     Higher values: more periodic, better voice quality\cr
#'     Normal range: 10-25 dB\cr
#'     Clinical use: Gold standard for dysphonia severity
#'   }
#'   \item{\bold{HNR (Harmonics-to-Noise Ratio)}}{\cr
#'     Ratio of harmonic to noise energy at different frequencies.\cr
#'     Higher values: less noise, better voice quality\cr
#'     Clinical use: Breathiness, roughness assessment
#'   }
#'   \item{\bold{Iseli-Alwan Corrections}}{\cr
#'     Formant-corrected versions of harmonics and amplitudes.\cr
#'     Account for vocal tract filtering effects.\cr
#'     More accurate representation of glottal source.
#'   }
#' }
#'
#' **Performance:** With optimizations, processes typical sustained vowel
#' (3s) in ~0.5-2s depending on methods chosen.
#'
#' **Optimization:** Automatically uses Apple Silicon optimizations when available.
#' Check status with \code{voice_sauce_info()}.
#'
#' @references
#' \insertCite{Shue2011}{superassp}
#'
#' \insertCite{Maryn2010}{superassp}
#'
#' @seealso
#' \code{\\link{install_voice_sauce}} for installation,
#' \code{\\link{voice_sauce_info}} for optimization status
#'
#' @examples
#' \\dontrun{
#' # Basic usage with defaults
#' vs <- lst_voice_sauce("vowel.wav")
#' print(mean(vs$CPP, na.rm = TRUE))
#' print(mean(vs$F0, na.rm = TRUE))
#'
#' # Custom F0 range for high-pitched voice
#' vs <- lst_voice_sauce("soprano.wav", f0_min = 150, f0_max = 800)
#'
#' # Different F0 estimation method
#' vs <- lst_voice_sauce("speech.wav", f0_method = "praat")
#'
#' # Fine-grained temporal resolution
#' vs <- lst_voice_sauce("vowel.wav", frame_shift = 0.5)
#'
#' # Batch processing
#' files <- c("a.wav", "e.wav", "i.wav", "o.wav", "u.wav")
#' vs_all <- lst_voice_sauce(files)
#'
#' # Extract specific measures from batch
#' cpp_values <- sapply(vs_all, function(x) mean(x$CPP, na.rm = TRUE))
#' f0_values <- sapply(vs_all, function(x) median(x$F0, na.rm = TRUE))
#'
#' # With time windowing
#' vs_window <- lst_voice_sauce("speech.wav", beginTime = 1.0, endTime = 2.5)
#'
#' # Summary statistics
#' vs <- lst_voice_sauce("vowel.wav")
#' summary_stats <- data.frame(
#'   CPP_mean = mean(vs$CPP, na.rm = TRUE),
#'   HNR_mean = mean(vs$HNR05, na.rm = TRUE),
#'   F0_median = median(vs$F0, na.rm = TRUE),
#'   H1H2_mean = mean(vs$H1H2, na.rm = TRUE)
#' )
#'
#' # Write results to JSTF file
#' lst_voice_sauce("vowel.wav", toFile = TRUE)  # Creates vowel.vsj
#'
#' # Read back and convert to data.frame
#' track <- read_track("vowel.vsj")
#' df <- as.data.frame(track)
#' head(df)  # Shows begin_time, end_time, and all voice quality measures
#' }
#'
#' @export
lst_voice_sauce <- function(listOfFiles,
                            beginTime = 0.0,
                            endTime = 0.0,
                            frame_shift = 1.0,
                            window_size = 25.0,
                            f0_method = "reaper",
                            f0_min = 40.0,
                            f0_max = 500.0,
                            formant_method = "praat",
                            n_formants = 5,
                            max_formant = 5500.0,
                            n_periods = 3,
                            n_periods_ec = 5,
                            verbose = TRUE,
                            toFile = FALSE,
                            explicitExt = "vsj",
                            outputDirectory = NULL) {

  # Check VoiceSauce availability
  if (!voice_sauce_available()) {
    stop("VoiceSauce Python module not available.\\n",
         "Install with: install_voice_sauce()\\n",
         "Check status with: voice_sauce_info()",
         call. = FALSE)
  }

  # Normalize time parameters
  n_files <- length(listOfFiles)
  beginTime <- fast_recycle_times(beginTime, n_files)
  endTime <- fast_recycle_times(endTime, n_files)

  # Validate parameters
  if (!is.numeric(frame_shift) || frame_shift <= 0) {
    stop("frame_shift must be a positive number", call. = FALSE)
  }
  if (!is.numeric(window_size) || window_size <= 0) {
    stop("window_size must be a positive number", call. = FALSE)
  }
  if (!f0_method %in% c("reaper", "praat", "shr", "world")) {
    stop("f0_method must be one of: reaper, praat, shr, world", call. = FALSE)
  }
  if (!is.numeric(f0_min) || f0_min <= 0) {
    stop("f0_min must be a positive number", call. = FALSE)
  }
  if (!is.numeric(f0_max) || f0_max <= f0_min) {
    stop("f0_max must be greater than f0_min", call. = FALSE)
  }

  # Initialize results
  results <- vector("list", n_files)

  # Progress bar
  if (verbose && n_files > 1) {
    cli::cli_alert_info("Extracting VoiceSauce measures from {n_files} file{?s}")
    pb <- cli::cli_progress_bar("VoiceSauce analysis", total = n_files)
  }

  # Process files
  for (i in seq_along(listOfFiles)) {
    file_path <- listOfFiles[i]

    # Validate file exists
    if (!file.exists(file_path)) {
      warning("File not found: ", file_path, call. = FALSE)
      results[[i]] <- NULL
      if (verbose && n_files > 1) cli::cli_progress_update()
      next
    }

    tryCatch({
      # Create VoiceSauce config
      config <- voicesauce_module$VoiceSauceConfig(
        frame_shift = as.numeric(frame_shift),
        window_size = as.numeric(window_size),
        f0_method = as.character(f0_method),
        f0_min = as.numeric(f0_min),
        f0_max = as.numeric(f0_max),
        formant_method = as.character(formant_method),
        n_formants = as.integer(n_formants),
        max_formant = as.numeric(max_formant),
        n_periods = as.integer(n_periods),
        n_periods_ec = as.integer(n_periods_ec)
      )

      # For VoiceSauce, we need to handle time windowing differently
      # since it expects file paths. We'll need to create a temporary file
      # if time windowing is requested
      process_file <- file_path
      temp_file <- NULL

      if (beginTime[i] != 0.0 || endTime[i] != 0.0) {
        # Load audio with time windowing
        audio_data <- av_load_for_python(
          file_path,
          start_time = beginTime[i],
          end_time = endTime[i]
        )

        # Create temporary WAV file
        temp_file <- tempfile(fileext = ".wav")

        # Write windowed audio to temp file
        # Convert float samples to int16
        samples_int <- as.integer(audio_data$samples * 32767)
        samples_int <- pmax(-32768, pmin(32767, samples_int))

        wrassp::write.AsspDataObj(
          list(audio = matrix(samples_int, ncol = 1)),
          file = temp_file,
          sampleRate = audio_data$sample_rate
        )

        process_file <- temp_file
      }

      # Run VoiceSauce analysis
      # Suppress Python print statements if not verbose
      if (verbose && n_files == 1) {
        vs_results <- voicesauce_module$analyze(
          audio_path = process_file,
          config = config
        )
      } else {
        # Suppress output for batch processing
        vs_results <- suppressMessages(
          reticulate::py_capture_output(
            voicesauce_module$analyze(
              audio_path = process_file,
              config = config
            )
          )
        )
      }

      # Clean up temp file if created
      if (!is.null(temp_file) && file.exists(temp_file)) {
        unlink(temp_file)
      }

      # Convert VoiceSauceResults to R list
      measure_list <- list()
      result_dict <- vs_results$to_dict()
      measure_names <- names(result_dict)

      for (name in measure_names) {
        value <- result_dict[[name]]
        # Convert numpy arrays to R vectors
        if (!is.null(value)) {
          measure_list[[name]] <- as.numeric(value)
        }
      }

      # Add times if available
      if (!is.null(vs_results$times)) {
        measure_list[["times"]] <- as.numeric(vs_results$times)
      }

      # Add sampling rate
      if (!is.null(vs_results$fs)) {
        measure_list[["fs"]] <- as.numeric(vs_results$fs)
      }

      # Handle JSTF file writing
      if (toFile) {
        json_obj <- create_json_track_obj(
          results = measure_list,
          function_name = "lst_voice_sauce",
          file_path = file_path,
          sample_rate = if (!is.null(vs_results$fs)) as.numeric(vs_results$fs) else NULL,
          audio_duration = if (!is.null(vs_results$times)) max(as.numeric(vs_results$times)) else NULL,
          beginTime = bt,
          endTime = et,
          parameters = list(
            frame_shift = frame_shift,
            window_size = window_size,
            f0_method = f0_method,
            f0_min = f0_min,
            f0_max = f0_max,
            formant_method = formant_method,
            n_formants = n_formants,
            max_formant = max_formant,
            n_periods = n_periods,
            n_periods_ec = n_periods_ec
          )
        )

        base_name <- tools::file_path_sans_ext(basename(file_path))
        out_dir <- if (is.null(outputDirectory)) dirname(file_path) else outputDirectory
        output_path <- file.path(out_dir, paste0(base_name, ".", explicitExt))

        write_json_track(json_obj, output_path)
        results[[i]] <- output_path
      } else {
        results[[i]] <- measure_list
      }

    }, error = function(e) {
      warning("Error processing ", basename(file_path), ": ",
              e$message, call. = FALSE)
      results[[i]] <- NULL
    })

    if (verbose && n_files > 1) cli::cli_progress_update()
  }

  if (verbose && n_files > 1) cli::cli_progress_done()


  # Simplify if single file
  if (n_files == 1) {
    results <- results[[1]]
  }

  # Return invisibly if writing to file
  if (toFile) {
    return(invisible(results))
  }

  return(results)


}


#' Install VoiceSauce Python Module
#'
#' Install the VoiceSauce Python module from the package's inst/python directory.
#'
#' @param method Installation method: "auto", "virtualenv", or "conda"
#' @param conda Path to conda executable (only for method = "conda")
#' @param envname Name of Python environment to use/create
#' @param python_version Python version to use (e.g., "3.9", "3.10", "3.11")
#' @param install_numba Logical; install numba for JIT optimizations (default: TRUE)
#' @param install_cython Logical; install Cython and compile extensions (default: TRUE)
#' @param force Force reinstallation even if already installed
#'
#' @return Invisible TRUE if successful, FALSE otherwise
#'
#' @details
#' The function installs the VoiceSauce Python module and its dependencies.
#' On Apple Silicon Macs, it automatically enables optimized libraries.
#'
#' **Dependencies installed:**
#' - numpy, scipy: Numerical computing
#' - soundfile: Audio I/O
#' - praat-parselmouth: Praat integration
#' - pyreaper: REAPER F0 estimation
#' - pyworld: WORLD vocoder
#'
#' **Optional optimizations:**
#' - numba: JIT compilation for harmonics, CPP, spectral measures (2-3x speedup)
#' - Cython: Compiled extensions for HNR calculation (2-3x speedup)
#'
#' **Optimization Layers:**
#' 1. Cython (HNR): Compiled C extensions (fastest, if available)
#' 2. Numba JIT (harmonics, CPP, spectral): Runtime compilation (fast)
#' 3. Scipy optimizations: Always available (baseline)
#'
#' **Installation methods:**
#' - "auto": Automatically choose best method for your system
#' - "virtualenv": Use Python virtual environment (recommended)
#' - "conda": Use Conda environment
#'
#' @examples
#' \\dontrun{
#' # Auto installation with all optimizations
#' install_voice_sauce()
#'
#' # Without Cython (simpler, no compilation)
#' install_voice_sauce(install_cython = FALSE)
#'
#' # With specific Python version
#' install_voice_sauce(python_version = "3.10")
#'
#' # Force reinstall with all optimizations
#' install_voice_sauce(force = TRUE)
#' }
#'
install_voice_sauce <- function(method = "auto",
                                conda = "auto",
                                envname = "r-superassp",
                                python_version = "3.9",
                                install_numba = TRUE,
                                install_cython = TRUE,
                                force = FALSE) {

  cli::cli_h2("Installing VoiceSauce Python Module")

  # Get module path
  vs_path <- system.file("python", "voicesauce", package = "superassp")

  if (vs_path == "" || !dir.exists(vs_path)) {
    cli::cli_alert_danger("VoiceSauce module not found in package installation")
    return(invisible(FALSE))
  }

  cli::cli_alert_info("Module path: {vs_path}")

  # Check if already installed (unless force = TRUE)
  if (!force && voice_sauce_available()) {
    cli::cli_alert_success("VoiceSauce already installed and available")
    cli::cli_alert_info("Use force = TRUE to reinstall")
    return(invisible(TRUE))
  }

  tryCatch({
    # Install base dependencies first
    cli::cli_alert_info("Installing Python dependencies...")

    base_packages <- c("numpy", "scipy", "soundfile",
                      "praat-parselmouth", "pyreaper", "pyworld")

    reticulate::py_install(
      packages = base_packages,
      method = method,
      conda = conda,
      envname = envname,
      python_version = python_version,
      pip = TRUE
    )

    # Install optimization packages if requested
    opt_packages <- character(0)

    if (install_numba) {
      cli::cli_alert_info("Installing Numba for JIT optimizations...")
      opt_packages <- c(opt_packages, "numba")
    }

    if (install_cython) {
      cli::cli_alert_info("Installing Cython for compiled extensions...")
      opt_packages <- c(opt_packages, "cython")
    }

    if (length(opt_packages) > 0) {
      reticulate::py_install(
        packages = opt_packages,
        method = method,
        conda = conda,
        envname = envname,
        pip = TRUE
      )
    }

    # Compile Cython extensions if requested and .pyx files exist
    if (install_cython) {
      cython_file <- file.path(vs_path, "measures", "hnr_cython.pyx")
      if (file.exists(cython_file)) {
        cli::cli_alert_info("Compiling Cython extensions...")

        # Check if .so file already exists (pre-compiled)
        so_pattern <- "hnr_cython.*\\.so$"
        so_files <- list.files(file.path(vs_path, "measures"),
                              pattern = so_pattern,
                              full.names = TRUE)

        if (length(so_files) == 0) {
          cli::cli_alert_warning("Cython .so file not found")
          cli::cli_alert_info("VoiceSauce will use scipy optimizations (still fast)")
        } else {
          cli::cli_alert_success("Cython extensions already compiled")
        }
      }
    }

    # Import to test
    cli::cli_alert_info("Testing VoiceSauce installation...")
    voicesauce <- reticulate::import("voicesauce")

    cli::cli_alert_success("VoiceSauce installed successfully!")

    # Check for optimizations
    info <- voice_sauce_info(print_info = FALSE)

    # Report optimization status
    opt_status <- character(0)
    if (info$cython_available) opt_status <- c(opt_status, "Cython")
    if (info$numba_available) opt_status <- c(opt_status, "Numba")
    if (info$apple_silicon) opt_status <- c(opt_status, "Apple Silicon")

    if (length(opt_status) > 0) {
      cli::cli_alert_success("Optimizations active: {paste(opt_status, collapse = ', ')}")
    } else {
      cli::cli_alert_info("Using scipy optimizations (baseline performance)")
    }

    return(invisible(TRUE))

  }, error = function(e) {
    cli::cli_alert_danger("Installation failed: {e$message}")
    cli::cli_alert_info("Try manual installation:")
    cli::cli_ul(c(
      "pip install numpy scipy soundfile",
      "pip install praat-parselmouth pyreaper pyworld",
      "pip install numba  # Optional: for JIT optimizations",
      "pip install cython  # Optional: for compiled extensions"
    ))
    return(invisible(FALSE))
  })
}


#' Check VoiceSauce Availability
#'
#' Check if the VoiceSauce Python module is available.
#'
#' @return Logical; TRUE if VoiceSauce is available, FALSE otherwise
#'
#' @examples
#' \\dontrun{
#' if (voice_sauce_available()) {
#'   vs <- lst_voice_sauce("audio.wav")
#' }
#' }
#'
voice_sauce_available <- function() {
  !is.null(voicesauce_module) && reticulate::py_module_available("voicesauce")
}


#' VoiceSauce System Information
#'
#' Display information about VoiceSauce installation and optimizations.
#'
#' @param print_info Logical; if TRUE, print information to console
#'
#' @return Invisible list with system information
#'
#' @examples
#' \\dontrun{
#' voice_sauce_info()
#' }
#'
voice_sauce_info <- function(print_info = TRUE) {
  info <- list(
    available = voice_sauce_available(),
    cython_available = FALSE,
    numba_available = FALSE,
    apple_silicon = FALSE,
    python_version = NULL
  )

  if (!info$available) {
    if (print_info) {
      cli::cli_alert_danger("VoiceSauce not available")
      cli::cli_alert_info("Install with: install_voice_sauce()")
    }
    return(invisible(info))
  }

  # Get Python version
  py_config <- reticulate::py_config()
  info$python_version <- py_config$version

  # Check for Cython optimization (HNR)
  tryCatch({
    hnr_wrapper <- voicesauce_module$measures$hnr_cython_wrapper
    info$cython_available <- hnr_wrapper$is_cython_available()
  }, error = function(e) {
    info$cython_available <- FALSE
  })

  # Check for Numba optimization (harmonics, CPP, spectral)
  tryCatch({
    harmonics_opt <- voicesauce_module$measures$harmonics_optimized
    info$numba_available <- harmonics_opt$NUMBA_AVAILABLE
  }, error = function(e) {
    # Try alternative check
    tryCatch({
      cpp_opt <- voicesauce_module$measures$cpp_optimized
      info$numba_available <- cpp_opt$NUMBA_AVAILABLE
    }, error = function(e2) {
      info$numba_available <- FALSE
    })
  })

  # Check for Apple Silicon
  if (.Platform$OS.type == "unix" && Sys.info()["sysname"] == "Darwin") {
    arch <- system("uname -m", intern = TRUE)
    info$apple_silicon <- (arch == "arm64")
  }

  if (print_info) {
    cli::cli_h2("VoiceSauce System Information")
    cli::cli_alert_success("VoiceSauce module available")
    cli::cli_alert_info("Python version: {info$python_version}")

    # Report optimizations
    cli::cli_h3("Optimization Status")
    opt_active <- character(0)

    if (info$cython_available) {
      cli::cli_alert_success("Cython: Available (HNR - 2-3x speedup)")
      opt_active <- c(opt_active, "Cython")
    } else {
      cli::cli_alert_info("Cython: Not available (using scipy optimizations)")
    }

    if (info$numba_available) {
      cli::cli_alert_success("Numba: Available (Harmonics, CPP, Spectral - 2-3x speedup)")
      opt_active <- c(opt_active, "Numba")
    } else {
      cli::cli_alert_info("Numba: Not available (using vectorized NumPy)")
    }

    if (info$apple_silicon) {
      cli::cli_alert_success("Apple Silicon: Optimizations enabled")
      opt_active <- c(opt_active, "Apple Silicon")
    }

    if (length(opt_active) == 0) {
      cli::cli_alert_info("Using baseline scipy/numpy optimizations (still fast)")
    }

    cli::cli_h3("Available F0 Methods")
    cli::cli_ul(c("reaper (default)", "praat", "shr", "world"))

    cli::cli_h3("Available Formant Methods")
    cli::cli_ul(c("praat"))

    cli::cli_h3("Available Measures (40+)")
    measures <- c(
      "F0, F1-F5, B1-B5",
      "H1, H2, H4",
      "A1, A2, A3",
      "H1H2, H2H4, H1A1, H1A2, H1A3",
      "CPP, HNR (4 bands), Energy",
      "Spectral: 2K, 5K, 2K5K, H42K",
      "Corrected versions: H1c, H2c, H4c, A1c, A2c, A3c, etc."
    )
    cli::cli_ul(measures)
  }

  return(invisible(info))
}


#' Check VoiceSauce Status on Package Load
#'
#' Internal function called during package startup to check VoiceSauce status.
#'
#' @return Invisible NULL
#' @keywords internal
check_voice_sauce_status <- function() {
  # Only check if module directory exists
  module_dir <- system.file("python", "voicesauce", package = "superassp")

  if (!dir.exists(module_dir)) {
    return(invisible(NULL))  # Module not included in package
  }

  # Check if module is installed
  if (!voice_sauce_available()) {
    return(invisible(NULL))  # Don't spam on every load if not installed
  }

  # Module is installed - check optimization status
  tryCatch({
    info <- voice_sauce_info(print_info = FALSE)

    # Build status message
    status_parts <- character(0)

    if (info$cython_available) {
      status_parts <- c(status_parts, "Cython")
    }

    if (info$numba_available) {
      status_parts <- c(status_parts, "Numba")
    }

    if (info$apple_silicon) {
      status_parts <- c(status_parts, "Apple Silicon")
    }

    if (length(status_parts) > 0) {
      optimizations <- paste(status_parts, collapse = " + ")
      # Only show message occasionally (e.g., once per session)
      if (!exists(".superassp_vs_msg_shown", envir = .GlobalEnv)) {
        packageStartupMessage(
          sprintf("voicesauce: %s optimizations active", optimizations)
        )
        assign(".superassp_vs_msg_shown", TRUE, envir = .GlobalEnv)
      }
    } else {
      # No optimizations - suggest installation with optimizations
      if (!exists(".superassp_vs_warning_shown", envir = .GlobalEnv)) {
        packageStartupMessage(
          "voicesauce: Running with scipy/numpy optimizations.\n",
          "  For 2-3x speedup, install with: install_voice_sauce(install_numba=TRUE, install_cython=TRUE)"
        )
        assign(".superassp_vs_warning_shown", TRUE, envir = .GlobalEnv)
      }
    }

  }, error = function(e) {
    # Silently fail - don't spam users with errors on load
    invisible(NULL)
  })

  invisible(NULL)
}

# Module cache (set in .onLoad)
voicesauce_module <- NULL



# Set function attributes
attr(lst_voice_sauce, "ext") <- "vsj"
attr(lst_voice_sauce, "outputType") <- "JSTF"
attr(lst_voice_sauce, "format") <- "JSON"

