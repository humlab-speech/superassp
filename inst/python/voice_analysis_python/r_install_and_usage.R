# Voice Analysis Toolbox - R Installation and Usage
# 
# This script helps install and use the Voice Analysis Toolbox from R
# via the reticulate package.

#' Install Voice Analysis Toolbox for R
#'
#' @param method Installation method ("auto", "virtualenv", or "conda")
#' @param python Path to Python executable (optional)
#' @param with_cython Install with Cython optimization (requires C compiler)
#' @param envname Name of Python environment to create
#'
#' @examples
#' # Basic installation
#' install_voice_analysis()
#'
#' # With Cython optimization (faster, requires compiler)
#' install_voice_analysis(with_cython = TRUE)
#'
#' # Using specific Python
#' install_voice_analysis(python = "/usr/bin/python3.10")
#'
install_voice_analysis <- function(
    method = "auto",
    python = NULL,
    with_cython = FALSE,
    envname = "voice-analysis"
) {
  if (!require("reticulate")) {
    stop("Please install reticulate: install.packages('reticulate')")
  }
  
  library(reticulate)
  
  # Configure Python
  if (!is.null(python)) {
    use_python(python, required = TRUE)
  }
  
  cat("Installing Voice Analysis Toolbox...\n")
  
  # Core dependencies
  packages <- c(
    "numpy>=1.20.0",
    "scipy>=1.7.0",
    "soundfile>=0.10.0",
    "pysptk>=0.2.0",
    "librosa>=0.9.0",
    "nolds>=0.5.0"
  )
  
  # Performance dependencies
  packages <- c(packages,
    "numba>=0.56.0",
    "joblib>=1.0.0",
    "pandas>=1.3.0"
  )
  
  # Optional dependencies
  packages <- c(packages,
    "PyWavelets>=1.3.0",
    "PyEMD>=0.3.0"
  )
  
  # Install packages
  py_install(
    packages = packages,
    method = method,
    pip = TRUE,
    envname = envname
  )
  
  # Install Cython if requested
  if (with_cython) {
    cat("\nInstalling Cython and compiling extensions...\n")
    cat("Note: This requires a C compiler (gcc/clang on macOS/Linux, Visual Studio on Windows)\n")
    
    py_install(
      packages = "cython>=0.29.0",
      method = method,
      pip = TRUE,
      envname = envname
    )
    
    # TODO: Compile Cython extensions
    # This would require the source code to be available
    cat("Note: For Cython optimization, install from source:\n")
    cat("  pip install --no-binary :all: voice-analysis\n")
  }
  
  # Verify installation
  cat("\nVerifying installation...\n")
  tryCatch({
    va <- import("voice_analysis")
    cat("✓ Voice Analysis Toolbox installed successfully!\n")
    
    # Get system info
    if (py_module_available("voice_analysis.r_interface")) {
      r_int <- import("voice_analysis.r_interface")
      info <- r_int$get_system_info()
      
      cat("\nSystem Information:\n")
      cat(sprintf("  Platform: %s %s\n", info$platform, info$machine))
      cat(sprintf("  CPU cores: %d physical, %d logical\n", 
                  info$cpu_count_physical, info$cpu_count))
      cat(sprintf("  Cython available: %s\n", 
                  ifelse(info$cython_available, "Yes", "No")))
      cat(sprintf("  Numba available: %s\n", 
                  ifelse(info$numba_available, "Yes", "No")))
      cat(sprintf("  Recommended workers: %d\n", info$recommended_workers))
    }
    
  }, error = function(e) {
    cat("✗ Installation verification failed:\n")
    cat(conditionMessage(e), "\n")
  })
  
  invisible(TRUE)
}


#' Analyze single audio file
#'
#' @param filepath Path to audio file
#' @param n_cores Number of CPU cores to use (default: 1)
#' @param verbose Print progress (default: TRUE)
#'
#' @return List with measures, F0, and metadata
#'
#' @examples
#' result <- analyze_voice_file("audio.wav", n_cores = 8)
#' measures <- result$measures
#' f0 <- result$f0
#'
analyze_voice_file <- function(
    filepath,
    n_cores = 1L,
    verbose = TRUE
) {
  if (!require("reticulate")) {
    stop("Please install reticulate: install.packages('reticulate')")
  }
  
  # Import R interface
  r_int <- import("voice_analysis.r_interface")
  
  # Analyze
  result <- r_int$analyze_for_r(
    audio_path_or_array = filepath,
    n_cores = as.integer(n_cores),
    verbose = verbose
  )
  
  if (!result$success) {
    warning(sprintf("Analysis failed: %s", result$error))
  }
  
  return(result)
}


#' Analyze multiple audio files in batch
#'
#' @param file_paths Character vector of file paths
#' @param n_cores Number of CPU cores to use (default: use all available - 1)
#' @param verbose Print progress (default: TRUE)
#' @param chunk_size Process in chunks for memory efficiency (default: 10)
#'
#' @return Data frame with measures for each file
#'
#' @examples
#' files <- c("audio1.wav", "audio2.wav", "audio3.wav")
#' results <- analyze_voice_batch(files, n_cores = 8)
#' print(head(results))
#'
analyze_voice_batch <- function(
    file_paths,
    n_cores = NULL,
    verbose = TRUE,
    chunk_size = 10L
) {
  if (!require("reticulate")) {
    stop("Please install reticulate: install.packages('reticulate')")
  }
  
  # Determine n_cores
  if (is.null(n_cores)) {
    r_int <- import("voice_analysis.r_interface")
    info <- r_int$get_system_info()
    n_cores <- info$recommended_workers
  }
  
  if (verbose) {
    cat(sprintf("Analyzing %d files with %d workers\n", 
                length(file_paths), n_cores))
  }
  
  # Import R interface
  r_int <- import("voice_analysis.r_interface")
  
  # Analyze batch
  results_dict <- r_int$analyze_batch_for_r(
    file_paths = file_paths,
    n_cores = as.integer(n_cores),
    verbose = verbose,
    return_dataframe = TRUE,
    chunk_size = as.integer(chunk_size)
  )
  
  # Convert to R data.frame
  df <- as.data.frame(results_dict, stringsAsFactors = FALSE)
  
  # Report failures
  if (any(!df$success)) {
    n_failed <- sum(!df$success)
    warning(sprintf("%d/%d files failed to process", n_failed, nrow(df)))
    
    if (verbose) {
      failed_files <- df$file[!df$success]
      cat("\nFailed files:\n")
      for (f in failed_files) {
        cat(sprintf("  - %s\n", f))
      }
    }
  }
  
  return(df)
}


#' Get available voice features
#'
#' @return Character vector of feature names
#'
get_available_features <- function() {
  r_int <- import("voice_analysis.r_interface")
  return(r_int$get_available_features())
}


#' Check system optimization status
#'
#' @return List with system information and recommendations
#'
check_voice_analysis_status <- function() {
  r_int <- import("voice_analysis.r_interface")
  info <- r_int$get_system_info()
  
  cat("Voice Analysis Toolbox Status\n")
  cat("==============================\n\n")
  
  cat(sprintf("Platform: %s %s\n", info$platform, info$machine))
  cat(sprintf("CPU cores: %d physical, %d logical\n", 
              info$cpu_count_physical, info$cpu_count))
  cat("\nOptimizations:\n")
  cat(sprintf("  Cython: %s\n", 
              ifelse(info$cython_available, "✓ Enabled", "✗ Not available")))
  cat(sprintf("  Numba:  %s\n", 
              ifelse(info$numba_available, "✓ Enabled", "✗ Not available")))
  
  cat("\nRecommendations:\n")
  cat(sprintf("  Single file: use n_cores = %d\n", 
              min(4, info$recommended_workers)))
  cat(sprintf("  Batch processing: use n_cores = %d\n", 
              info$recommended_workers))
  
  if (!info$cython_available) {
    cat("\n⚠ For best performance, install with Cython:\n")
    cat("  install_voice_analysis(with_cython = TRUE)\n")
  }
  
  invisible(info)
}


# Example usage script
example_usage <- function() {
  cat("Voice Analysis Toolbox - R Examples\n")
  cat("====================================\n\n")
  
  cat("# Installation\n")
  cat("install_voice_analysis()\n\n")
  
  cat("# Check status\n")
  cat("check_voice_analysis_status()\n\n")
  
  cat("# Single file analysis\n")
  cat('result <- analyze_voice_file("audio.wav", n_cores = 8)\n')
  cat('print(result$measures)\n\n')
  
  cat("# Batch processing\n")
  cat('files <- list.files(pattern = "\\\\.wav$", full.names = TRUE)\n')
  cat('df <- analyze_voice_batch(files, n_cores = 30)\n')
  cat('print(head(df))\n\n')
  
  cat("# Export results\n")
  cat('write.csv(df, "voice_analysis_results.csv", row.names = FALSE)\n\n')
}

# Print usage on source
if (!interactive()) {
  example_usage()
}
