#' Voice Analysis Toolbox - Compute 132 Dysphonia Measures
#'
#' This function computes comprehensive voice quality measures from sustained vowel
#' recordings using the Python implementation of the Voice Analysis Toolbox, a faithful
#' reimplementation of the MATLAB toolbox by Athanasios Tsanas.
#'
#' The function computes 132 dysphonia measures across multiple feature categories:
#' \itemize{
#'   \item Jitter measures (22-25): Local, RAP, PPQ5, DDP, etc.
#'   \item Shimmer measures (22-25): Local, APQ3, APQ5, APQ11, DDA, etc.
#'   \item Harmonic-to-Noise measures (4): HNR, NHR at different frequency bands
#'   \item Nonlinear dynamics (3): DFA, RPDE, PPE
#'   \item Glottal measures (9): GNE variants, Glottal Quotient, VFER
#'   \item MFCCs (84): Delta and delta-delta coefficients
#'   \item Wavelet features (~50): Energy and entropy at multiple scales
#'   \item EMD features (6): Based on Empirical Mode Decomposition
#' }
#'
#' The function uses an in-memory processing pipeline: audio is loaded via the av
#' package and passed directly to Python as a numpy array, eliminating disk I/O.
#' This enables processing of any media format supported by FFmpeg (WAV, MP3, MP4, etc.)
#' without intermediate file conversion.
#'
#' @param listOfFiles Character scalar or vector of input file path(s). Can be any
#'   media format supported by av package (WAV, MP3, MP4, MOV, etc.)
#' @param beginTime Numeric scalar or vector of start time(s) in seconds. Default 0.
#'   If vector, must match length of listOfFiles.
#' @param endTime Numeric scalar or vector of end time(s) in seconds. Default NULL
#'   (use full file duration). If vector, must match length of listOfFiles.
#' @param f0_min Minimum fundamental frequency in Hz. Default 50.
#' @param f0_max Maximum fundamental frequency in Hz. Default 500.
#' @param f0_algorithm F0 estimation algorithm. One of "SWIPE" (default) or "PRAAT".
#'   SWIPE is generally faster and more robust.
#' @param use_thesis_mode Logical. If TRUE, use thesis-compliant implementations
#'   (Tsanas 2012) which differ slightly from the MATLAB code. Default FALSE uses
#'   MATLAB-compatible implementations. Affects PPE (semitone vs log), shimmer dB,
#'   and inclusion of AR jitter & NMSP measures.
#' @param n_cores Integer. Number of CPU cores to use for feature computation.
#'   Default 1 (sequential). Values > 1 enable parallel processing within features.
#'   Use \code{voice_analysis_info()$recommended_workers} for optimal performance.
#' @param use_cython Logical. Use Cython-optimized functions if available (faster).
#'   Default TRUE. Falls back to Numba automatically if Cython not available.
#' @param timeout Numeric. Maximum time in seconds for analysis. Default NULL (no timeout).
#' @param verbose Logical. Print progress messages. Default TRUE.
#' @param return_f0 Logical. If TRUE, include F0 contour in output. Default FALSE.
#' @param toFile Logical. If TRUE, write results to JSTF file. Default FALSE.
#' @param explicitExt Character. File extension for output. Default "vat".
#' @param outputDirectory Character. Output directory path. Default NULL (use input directory).
#'
#' @return If \code{toFile=FALSE} (default), a list with the following components:
#'   \describe{
#'     \item{measures}{Named numeric vector of 132 dysphonia measures}
#'     \item{f0}{Numeric vector of F0 contour (only if return_f0=TRUE)}
#'   }
#'   If \code{toFile=TRUE}, invisibly returns the path(s) to the written JSTF file(s).
#'     \item{fs}{Sample rate in Hz}
#'     \item{success}{Logical indicating if analysis succeeded}
#'     \item{error}{Character string with error message (NULL if success=TRUE)}
#'     \item{file}{Original file path}
#'   }
#'
#'   If multiple files are provided, returns a list of results (one per file).
#'
#' @details
#' ## Computational Requirements
#'
#' Voice analysis is computationally intensive. For a typical 3-second sustained
#' vowel recording:
#' - Sequential processing: 5-15 seconds
#' - Parallel processing (8 cores): 2-5 seconds
#' - With Cython optimizations: 2-3x faster
#'
#' ## Feature Categories
#'
#' The 132 measures are organized into feature groups:
#' 1. **Jitter (22-25 measures)**: Fundamental frequency perturbation
#'    - Local jitter, RAP, PPQ5, DDP
#'    - Absolute and relative variants
#'    - Optional: AR-based jitter (thesis mode)
#'
#' 2. **Shimmer (22-25 measures)**: Amplitude perturbation
#'    - Local shimmer, APQ3, APQ5, APQ11, DDA
#'    - dB and percentage variants
#'    - Optional: NMSP (thesis mode)
#'
#' 3. **Harmonic-to-Noise (4 measures)**:
#'    - HNR (0-500 Hz, 0-1500 Hz)
#'    - NHR (0-500 Hz, 0-1500 Hz)
#'
#' 4. **Nonlinear Dynamics (3 measures)**:
#'    - DFA (Detrended Fluctuation Analysis)
#'    - RPDE (Recurrence Period Density Entropy)
#'    - PPE (Pitch Period Entropy)
#'
#' 5. **Glottal Measures (9 measures)**:
#'    - GNE (Glottal-to-Noise Excitation, 6 variants)
#'    - Glottal Quotient (3 variants)
#'
#' 6. **VFER (7 measures)**: Vocal Fold Excitation Ratio
#'
#' 7. **MFCCs (84 measures)**: Mel-Frequency Cepstral Coefficients
#'    - 12 MFCCs with delta and delta-delta (36 total)
#'    - Statistics: mean, std, min, max, range, skewness, kurtosis
#'
#' 8. **Wavelet (~50 measures)**: Multi-scale wavelet decomposition
#'    - Energy and entropy at multiple scales
#'
#' 9. **EMD (6 measures)**: Empirical Mode Decomposition
#'    - Energy and frequency of intrinsic mode functions
#'
#' ## Performance Tips
#'
#' - Install Cython extensions for 2-3x speedup: \code{install_voice_analysis(method="cython")}
#' - Use parallel processing for faster analysis: \code{n_cores=8}
#' - Process multiple files in batch for efficiency
#' - SWIPE algorithm is typically faster than PRAAT for F0 estimation
#'
#' @references
#' Tsanas, A., Little, M., McSharry, P., & Ramig, L. (2011).
#' Nonlinear speech analysis algorithms mapped to a standard metric achieve
#' clinically useful quantification of average Parkinson's disease symptom severity.
#' Journal of the Royal Society Interface, 8(59), 842-855.
#'
#' Tsanas, A. (2012). Accurate telemonitoring of Parkinson's disease symptom severity
#' using nonlinear speech signal processing and statistical machine learning.
#' D.Phil. thesis, University of Oxford.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # First, install the voice_analysis module
#' install_voice_analysis()
#'
#' # Basic usage - single file
#' result <- lst_vat("sustained_vowel.wav")
#' print(names(result$measures))  # Show all 132 measures
#' print(result$measures["DFA"])  # Detrended Fluctuation Analysis
#'
#' # With time windowing - extract 1.0 to 3.0 seconds
#' result <- lst_vat("recording.wav", beginTime = 1.0, endTime = 3.0)
#'
#' # Parallel processing for speed
#' result <- lst_vat("vowel.wav", n_cores = 8)
#'
#' # Different F0 estimation algorithm
#' result <- lst_vat("vowel.wav", f0_algorithm = "PRAAT")
#'
#' # Custom F0 range (e.g., for male voice)
#' result <- lst_vat("male_vowel.wav", f0_min = 75, f0_max = 300)
#'
#' # Return F0 contour for inspection
#' result <- lst_vat("vowel.wav", return_f0 = TRUE)
#' plot(result$f0, type = "l", main = "F0 Contour")
#'
#' # Batch processing - multiple files
#' files <- c("vowel1.wav", "vowel2.wav", "vowel3.wav")
#' results <- lst_vat(files)
#' # Extract DFA from all files
#' dfa_values <- sapply(results, function(r) r$measures["DFA"])
#'
#' # Process video file (extracts audio automatically)
#' result <- lst_vat("interview.mp4", beginTime = 10, endTime = 13)
#'
#' # Thesis mode (for replication of Tsanas 2012)
#' result <- lst_vat("vowel.wav", use_thesis_mode = TRUE)
#'
#' # Write results to JSTF file
#' lst_vat("vowel.wav", toFile = TRUE)  # Creates vowel.vat
#'
#' # Read back and convert to data.frame
#' track <- read_track("vowel.vat")
#' df <- as.data.frame(track)
#' head(df)  # Shows begin_time, end_time, and all 132 measures
#' }
lst_vat <- function(listOfFiles,
                    beginTime = 0.0,
                    endTime = NULL,
                    f0_min = 50,
                    f0_max = 500,
                    f0_algorithm = c("SWIPE", "PRAAT"),
                    use_thesis_mode = FALSE,
                    n_cores = 1L,
                    use_cython = TRUE,
                    timeout = NULL,
                    verbose = TRUE,
                    return_f0 = FALSE,
                    toFile = FALSE,
                    explicitExt = "vat",
                    outputDirectory = NULL) {

  # Validate inputs
  f0_algorithm <- match.arg(f0_algorithm)

  if (!voice_analysis_available()) {
    stop("voice_analysis module not available.\n",
         "Install with: install_voice_analysis()\n",
         "See ?install_voice_analysis for details.")
  }

  # Handle file inputs
  n_files <- length(listOfFiles)

  if (n_files == 0) {
    stop("No input files provided")
  }

  # Normalize paths
  listOfFiles <- normalizePath(path.expand(listOfFiles), mustWork = FALSE)

  # Check files exist
  missing_files <- !file.exists(listOfFiles)
  if (any(missing_files)) {
    stop("File(s) not found: ", paste(listOfFiles[missing_files], collapse = ", "))
  }

  # Handle time parameters
  if (is.null(beginTime)) {
    beginTime <- 0.0
  }

  if (length(beginTime) == 1) {
    beginTime <- rep(beginTime, n_files)
  } else if (length(beginTime) != n_files) {
    stop("beginTime must be scalar or match length of listOfFiles")
  }

  if (!is.null(endTime)) {
    if (length(endTime) == 1) {
      endTime <- rep(endTime, n_files)
    } else if (length(endTime) != n_files) {
      stop("endTime must be scalar or match length of listOfFiles")
    }
  } else {
    endTime <- rep(NULL, n_files)
  }

  # Import voice_analysis R interface
  va <- reticulate::import("voice_analysis.r_interface")

  # Process files
  results <- vector("list", n_files)

  if (verbose && n_files > 1) {
    message("Processing ", n_files, " files...")
  }

  for (i in seq_len(n_files)) {
    file_path <- listOfFiles[i]
    bt <- beginTime[i]
    et <- endTime[[i]]

    if (verbose && n_files > 1) {
      message(sprintf("[%d/%d] %s", i, n_files, basename(file_path)))
    }

    # Load audio with av
    audio_data <- av_load_for_python(
      file_path = file_path,
      start_time = bt,
      end_time = et,
      target_sample_rate = NULL
    )

    # Perform voice analysis via Python
    result <- va$analyze_for_r(
      audio_path_or_array = audio_data$audio_np,
      fs = as.integer(audio_data$sample_rate),
      features = 'all',
      n_cores = as.integer(n_cores),
      verbose = verbose && (n_files == 1),
      use_cython = use_cython,
      timeout = timeout
    )

    # Convert to R list
    result <- as.list(result)

    # Add file path
    result$file <- file_path

    # Remove F0 if not requested
    if (!return_f0 && !is.null(result$f0)) {
      result$f0 <- NULL
    }

    # Handle JSTF file writing
    if (toFile) {
      json_obj <- create_json_track_obj(
        results = result,
        function_name = "lst_vat",
        file_path = file_path,
        sample_rate = audio_data$sample_rate,
        audio_duration = audio_data$duration,
        beginTime = bt,
        endTime = if (!is.null(et) && et > 0) et else audio_data$duration,
        parameters = list(
          f0_min = f0_min,
          f0_max = f0_max,
          f0_algorithm = f0_algorithm,
          use_thesis_mode = use_thesis_mode,
          n_cores = n_cores,
          use_cython = use_cython
        )
      )

      base_name <- tools::file_path_sans_ext(basename(file_path))
      out_dir <- if (is.null(outputDirectory)) dirname(file_path) else outputDirectory
      output_path <- file.path(out_dir, paste0(base_name, ".", explicitExt))

      write_json_track(json_obj, output_path)
      results[[i]] <- output_path
    } else {
      results[[i]] <- result
    }

  }

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

# Set function attributes


attr(lst_vat, "ext") <- "vat"
attr(lst_vat, "outputType") <- "JSTF"
attr(lst_vat, "format") <- "JSON"
attr(lst_vat, "tracks") <- c(
  "jitter", "shimmer", "hnr", "nhr", "dfa", "rpde", "ppe",
  "gne", "glottal_quotient", "vfer", "mfcc", "wavelet", "emd"
)
