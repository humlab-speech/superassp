#' STRAIGHT Speech Synthesis (Legacy Vocoder)
#'
#' Synthesizes speech from STRAIGHT parameters (F0, spectral envelope, aperiodicity).
#' Provides high-quality speech resynthesis for vocoding and voice conversion.
#'
#' @param f0 Numeric vector; F0 contour in Hz
#' @param spec Numeric matrix; spectral envelope \[frames x frequencies\]
#' @param ap Numeric matrix (optional); aperiodicity \[frames x frequencies\]
#' @param sample_rate Numeric; sampling rate in Hz (default: 22050)
#' @param frame_shift Numeric; frame shift in milliseconds (default: 1.0)
#' @param output_file Character (optional); path to save synthesized audio
#' @param verbose Logical; print progress messages (default: TRUE)
#'
#' @return Numeric vector of synthesized audio samples (range: \[-1, 1\])
#'
#' @details
#' The STRAIGHT synthesis algorithm reconstructs speech from vocoder parameters
#' using source-filter decomposition with high-quality pitch-adaptive synthesis.
#'
#' **Input Requirements**:
#' - F0: Fundamental frequency contour (use `trk_straight_f0()`)
#' - Spectral envelope: Pitch-adaptive smoothed spectrum (use `trk_straight_spec()`)
#' - Aperiodicity (optional): Noise component (default: all periodic)
#'
#' **Quality**: Perceptually identical to MATLAB STRAIGHT synthesis
#'
#' @examples
#' \dontrun{
#' # Check availability
#' if (!straight_available()) {
#'   install_legacy_straight()
#' }
#'
#' # Full analysis-synthesis cycle
#' wav_file <- system.file("samples", "sustained", "a1.wav", package = "superassp")
#'
#' # Extract parameters
#' f0_data <- trk_straight_f0(wav_file, toFile = FALSE)
#' spec_data <- trk_straight_spec(wav_file, toFile = FALSE)
#'
#' # Synthesize
#' audio_synth <- straight_synth(
#'   f0 = f0_data$f0[,1],
#'   spec = spec_data$spec,
#'   sample_rate = 22050,
#'   output_file = "resynthesized.wav"
#' )
#'
#' # Play synthesized audio
#' av::av_audio_convert(audio_synth, "resynthesized.wav",
#'                      sample_rate = 22050, channels = 1)
#' }
#'
#' @references
#' \insertRef{Kawahara.1999.SpeechCommunication}{superassp}
#'
#' @seealso
#' \code{\link{trk_straight_f0}}, \code{\link{trk_straight_spec}},
#' \code{\link{install_legacy_straight}}
#'
straight_synth <- function(f0, spec, ap = NULL, sample_rate = 22050,
                          frame_shift = 1.0, output_file = NULL,
                          verbose = TRUE) {
  
  .Deprecated(msg = "straight_synth() is deprecated. WORLD C++ vocoder is the recommended replacement.")

  # Check if STRAIGHT is available
  if (!straight_available()) {
    stop("Legacy STRAIGHT not available. Install with: install_legacy_straight()",
         call. = FALSE)
  }

  # Validate inputs
  if (!is.numeric(f0) || !is.numeric(spec)) {
    stop("f0 and spec must be numeric", call. = FALSE)
  }
  
  if (!is.matrix(spec)) {
    stop("spec must be a matrix [frames x frequencies]", call. = FALSE)
  }
  
  n_frames <- length(f0)
  if (nrow(spec) != n_frames) {
    stop(sprintf("Mismatch: f0 has %d frames, spec has %d frames",
                n_frames, nrow(spec)), call. = FALSE)
  }
  
  # Import Python modules
  .setup_straight_path <- get(".setup_straight_path", envir = asNamespace("superassp"))
  .setup_straight_path()
  
  py_synth <- reticulate::import("legacy_STRAIGHT.synthesis")
  np <- reticulate::import("numpy", convert = FALSE)
  
  # Convert to numpy arrays
  f0_np <- np$array(as.numeric(f0), dtype = "float64")
  spec_np <- np$array(t(spec), dtype = "float64")  # Transpose to [freqs x frames]
  
  # Handle aperiodicity
  if (is.null(ap)) {
    # Default: all periodic (zeros)
    ap_np <- np$zeros_like(spec_np)
  } else {
    ap_np <- np$array(t(ap), dtype = "float64")
  }
  
  if (verbose) {
    message("Running STRAIGHT synthesis...")
  }
  
  # Synthesize
  audio_np <- py_synth$exstraightsynth(
    f0 = f0_np,
    spec = spec_np,
    ap = ap_np,
    fs = as.integer(sample_rate),
    frame_shift = frame_shift / 1000.0
  )
  
  # Convert to R vector
  audio_synth <- as.numeric(reticulate::py_to_r(audio_np))
  
  # Save to file if requested
  if (!is.null(output_file)) {
    if (verbose) {
      message("Saving to: ", output_file)
    }
    
    # Convert to int16 PCM
    audio_int16 <- as.integer(audio_synth * 32767)
    audio_int16 <- pmax(pmin(audio_int16, 32767L), -32768L)  # Clip
    
    # Write using soundfile through Python
    sf <- reticulate::import("soundfile")
    sf$write(output_file, audio_np, as.integer(sample_rate))
  }
  
  if (verbose) {
    message(sprintf("Synthesis complete! Duration: %.2f seconds",
                   length(audio_synth) / sample_rate))
  }
  
  return(audio_synth)
}


#' STRAIGHT Analysis-Synthesis Pipeline
#'
#' Complete analysis-synthesis pipeline combining F0 extraction, spectral
#' analysis, and speech synthesis in one function.
#'
#' @param input_file Character; path to input audio file
#' @param output_file Character (optional); path to save resynthesized audio
#' @param f0_floor Numeric; minimum F0 in Hz (default: 40)
#' @param f0_ceil Numeric; maximum F0 in Hz (default: 800)
#' @param fft_size Numeric; FFT size (default: 2048)
#' @param frame_shift Numeric; frame shift in ms (default: 1.0)
#' @param verbose Logical; print progress (default: TRUE)
#'
#' @return List with:
#'   - `audio`: Synthesized audio vector
#'   - `f0`: F0 data (AsspDataObj)
#'   - `spec`: Spectral data (AsspDataObj)
#'   - `sample_rate`: Sample rate
#'
#' @examples
#' \dontrun{
#' # Full pipeline
#' wav_file <- system.file("samples", "sustained", "a1.wav", package = "superassp")
#' result <- straight_pipeline(wav_file, output_file = "resynthesis.wav")
#'
#' # Access components
#' plot(result$f0$f0[,1], type = "l", ylab = "F0 (Hz)")
#' }
#'
straight_pipeline <- function(input_file, output_file = NULL,
                             f0_floor = 40, f0_ceil = 800,
                             fft_size = 2048, frame_shift = 1.0,
                             verbose = TRUE) {
  
  if (verbose) {
    message("=== STRAIGHT Analysis-Synthesis Pipeline ===\n")
  }
  
  # Step 1: F0 extraction
  if (verbose) message("Step 1: F0 extraction...")
  f0_data <- trk_straight_f0(
    input_file,
    f0_floor = f0_floor,
    f0_ceil = f0_ceil,
    frame_shift = frame_shift,
    toFile = FALSE,
    verbose = FALSE
  )
  
  # Step 2: Spectral analysis
  if (verbose) message("Step 2: Spectral analysis...")
  spec_data <- trk_straight_spec(
    input_file,
    f0_floor = f0_floor,
    f0_ceil = f0_ceil,
    fft_size = fft_size,
    frame_shift = frame_shift,
    toFile = FALSE,
    verbose = FALSE
  )
  
  # Step 3: Synthesis
  if (verbose) message("Step 3: Speech synthesis...")
  audio_synth <- straight_synth(
    f0 = f0_data$f0[,1],
    spec = spec_data$spec,
    sample_rate = attr(f0_data, "sampleRate"),
    frame_shift = frame_shift,
    output_file = output_file,
    verbose = FALSE
  )
  
  if (verbose) {
    message("\n✓ Pipeline complete!")
    if (!is.null(output_file)) {
      message("✓ Saved: ", output_file)
    }
  }
  
  return(list(
    audio = audio_synth,
    f0 = f0_data,
    spec = spec_data,
    sample_rate = attr(f0_data, "sampleRate")
  ))
}
