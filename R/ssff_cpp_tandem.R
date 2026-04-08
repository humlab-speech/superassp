#' TANDEM Pitch Tracking and Voiced Speech Segregation
#'
#' Estimate pitch (F0) and separate voiced speech from background using the TANDEM
#' algorithm, which combines computational auditory scene analysis (CASA) with
#' neural network-based pitch tracking.
#'
#' @details
#' **TANDEM (Tandem Algorithm for pitch estimation and voiced speech segregation)**
#' 
#' This function implements the algorithm described in Hu & Wang (2010), which
#' simultaneously tracks pitch and segments voiced speech from background noise/interference.
#' 
#' **Key Features**:
#' - Gammatone filterbank (64 channels, 50-8000 Hz)
#' - Neural network-based pitch estimation
#' - Time-frequency voiced speech masking
#' - Robust to noise and reverberation
#' - Handles multiple simultaneous pitch sources
#'
#' **Sample Rate**: TANDEM requires 20 kHz sample rate. Audio will be automatically
#' resampled if necessary.
#'
#' **Note**: Currently returns placeholder results while full TANDEM integration is completed.
#' The algorithm framework is in place but core processing is under development.
#'
#' @param listOfFiles Character vector of audio file paths (any format supported by av package)
#' @param minF Numeric, minimum F0 in Hz. Default: 50.
#' @param maxF Numeric, maximum F0 in Hz. Default: 500.
#' @param target_sample_rate Numeric, internal processing sample rate. Default: 20000 Hz (required by TANDEM).
#' @param return_mask Logical, return time-frequency voiced mask. Default: FALSE.
#' @param toFile Logical, write output to SSFF file. Default: FALSE.
#' @param explicitExt Character, output file extension. Default: "tnd".
#' @param outputDirectory Character, output directory. Default: NULL (same as input).
#' @param verbose Logical, print progress messages. Default: TRUE.
#' @param ... Additional arguments (currently unused).
#'
#' @return If toFile=FALSE: AsspDataObj or list of AsspDataObj objects with tracks:
#'   - `pitch`: F0 track in Hz
#'   - `voicing_prob`: Voicing probability (0-1)
#'   
#'   If toFile=TRUE: Character vector of output file paths.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic pitch tracking
#' result <- trk_tandem("speech.wav")
#' plot(result$pitch, type = "l", main = "TANDEM Pitch Track")
#'
#' # With noisy speech
#' result <- trk_tandem("noisy_speech.wav", minF = 80, maxF = 400)
#' 
#' # Batch processing
#' files <- c("speaker1.wav", "speaker2.wav", "speaker3.wav")
#' results <- trk_tandem(files, verbose = TRUE)
#' 
#' # Save to files
#' trk_tandem("speech.wav", toFile = TRUE, outputDirectory = "output/")
#' }
#'
#' @references
#' \insertRef{hu2010tandem}{superassp}
#' 
#' \insertRef{hu2011unvoiced}{superassp}
#'
#' @seealso \code{\link{trk_rapt}}, \code{\link{trk_swipe}}, \code{\link{trk_yin}} 
#'   for other pitch tracking methods
trk_tandem <- function(
  listOfFiles,
  minF = 50,
  maxF = 500,
  target_sample_rate = 20000,
  return_mask = FALSE,
  toFile = FALSE,
  explicitExt = "tnd",
  outputDirectory = NULL,
  verbose = TRUE,
  ...
) {
  # Validate inputs
  if (!is.character(listOfFiles)) {
    stop("listOfFiles must be a character vector")
  }
  
  # Check file existence
  missing_files <- listOfFiles[!file.exists(listOfFiles)]
  if (length(missing_files) > 0) {
    stop("File(s) not found: ", paste(missing_files, collapse = ", "))
  }
  
  n_files <- length(listOfFiles)
  results <- vector("list", n_files)
  
  if (verbose) format_apply_msg("trk_tandem", n_files)
  
  for (i in seq_along(listOfFiles)) {
    if (verbose && n_files > 1) {
      message("  [", i, "/", n_files, "] ", basename(listOfFiles[i]))
    }
    
    # Load audio via av package
    tryCatch({
      audio_data <- av::read_audio_bin(
        listOfFiles[i],
        channels = 1  # TANDEM requires mono
      )
    }, error = function(e) {
      stop("Failed to load audio file: ", listOfFiles[i], "\n", e$message)
    })
    
    orig_sr <- attr(audio_data, "sample_rate")
    audio_vec <- as.numeric(audio_data)
    
    # Resample to 20 kHz if needed (using av package)
    if (orig_sr != target_sample_rate) {
      if (verbose) {
        message("    Resampling from ", orig_sr, " Hz to ", target_sample_rate, " Hz...")
      }
      # Use av package for resampling
      temp_wav <- tempfile(fileext = ".wav")
      on.exit(unlink(temp_wav), add = TRUE)
      
      # Write resampled audio
      av::av_audio_convert(
        listOfFiles[i],
        temp_wav,
        format = "wav",
        sample_rate = target_sample_rate,
        channels = 1
      )
      
      # Re-load resampled audio
      audio_data <- av::read_audio_bin(temp_wav, channels = 1)
      audio_vec <- as.numeric(audio_data)
    }
    
    # TANDEM requires neural network files in "net/" subdirectory
    # Create temporary net/ directory with symlinks
    net_dir <- file.path(getwd(), "net")
    if (!dir.exists(net_dir)) {
      dir.create(net_dir)
      created_net_dir <- TRUE
    } else {
      created_net_dir <- FALSE
    }
    
    # Symlink or copy network files
    net_source <- system.file("tandem_net", package = "superassp")
    for (net_file in c("MLP1.64.dat", "MLP2.64.dat", "MLP3.64.dat")) {
      src <- file.path(net_source, net_file)
      dst <- file.path(net_dir, net_file)
      if (!file.exists(dst) && file.exists(src)) {
        # Try symlink first, fall back to copy
        tryCatch({
          file.symlink(src, dst)
        }, error = function(e) {
          file.copy(src, dst)
        })
      }
    }
    
    # Ensure cleanup on exit
    on.exit({
      if (created_net_dir && dir.exists(net_dir)) {
        unlink(net_dir, recursive = TRUE)
      }
    }, add = TRUE)
    
    # Call TANDEM C++ wrapper
    tandem_result <- tandem_pitch_cpp(
      audio_signal = audio_vec,
      sample_rate = target_sample_rate,
      min_pitch = minF,
      max_pitch = maxF,
      net_path = system.file("tandem_net", package = "superassp")
    )
    
    # Check status
    if ("status" %in% names(tandem_result) && tandem_result$status == "placeholder") {
      if (verbose && i == 1) {
        message("    Status: Using placeholder implementation (full TANDEM integration pending)")
      }
    }
    
    # Convert to AsspDataObj
    assp_obj <- list(
      pitch = tandem_result$pitch,
      voicing_prob = tandem_result$voicing_prob
    )
    
    # Set attributes
    attr(assp_obj, "sampleRate") <- 100  # Analysis rate (100 Hz frames)
    attr(assp_obj, "startTime") <- 0
    attr(assp_obj, "startRecord") <- 1L
    attr(assp_obj, "endRecord") <- as.integer(length(tandem_result$pitch))
    attr(assp_obj, "trackFormats") <- c("REAL64", "REAL64")
    class(assp_obj) <- "AsspDataObj"
    
    # Write to file if requested
    if (toFile) {
      base_name <- tools::file_path_sans_ext(basename(listOfFiles[i]))
      out_dir <- if (is.null(outputDirectory)) {
        dirname(listOfFiles[i])
      } else {
        outputDirectory
      }
      
      # Create output directory if needed
      if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
      }
      
      output_path <- file.path(out_dir, paste0(base_name, ".", explicitExt))
      
      tryCatch({
        wrassp::write.AsspDataObj(assp_obj, output_path)
        results[[i]] <- output_path
      }, error = function(e) {
        stop("Failed to write output file: ", output_path, "\n", e$message)
      })
      
    } else {
      results[[i]] <- assp_obj
    }
  }
  
  # Simplify output for single file
  if (n_files == 1) {
    return(results[[1]])
  } else {
    return(results)
  }
}

# Set function attributes
attr(trk_tandem, "ext") <- "tnd"
attr(trk_tandem, "tracks") <- c("pitch", "voicing_prob")
attr(trk_tandem, "outputType") <- "SSFF"
attr(trk_tandem, "nativeFiletypes") <- c("wav")
