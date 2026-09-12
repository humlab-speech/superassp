##' Harvest Pitch Tracking (C++ implementation)
##'
##' @description Extract F0 (fundamental frequency) using the Harvest algorithm
##'   from WORLD vocoder (via SPTK). This is a high-performance C++
##'   implementation that is 2-3x faster than the Python version and requires
##'   no Python dependencies.
##'
##'   Harvest is designed to be robust and accurate for speech analysis, with
##'   good performance even on noisy signals.
##'
##'   All input media formats are supported via the av package, including video
##'   files from which audio will be automatically extracted.
##'
##' @inheritParams trk_acf
##' @param windowShift Frame shift in milliseconds (default: 10.0)
##' @param minF Minimum F0 in Hz (default: 60.0)
##' @param maxF Maximum F0 in Hz (default: 400.0)
##' @param voicing_threshold Voicing threshold (default: 0.1, valid range: 0.02-0.2 for WORLD/Harvest)
##' @param toFile Write results to file (default: TRUE)
##' @param explicitExt Output file extension (default: "f0")
##' @param parallel Logical. Use parallel processing for multiple files. \code{NULL}
##'   (default) enables automatically for 2+ files.
##' @param n_cores Integer. Number of cores for parallel processing. \code{NULL}
##'   (default) uses \code{detectCores() - 1}.
##'
##' @return If toFile=TRUE, returns the number of successfully processed files.
##'   If toFile=FALSE, returns AsspDataObj or list of AsspDataObj objects.
##'
##' @usage trk_pitch_harvest(
##'   listOfFiles,
##'   beginTime = 0,
##'   endTime = 0,
##'   windowShift = 10,
##'   minF = 60,
##'   maxF = 400,
##'   voicing_threshold = 0.1,
##'   toFile = FALSE,
##'   explicitExt = "f0",
##'   outputDirectory = NULL,
##'   verbose = TRUE,
##'   parallel = NULL,
##'   n_cores = NULL
##' )
##' @param beginTime Start time for the extracted portion in seconds. Default: NULL (beginning of signal). Note: uses `beginTime`/`endTime` (seconds) matching DSP function conventions, unlike [read_audio()] which uses `begin`/`end`.
##' @param endTime The end time of the section of the sound files that should be analysed (in seconds). Use 0 for end of file.
##' @param outputDirectory The directory where the slice file should be stored. If not defiled (NULL), the sparse slice file will placed in the same folder as the media file.
##' @param verbose Logical. Show a progress bar (sequential path) or a progress-aware parallel apply (`pbapply`/`pbmcapply`, if installed).
##' @export
##' @examples
##' wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
##'
##' \donttest{
##' f0 <- trk_pitch_harvest(wav, toFile = FALSE, verbose = FALSE)
##' head(as.data.frame(f0))
##'
##' # Custom F0 range
##' trk_pitch_harvest(wav, minF = 75, maxF = 300, toFile = FALSE, verbose = FALSE)
##' }
trk_pitch_harvest <- function(listOfFiles,
                    beginTime = 0.0,
                    endTime = 0.0,
                    windowShift = 10.0,
                    minF = 60.0,
                    maxF = 400.0,
                    voicing_threshold = 0.1,
                    toFile = FALSE,
                    explicitExt = "f0",
                    outputDirectory = NULL,
                    verbose = TRUE,
                    parallel = NULL,
                    n_cores = NULL) {

  # Validate inputs
  if (is.null(listOfFiles) || length(listOfFiles) == 0) {
    cli::cli_abort("No input files specified in {.arg listOfFiles}")
  }

  # Normalize paths
  listOfFiles <- fast_strip_file_protocol(listOfFiles)
  listOfFiles <- normalizePath(path.expand(listOfFiles), mustWork = FALSE)

  # Check file existence
  validate_file_paths(listOfFiles, function_name = "trk_pitch_harvest")

  n_files <- length(listOfFiles)

  # Normalize time parameters
  beginTime <- if (is.null(beginTime)) 0.0 else beginTime
  endTime <- if (is.null(endTime)) 0.0 else endTime

  # Recycle time parameters
  if (length(beginTime) == 1) beginTime <- rep(beginTime, n_files)
  if (length(endTime) == 1) endTime <- rep(endTime, n_files)

  # Setup output directory
  makeOutputDirectory(outputDirectory, FALSE, "trk_pitch_harvest")

  if (verbose) format_apply_msg("trk_pitch_harvest", n_files, beginTime, endTime)

  # Process each file (parallel for 2+ files, sequential otherwise)
  process_single_file <- function(i) {
    file_path <- listOfFiles[i]
    bt <- beginTime[i]
    et <- endTime[i]

    tryCatch({
      audio_obj <- assp_load_audio_for_dsp(
        file_path,
        begin = bt,
        end   = et,
        framework = "raw"
      )

      harvest_result <- harvest_cpp(
        audio_obj = audio_obj,
        minF = minF,
        maxF = maxF,
        windowShift = windowShift,
        voicing_threshold = voicing_threshold,
        verbose = FALSE
      )

      out_obj <- create_f0_asspobj(harvest_result, windowShift)

      if (toFile) {
        out_file <- generate_output_path(file_path, explicitExt, outputDirectory)
        write.AsspDataObj(out_obj, out_file)
        TRUE
      } else {
        out_obj
      }
    }, error = function(e) {
      cli::cli_warn("Error processing {.file {basename(file_path)}}: {conditionMessage(e)}")
      if (toFile) FALSE else NULL
    })
  }

  results <- run_parallel_files(
    n_files = n_files,
    process_single_file = process_single_file,
    parallel = parallel,
    n_cores = n_cores,
    verbose = verbose,
    export_vars = c("listOfFiles", "beginTime", "endTime", "minF", "maxF",
                     "windowShift", "voicing_threshold", "toFile",
                     "explicitExt", "outputDirectory"),
    export_env = environment()
  )

  # Return results
  if (toFile) {
    n_success <- sum(unlist(results), na.rm = TRUE)
    if (verbose) {
      cli::cli_inform("Successfully processed {n_success} of {n_files} file{?s}")
    }
    return(invisible(n_success))
  } else {
    results <- results[!sapply(results, is.null)]
    if (length(results) == 1) {
      return(results[[1]])
    } else {
      return(results)
    }
  }
}

attr(trk_pitch_harvest, "ext") <- "f0"
attr(trk_pitch_harvest, "tracks") <- c("f0")
attr(trk_pitch_harvest, "outputType") <- "SSFF"
attr(trk_pitch_harvest, "nativeFiletypes") <- c("wav", "flac", "mp3", "mp4", "mkv", "avi")  # Via av
attr(trk_pitch_harvest, "suggestCaching") <- FALSE
