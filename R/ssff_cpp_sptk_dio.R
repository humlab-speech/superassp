##' DIO Pitch Tracking (C++ implementation)
##'
##' @description Extract F0 using the DIO algorithm from the WORLD vocoder (via SPTK).
##'   DIO is designed for high-quality pitch extraction for speech synthesis applications.
##'
##' @inheritParams trk_pitch_rapt
##' @param voicing_threshold Voicing threshold (default: 0.1, valid range: 0.02-0.2 for WORLD/DIO)
##' @param parallel Logical. Use parallel processing for multiple files. \code{NULL}
##'   (default) enables automatically for 2+ files.
##' @param n_cores Integer. Number of cores for parallel processing. \code{NULL}
##'   (default) uses \code{detectCores() - 1}.
##'
##' @return If toFile=TRUE, returns the number of successfully processed files.
##'   If toFile=FALSE, returns AsspDataObj or list of AsspDataObj objects.
##'
##' @export
##' @examples
##' wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
##'
##' \donttest{
##' f0 <- trk_pitch_dio(wav, toFile = FALSE, verbose = FALSE)
##' head(as.data.frame(f0))
##'
##' # Custom F0 range
##' trk_pitch_dio(wav, minF = 80, maxF = 350, toFile = FALSE, verbose = FALSE)
##' }
trk_pitch_dio <- function(listOfFiles,
                beginTime = 0.0,
                endTime = 0.0,
                windowShift = 10.0,
                minF = 60.0,
                maxF = 400.0,
                voicing_threshold = 0.1,
                toFile = TRUE,
                explicitExt = "f0",
                outputDirectory = NULL,
                verbose = TRUE,
                parallel = NULL,
                n_cores = NULL) {

  if (is.null(listOfFiles) || length(listOfFiles) == 0) {
    cli::cli_abort("No input files specified in {.arg listOfFiles}")
  }

  listOfFiles <- fast_strip_file_protocol(listOfFiles)
  listOfFiles <- normalizePath(path.expand(listOfFiles), mustWork = FALSE)

  validate_file_paths(listOfFiles, function_name = "trk_pitch_dio")

  n_files <- length(listOfFiles)

  beginTime <- if (is.null(beginTime)) 0.0 else beginTime
  endTime <- if (is.null(endTime)) 0.0 else endTime

  if (length(beginTime) == 1) beginTime <- rep(beginTime, n_files)
  if (length(endTime) == 1) endTime <- rep(endTime, n_files)

  makeOutputDirectory(outputDirectory, FALSE, "trk_pitch_dio")

  if (verbose) format_apply_msg("trk_pitch_dio", n_files, beginTime, endTime)

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

      dio_result <- dio_cpp(
        audio_obj = audio_obj,
        minF = minF,
        maxF = maxF,
        windowShift = windowShift,
        voicing_threshold = voicing_threshold,
        verbose = FALSE
      )

      out_obj <- create_f0_asspobj(dio_result, windowShift)

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

attr(trk_pitch_dio, "ext") <- "f0"
attr(trk_pitch_dio, "tracks") <- c("f0")
attr(trk_pitch_dio, "outputType") <- "SSFF"
attr(trk_pitch_dio, "nativeFiletypes") <- c("wav", "flac", "mp3", "mp4", "mkv", "avi")
attr(trk_pitch_dio, "suggestCaching") <- FALSE
