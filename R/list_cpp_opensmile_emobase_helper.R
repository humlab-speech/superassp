#' Compute emobase Features via SMILExtract (C++ Implementation)
#'
#' Extracts emobase features using the SMILExtract command-line tool.
#' This approach is used because emobase's frameMode=full is incompatible
#' with the external audio source architecture used by other feature sets.
#'
#' @param file Path to audio file
#' @param beginTime Start time in seconds (default: 0)
#' @param endTime End time in seconds (default: 0 = end of file)
#' @param verbose Print processing information (default: FALSE)
#' @return Named list with 988 emobase features
#' @keywords internal
lst_emobase_cpp <- function(file, beginTime = 0, endTime = 0, verbose = FALSE) {
  
  # Find SMILExtract binary (built by ./configure at install time)
  smile_bin <- system.file("opensmile", "bin", "SMILExtract",
                          package = "superassp")
  if (smile_bin == "" || !file.exists(smile_bin)) {
    smile_bin <- system.file("opensmile", "bin", "SMILExtract.exe",
                            package = "superassp")
  }

  if (smile_bin == "" || !file.exists(smile_bin)) {
    cli::cli_abort("SMILExtract binary not found. Please reinstall the package.")
  }
  
  # Find emobase config
  config_file <- system.file("opensmile", "config", "emobase", "emobase.conf",
                            package = "superassp")
  
  # Development fallback
  if (config_file == "" || !file.exists(config_file)) {
    pkg_root <- system.file(package = "superassp")
    config_file <- file.path(dirname(dirname(pkg_root)), "src",
                             "opensmile", "config", "emobase", "emobase.conf")
  }
  
  if (!file.exists(config_file)) {
    cli::cli_abort("emobase config file not found")
  }
  
  # SMILExtract's wave reader only takes 16-bit PCM RIFF WAV. Whole 16-bit
  # WAV files go in unchanged; anything else (other containers or codecs,
  # 24/32-bit PCM, a time window) is decoded with the package loader -- the
  # INT16 samples every other DSP wrapper sees -- and written to a temporary
  # 16-bit WAV at the original rate and channel count (openSMILE mixes
  # channels down itself, as it would for the original file).
  input_file <- file
  cleanup_temp <- FALSE
  whole_s16_wav <- beginTime == 0 && endTime == 0 &&
    tolower(tools::file_ext(file)) == "wav" &&
    identical(media_info(file)$audio$codec, "pcm_s16le")

  if (!whole_s16_wav) {
    audio_obj <- assp_load_audio_for_dsp(file, begin = beginTime,
                                         end = endTime, framework = "assp")
    attr(audio_obj, "fileInfo") <- c(21L, 2L)   # RIFF WAV, not the source container
    input_file <- tempfile(fileext = ".wav")
    cleanup_temp <- TRUE
    write_ssff(audio_obj, input_file)
  }

  # Create temp output file
  output_arff <- tempfile(fileext = ".arff")
  
  # Build command (use -O for ARFF output, which is default)
  cmd <- sprintf("%s -C %s -I %s -O %s -instname off -l 0",
                shQuote(smile_bin),
                shQuote(config_file),
                shQuote(input_file),
                shQuote(output_arff))
  
  if (verbose) {
    cat("Running SMILExtract...\n")
    cat("Command:", cmd, "\n")
  }
  
  # Run SMILExtract
  exit_code <- system(cmd, ignore.stdout = !verbose, ignore.stderr = !verbose)
  
  # Cleanup temp input if created
  if (cleanup_temp && file.exists(input_file)) {
    unlink(input_file)
  }
  
  if (exit_code != 0) {
    if (file.exists(output_arff)) unlink(output_arff)
    cli::cli_abort("SMILExtract failed with exit code {exit_code}.")
  }
  
  if (!file.exists(output_arff)) {
    cli::cli_abort("SMILExtract did not produce output file")
  }
  
  # Read CSV output
  if (!requireNamespace("readr", quietly = TRUE)) {
    cli::cli_abort("readr package required for CSV parsing")
  }
  
  # Read ARFF format (openSMILE default)
  lines <- readLines(output_arff)
  
  # Find @data line
  data_idx <- which(lines == "@data")
  if (length(data_idx) == 0) {
    unlink(output_arff)
    cli::cli_abort("Invalid ARFF format: @data section not found")
  }
  
  # Extract attribute names
  attr_lines <- lines[grep("^@attribute", lines)]
  feature_names <- sub("^@attribute\\s+(\\S+)\\s+.*$", "\\1", attr_lines)
  
  # Find first non-empty line after @data
  data_line_idx <- data_idx + 1
  while (data_line_idx <= length(lines) && nchar(trimws(lines[data_line_idx])) == 0) {
    data_line_idx <- data_line_idx + 1
  }
  
  if (data_line_idx > length(lines)) {
    unlink(output_arff)
    cli::cli_abort("No data found after @data section")
  }
  
  data_line <- lines[data_line_idx]
  
  # ARFF rows hold one field per @attribute, in order: the instance name
  # ('off' with -instname off), the features, and a trailing 'class' ('?').
  # Pair fields with attributes by position and keep the features, so a
  # name can never slide onto its neighbour's value.
  data_line <- sub('^"[^"]*",', 'name,', data_line)   # quoted instance name
  parts <- strsplit(data_line, ",")[[1]]
  if (length(parts) != length(feature_names)) {
    unlink(output_arff)
    cli::cli_abort(c(
      "Malformed emobase ARFF row.",
      "i" = "{length(feature_names)} attribute{?s} but {length(parts)} field{?s}."
    ))
  }
  keep <- !feature_names %in% c("name", "frameTime", "class")
  result <- as.list(as.numeric(parts[keep]))
  names(result) <- feature_names[keep]

  # Cleanup
  unlink(output_arff)
  
  if (verbose) {
    cat("Extracted", length(result), "features\n")
  }
  
  return(result)
}
