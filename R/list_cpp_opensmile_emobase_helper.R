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
  
  # Find SMILExtract binary
  smile_bin <- system.file("opensmile", "bin", "SMILExtract",
                          package = "superassp")
  
  # If not found in inst, try build location (development)
  if (smile_bin == "" || !file.exists(smile_bin)) {
    # Try relative path from package root
    pkg_root <- system.file(package = "superassp")
    smile_bin <- file.path(dirname(dirname(pkg_root)), "src", 
                           "opensmile", "build_r", "progsrc", 
                           "smilextract", "SMILExtract")
  }
  
  if (!file.exists(smile_bin)) {
    stop("SMILExtract binary not found. Please rebuild package.")
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
    stop("emobase config file not found")
  }
  
  # Handle time windowing if needed
  input_file <- file
  cleanup_temp <- FALSE
  
  if (beginTime > 0 || endTime > 0) {
    # Need to extract audio segment first
    if (!requireNamespace("av", quietly = TRUE)) {
      stop("av package required for time windowing")
    }
    
    # Read audio segment
    audio_data <- av::read_audio_bin(file, start_time = beginTime,
                                     end_time = if(endTime > 0) endTime else NULL,
                                     channels = 1)
    sample_rate <- attr(audio_data, "sample_rate")
    
    # Write to temp WAV file
    temp_wav <- tempfile(fileext = ".wav")
    cleanup_temp <- TRUE
    
    # Convert to 16-bit PCM
    audio_int16 <- as.integer(audio_data)
    
    # Write WAV file
    if (!requireNamespace("wrassp", quietly = TRUE)) {
      stop("wrassp package required for WAV writing")
    }
    
    # Create AsspDataObj
    audio_obj <- list(audio = matrix(audio_int16, ncol = 1))
    attr(audio_obj, "sampleRate") <- sample_rate
    attr(audio_obj, "startTime") <- 0
    attr(audio_obj, "trackFormats") <- "INT16"
    class(audio_obj) <- "AsspDataObj"
    
    wrassp::write.AsspDataObj(audio_obj, file = temp_wav)
    input_file <- temp_wav
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
    stop("SMILExtract failed with exit code ", exit_code)
  }
  
  if (!file.exists(output_arff)) {
    stop("SMILExtract did not produce output file")
  }
  
  # Read CSV output
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop("readr package required for CSV parsing")
  }
  
  # Read ARFF format (openSMILE default)
  lines <- readLines(output_arff)
  
  # Find @data line
  data_idx <- which(lines == "@data")
  if (length(data_idx) == 0) {
    unlink(output_arff)
    stop("Invalid ARFF format: @data section not found")
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
    stop("No data found after @data section")
  }
  
  data_line <- lines[data_line_idx]
  
  # Remove quoted filename at the start (handles filenames with commas)
  data_line <- sub('^"[^"]*",', '', data_line)
  
  # Split by comma
  parts <- strsplit(data_line, ",")[[1]]
  
  # With -instname off, first field is 'off', then features
  # Skip first field
  if (length(parts) < 2) {
    unlink(output_arff)
    stop("Invalid data line format: insufficient fields")
  }
  
  values <- as.numeric(parts[2:length(parts)])
  
  # Remove NA values if any
  values <- values[!is.na(values)]
  
  # Adjust feature names - skip 'name' and 'frameTime' attributes from header
  # (we use -instname off so frameTime becomes the instance name field)
  if (length(feature_names) >= 2) {
    # Skip first two attributes (name, frameTime or similar)
    feature_names <- feature_names[3:length(feature_names)]
  }
  
  # Match lengths
  min_len <- min(length(values), length(feature_names))
  values <- values[1:min_len]
  feature_names <- feature_names[1:min_len]
  
  # Create named list
  result <- as.list(values)
  names(result) <- feature_names
  
  # Cleanup
  unlink(output_arff)
  
  if (verbose) {
    cat("Extracted", length(result), "features\n")
  }
  
  return(result)
}
