# superassp Code Examples & Patterns

## 1. TRK_* FUNCTION WRAPPER PATTERN

### Minimal Wrapper (Template)
```r
# From ssff_cpp_sptk_dio.R
trk_dio <- function(listOfFiles,
                    beginTime = 0.0,
                    endTime = 0.0,
                    windowShift = 10.0,
                    minF = 60.0,
                    maxF = 400.0,
                    voicing_threshold = 0.85,
                    toFile = TRUE,
                    explicitExt = "f0",
                    outputDirectory = NULL,
                    verbose = TRUE) {

  # Validation
  if (is.null(listOfFiles) || length(listOfFiles) == 0) {
    cli::cli_abort("No input files specified")
  }

  # Normalize paths
  listOfFiles <- normalizePath(path.expand(listOfFiles), mustWork = FALSE)
  files_exist <- file.exists(listOfFiles)
  if (!all(files_exist)) {
    cli::cli_abort(c("!", "Some files not found"))
  }

  # Setup
  n_files <- length(listOfFiles)
  beginTime <- if (length(beginTime) == 1) rep(beginTime, n_files) else beginTime
  endTime <- if (length(endTime) == 1) rep(endTime, n_files) else endTime

  makeOutputDirectory(outputDirectory, FALSE, "trk_dio")

  if (verbose) {
    cli::cli_inform("Processing {n_files} file{?s}")
  }

  # Process
  results <- vector("list", n_files)

  for (i in seq_len(n_files)) {
    tryCatch({
      # Load audio using av (universal format support)
      audio_obj <- av_to_asspDataObj(
        file_path = listOfFiles[i],
        start_time = beginTime[i],
        end_time = if (endTime[i] == 0.0) NULL else endTime[i]
      )

      # Call backend function
      dio_result <- dio_cpp(
        audio_obj = audio_obj,
        minF = minF,
        maxF = maxF,
        windowShift = windowShift,
        voicing_threshold = voicing_threshold,
        verbose = FALSE
      )

      # Create AsspDataObj
      out_obj <- create_f0_asspobj(dio_result, windowShift)

      if (toFile) {
        out_file <- generate_output_path(listOfFiles[i], explicitExt, outputDirectory)
        write.AsspDataObj(out_obj, out_file)
        results[[i]] <- TRUE
      } else {
        results[[i]] <- out_obj
      }

    }, error = function(e) {
      cli::cli_warn("Error in {basename(listOfFiles[i])}: {e$message}")
      results[[i]] <- if (toFile) FALSE else NULL
    })
  }

  # Return
  if (toFile) {
    n_success <- sum(unlist(results), na.rm = TRUE)
    if (verbose) cli::cli_inform("Processed {n_success}/{n_files} files")
    return(invisible(n_success))
  } else {
    results <- Filter(Negate(is.null), results)
    if (length(results) == 1) return(results[[1]]) else return(results)
  }
}

# Add function attributes
attr(trk_dio, "ext") <- "f0"
attr(trk_dio, "tracks") <- c("f0")
attr(trk_dio, "outputType") <- "SSFF"
attr(trk_dio, "nativeFiletypes") <- c("wav", "flac", "mp3", "mp4", "mkv", "avi")
attr(trk_dio, "suggestCaching") <- FALSE
```

### Python-Based Wrapper Pattern
```r
# From ssff_python_crepe.R
trk_crepe <- function(listOfFiles,
                      beginTime = 0,
                      endTime = 0,
                      windowShift = 5,
                      minF = 70,
                      maxF = 200,
                      voicing.threshold = 0.21,
                      model = c("tiny", "full"),
                      explicitExt = "crp",
                      outputDirectory = NULL,
                      toFile = TRUE) {

  # Validation
  if (length(listOfFiles) > 1 && !toFile) {
    stop("toFile=FALSE only allowed for single files")
  }

  # Load audio via av (in-memory, supports all formats)
  audio_data <- av::read_audio_bin(
    audio = origSoundFile,
    start_time = if (beginTime > 0) beginTime else NULL,
    end_time = if (endTime > 0) endTime else NULL,
    channels = 1
  )

  # Get sample rate
  fs <- attr(audio_data, "sample_rate")

  # Convert to numpy
  audio_float <- as.numeric(audio_data) / 2147483647.0
  np <- reticulate::import("numpy", convert = FALSE)
  audio_array <- np$array(audio_float, dtype = "float32")

  # Call Python
  py <- reticulate::import_main()
  py$waveform <- audio_array
  py$fs <- as.integer(fs)
  py$minF <- as.numeric(minF)
  py$maxF <- as.numeric(maxF)

  reticulate::py_run_string("
import torchcrepe
import numpy as np

# Extract pitch using CREPE
pitch, periodicity = torchcrepe.predict(
    waveform,
    sample_rate=fs,
    hop_length=int(fs * windowShift / 1000),
    fmin=minF,
    fmax=maxF,
    model='tiny',
    batch_size=2048,
    device='cpu',
    return_periodicity=True
)

# Convert to numpy
pitch = pitch.numpy()
periodicity = periodicity.numpy()
")

  # Extract results
  pitch_results <- py$pitch
  periodicity_results <- py$periodicity

  # Create AsspDataObj with equal-interval tracks
  inTable <- data.frame(
    f0 = pitch_results,
    periodicity = periodicity_results
  )

  # Set up AsspDataObj
  sampleRate <- 1000 / windowShift  # Frames per second
  outDataObj <- list()
  attr(outDataObj, "trackFormats") <- c("REAL32", "REAL32")
  attr(outDataObj, "sampleRate") <- sampleRate
  attr(outDataObj, "origFreq") <- as.numeric(fs)
  attr(outDataObj, "startTime") <- 0.0
  attr(outDataObj, "startRecord") <- 1
  attr(outDataObj, "endRecord") <- nrow(inTable)
  class(outDataObj) <- "AsspDataObj"

  AsspFileFormat(outDataObj) <- "SSFF"
  AsspDataFormat(outDataObj) <- as.integer(2)

  # Add tracks
  outDataObj <- addTrack(outDataObj, "fo[Hz]", 
                        as.matrix(inTable$f0), "REAL32")
  outDataObj <- addTrack(outDataObj, "periodicity",
                        as.matrix(inTable$periodicity), "REAL32")

  # Write to file
  if (toFile) {
    ssff_file <- sub("wav$", explicitExt, origSoundFile)
    attr(outDataObj, "filePath") <- ssff_file
    write.AsspDataObj(dobj = outDataObj, file = ssff_file)
  }

  return(outDataObj)
}
```

---

## 2. ASPDATAOBJ CREATION PATTERNS

### Minimal AsspDataObj
```r
# Single track, single column
create_minimal_assp <- function(data_vector, 
                                track_name = "data",
                                sample_rate = 100,
                                orig_freq = 16000,
                                format = "REAL32") {
  
  n_frames <- length(data_vector)
  
  obj <- list()
  attr(obj, "trackFormats") <- format
  attr(obj, "sampleRate") <- sample_rate
  attr(obj, "origFreq") <- orig_freq
  attr(obj, "startTime") <- 0.0
  attr(obj, "startRecord") <- 1
  attr(obj, "endRecord") <- n_frames
  attr(obj, "filePath") <- NA_character_
  
  class(obj) <- "AsspDataObj"
  AsspFileFormat(obj) <- "SSFF"
  AsspDataFormat(obj) <- as.integer(2)
  
  # Add track
  obj <- addTrack(obj, track_name, as.matrix(data_vector), format)
  
  return(obj)
}
```

### Multi-Column AsspDataObj (Formants)
```r
# Multiple formant tracks
create_formant_assp <- function(f1_hz, f2_hz, f3_hz,
                                b1_hz, b2_hz, b3_hz,
                                sample_rate = 100,
                                orig_freq = 16000) {
  
  n_frames <- length(f1_hz)
  
  obj <- list()
  attr(obj, "trackFormats") <- rep("REAL32", 6)
  attr(obj, "sampleRate") <- sample_rate
  attr(obj, "origFreq") <- orig_freq
  attr(obj, "startTime") <- 0.0
  attr(obj, "startRecord") <- 1
  attr(obj, "endRecord") <- n_frames
  
  class(obj) <- "AsspDataObj"
  AsspFileFormat(obj) <- "SSFF"
  AsspDataFormat(obj) <- as.integer(2)
  
  # Add tracks with template names
  obj <- addTrack(obj, "F1[Hz]", as.matrix(f1_hz), "REAL32")
  obj <- addTrack(obj, "F2[Hz]", as.matrix(f2_hz), "REAL32")
  obj <- addTrack(obj, "F3[Hz]", as.matrix(f3_hz), "REAL32")
  obj <- addTrack(obj, "B1[Hz]", as.matrix(b1_hz), "REAL32")
  obj <- addTrack(obj, "B2[Hz]", as.matrix(b2_hz), "REAL32")
  obj <- addTrack(obj, "B3[Hz]", as.matrix(b3_hz), "REAL32")
  
  return(obj)
}
```

### Multi-Column Matrix Track
```r
# Single track with multiple columns
create_multicolumn_assp <- function(data_matrix,
                                    track_name = "Fi[Hz]",
                                    sample_rate = 100,
                                    orig_freq = 16000) {
  
  n_frames <- nrow(data_matrix)
  n_cols <- ncol(data_matrix)
  
  obj <- list()
  attr(obj, "trackFormats") <- "REAL32"
  attr(obj, "sampleRate") <- sample_rate
  attr(obj, "origFreq") <- orig_freq
  attr(obj, "startTime") <- 0.0
  attr(obj, "startRecord") <- 1
  attr(obj, "endRecord") <- n_frames
  
  class(obj) <- "AsspDataObj"
  AsspFileFormat(obj) <- "SSFF"
  AsspDataFormat(obj) <- as.integer(2)
  
  # Add multi-column track
  obj <- addTrack(obj, track_name, data_matrix, "REAL32")
  
  return(obj)
}
```

---

## 3. EQUAL-INTERVAL PROCESSING IN PYTHON

### Basic Pattern
```python
import numpy as np

def analyze_with_equal_intervals(audio_array, sample_rate, frame_shift_ms=10.0):
    """
    Perform analysis returning results at equal time intervals.
    
    Args:
        audio_array: Input audio (numpy array, float64)
        sample_rate: Sample rate in Hz
        frame_shift_ms: Frame shift in milliseconds
    
    Returns:
        dict with equal-interval results
    """
    
    # Perform cycle-based analysis (may have irregular timing)
    cycle_times, f0_values = detect_f0_cycles(audio_array, sample_rate)
    
    # Create equal-interval grid
    duration = len(audio_array) / sample_rate
    frame_shift_sec = frame_shift_ms / 1000.0
    n_frames = int(np.ceil(duration / frame_shift_sec))
    frame_times = np.arange(n_frames) * frame_shift_sec
    
    # Interpolate to equal intervals
    f0_interp = np.zeros(n_frames, dtype=np.float32)
    voicing = np.zeros(n_frames, dtype=np.int16)
    
    for i, frame_time in enumerate(frame_times):
        # Find nearest analysis point
        half_shift = frame_shift_sec / 2.0
        nearby = np.where(np.abs(cycle_times - frame_time) <= half_shift)[0]
        
        if len(nearby) > 0:
            # Use nearest
            nearest_idx = nearby[np.argmin(np.abs(cycle_times[nearby] - frame_time))]
            f0_interp[i] = f0_values[nearest_idx]
            voicing[i] = 1
        else:
            f0_interp[i] = 0.0
            voicing[i] = 0
    
    return {
        'f0': f0_interp,
        'voicing': voicing,
        'times': frame_times,
        'raw_f0': f0_values,
        'raw_times': cycle_times,
        'sample_rate': sample_rate,
        'n_cycles': len(f0_values)
    }
```

### Linear Interpolation Pattern
```python
def analyze_with_linear_interpolation(audio_array, sample_rate, frame_shift_ms=10.0):
    """
    Interpolate analysis results to equal intervals using linear interpolation.
    """
    
    # Get cycle-based results
    cycle_times, cycle_values = detect_feature(audio_array, sample_rate)
    
    # Create grid
    duration = len(audio_array) / sample_rate
    frame_shift_sec = frame_shift_ms / 1000.0
    n_frames = int(np.ceil(duration / frame_shift_sec))
    frame_times = np.arange(n_frames) * frame_shift_sec
    
    # Linear interpolation
    if len(cycle_times) > 1:
        # Use scipy.interpolate.interp1d for smooth interpolation
        from scipy.interpolate import interp1d
        
        # Ensure monotonic times
        sort_idx = np.argsort(cycle_times)
        times_sorted = cycle_times[sort_idx]
        values_sorted = cycle_values[sort_idx]
        
        # Create interpolator
        f = interp1d(times_sorted, values_sorted, 
                    kind='linear', 
                    bounds_error=False,
                    fill_value=0.0)
        
        # Evaluate at grid points
        interp_values = f(frame_times)
    else:
        interp_values = np.zeros(n_frames, dtype=np.float32)
    
    return {
        'values': interp_values,
        'times': frame_times,
        'sample_rate': sample_rate
    }
```

---

## 4. AUDIO LOADING & CONVERSION

### Load Audio with av Package
```r
# Load audio from any format (in-memory)
load_audio_universal <- function(file_path, 
                                 start_time = 0,
                                 end_time = NULL,
                                 channels = 1) {
  
  # Get media info
  info <- av::av_media_info(file_path)
  
  if (length(info$audio) == 0) {
    stop("No audio stream found")
  }
  
  # Read audio
  audio_data <- av::read_audio_bin(
    audio = file_path,
    start_time = if (start_time > 0) start_time else NULL,
    end_time = end_time,
    channels = channels
  )
  
  # Extract properties
  sample_rate <- info$audio$sample_rate
  n_channels <- info$audio$channels
  
  # Normalize to float64 [-1, 1]
  audio_float <- as.numeric(audio_data) / 2147483648.0
  
  # De-interleave if needed
  if (n_channels > 1 && channels == 1) {
    n_frames <- length(audio_float) / n_channels
    audio_float <- audio_float[seq(1, length(audio_float), by = n_channels)]
  }
  
  return(list(
    samples = audio_float,
    sample_rate = sample_rate,
    channels = channels,
    duration = length(audio_float) / sample_rate
  ))
}
```

### Convert av Audio to Python NumPy
```r
av_to_numpy <- function(file_path,
                        start_time = 0,
                        end_time = NULL,
                        channels = 1) {
  
  # Load audio
  audio_data <- av::read_audio_bin(
    audio = file_path,
    start_time = if (start_time > 0) start_time else NULL,
    end_time = end_time,
    channels = channels
  )
  
  sample_rate <- attr(audio_data, "sample_rate")
  
  # Normalize to float
  audio_float <- as.numeric(audio_data) / 2147483648.0
  
  # Convert to numpy array
  np <- reticulate::import("numpy", convert = FALSE)
  audio_np <- np$array(audio_float, dtype = np$float64)
  
  return(list(
    array = audio_np,
    sample_rate = sample_rate
  ))
}
```

---

## 5. TRACK NAMING & DATA FRAME CONVERSION

### Expand Track Names
```r
# From track_helpers.R
expand_track_names <- function(track_names, n_cols) {
  """
  Expand template names like 'Fi[Hz]' to 'F1[Hz]', 'F2[Hz]', etc.
  """
  
  expanded <- character()
  
  for (name in track_names) {
    if (grepl("[A-Z]i(\\[|$)", name)) {
      # Has placeholder 'i'
      for (i in 1:n_cols) {
        expanded_name <- sub("([A-Z])i(\\[|$)", paste0("\\1", i, "\\2"), name)
        expanded <- c(expanded, expanded_name)
      }
    } else {
      # No placeholder
      expanded <- c(expanded, name)
    }
  }
  
  return(expanded)
}

# Usage
expand_track_names(c("Fi[Hz]", "Bi[Hz]", "RMS[dB]"), 3)
# Result: c("F1[Hz]", "F2[Hz]", "F3[Hz]", "B1[Hz]", "B2[Hz]", "B3[Hz]", "RMS[dB]")
```

### Convert AsspDataObj to Data Frame
```r
assp_to_dataframe <- function(assp_obj, clean_names = TRUE) {
  
  # Get sampling info
  sample_rate <- attr(assp_obj, "sampleRate")
  start_time <- attr(assp_obj, "startTime")
  n_frames <- attr(assp_obj, "endRecord") - attr(assp_obj, "startRecord") + 1
  
  # Create time column
  frame_times <- start_time + (0:(n_frames - 1)) / sample_rate
  
  df <- data.frame(time = frame_times)
  
  # Get track names
  track_names <- names(assp_obj)
  
  # Add tracks
  for (i in seq_along(track_names)) {
    track_data <- assp_obj[[track_names[i]]]
    
    if (is.matrix(track_data) && ncol(track_data) > 1) {
      # Multi-column track
      col_names <- .expand_track_template(track_names[i], ncol(track_data))
      
      for (j in seq_len(ncol(track_data))) {
        col_name <- col_names[j]
        if (clean_names) {
          col_name <- gsub("[\\[\\]]", "_", col_name)
          col_name <- gsub("-", "_", col_name)
        }
        df[[col_name]] <- track_data[, j]
      }
    } else {
      # Single-column track
      col_name <- track_names[i]
      if (clean_names) {
        col_name <- gsub("[\\[\\]]", "_", col_name)
        col_name <- gsub("-", "_", col_name)
      }
      df[[col_name]] <- as.vector(track_data)
    }
  }
  
  return(df)
}
```

---

## 6. WORKING WITH MULTIPLE FILES

### Batch Processing
```r
process_batch <- function(file_list,
                         trk_func,
                         output_dir = NULL,
                         verbose = TRUE,
                         ...) {
  
  n_files <- length(file_list)
  results <- vector("list", n_files)
  success_count <- 0
  
  if (verbose && n_files > 1) {
    cli::cli_progress_bar("Processing", total = n_files)
  }
  
  for (i in seq_len(n_files)) {
    tryCatch({
      result <- trk_func(
        file_list[i],
        outputDirectory = output_dir,
        verbose = FALSE,
        ...
      )
      
      if (is.null(result)) {
        results[[i]] <- FALSE
      } else {
        results[[i]] <- TRUE
        success_count <- success_count + 1
      }
      
    }, error = function(e) {
      if (verbose) {
        cli::cli_warn("Error in {basename(file_list[i])}: {e$message}")
      }
      results[[i]] <- FALSE
    })
    
    if (verbose && n_files > 1) {
      cli::cli_progress_update()
    }
  }
  
  if (verbose) {
    cli::cli_inform("Successfully processed {success_count}/{n_files} files")
  }
  
  return(invisible(success_count))
}

# Usage
success <- process_batch(
  file_list = c("file1.wav", "file2.wav", "file3.wav"),
  trk_func = trk_dio,
  output_dir = "output/",
  minF = 70,
  maxF = 300
)
```

---

## 7. VALIDATION & ERROR HANDLING

### Validate AsspDataObj
```r
validate_assp <- function(obj) {
  
  # Check class
  if (!inherits(obj, "AsspDataObj")) {
    stop("Object is not an AsspDataObj")
  }
  
  # Check required attributes
  required_attrs <- c("trackFormats", "sampleRate", "startTime", "startRecord", "endRecord")
  missing <- setdiff(required_attrs, names(attributes(obj)))
  if (length(missing) > 0) {
    stop(paste("Missing attributes:", paste(missing, collapse = ", ")))
  }
  
  # Check consistency
  n_tracks <- length(obj)
  n_formats <- length(attr(obj, "trackFormats"))
  if (n_tracks != n_formats) {
    stop(paste("Number of tracks", n_tracks, "doesn't match formats", n_formats))
  }
  
  # Check frame counts
  n_frames_expected <- attr(obj, "endRecord") - attr(obj, "startRecord") + 1
  for (i in seq_len(n_tracks)) {
    track_data <- obj[[i]]
    if (is.matrix(track_data)) {
      n_frames <- nrow(track_data)
    } else {
      n_frames <- length(track_data)
    }
    if (n_frames != n_frames_expected) {
      stop(paste("Track", i, "has", n_frames, "frames, expected", n_frames_expected))
    }
  }
  
  return(TRUE)
}
```

### Safe File Processing
```r
safe_process_file <- function(file_path, process_func) {
  
  tryCatch({
    # Check file exists
    if (!file.exists(file_path)) {
      warning("File not found: ", file_path)
      return(NULL)
    }
    
    # Check file is readable
    if (!file.access(file_path, mode = 4) == 0) {
      warning("Cannot read file: ", file_path)
      return(NULL)
    }
    
    # Process
    result <- process_func(file_path)
    
    # Validate result
    if (!is.null(result) && !is.AsspDataObj(result)) {
      warning("Invalid result for: ", file_path)
      return(NULL)
    }
    
    return(result)
    
  }, error = function(e) {
    warning("Error processing ", file_path, ": ", e$message)
    return(NULL)
  }, warning = function(w) {
    warning("Warning processing ", file_path, ": ", w$message)
    return(NULL)
  })
}
```

---

## 8. PERFORMANCE TIPS

### Vectorized Operations
```r
# Instead of loop
bad_way <- function(files) {
  results <- list()
  for (i in seq_along(files)) {
    audio <- av::read_audio_bin(files[i], channels = 1)
    results[[i]] <- list(samples = audio, sr = attr(audio, "sample_rate"))
  }
  return(results)
}

# Better: let trk_* handle batch
good_way <- function(files) {
  # trk_* functions already handle batches efficiently
  result <- trk_dio(files, toFile = TRUE)
  return(result)
}
```

### Parallel Processing
```r
# Use pbmcapply for progress + parallel
if (requireNamespace("pbmcapply", quietly = TRUE)) {
  results <- pbmcapply::pbmclapply(
    file_list,
    function(file) {
      # Process single file
      trk_dio(file, outputDirectory = "output/", verbose = FALSE)
    },
    mc.cores = parallel::detectCores()
  )
}
```

### Memory Efficiency
```r
# Process large files in chunks
process_large_file <- function(file_path, chunk_duration = 10) {
  
  info <- av::av_media_info(file_path)
  total_duration <- info$duration
  
  n_chunks <- ceiling(total_duration / chunk_duration)
  results <- list()
  
  for (i in 1:n_chunks) {
    start_time <- (i - 1) * chunk_duration
    end_time <- min(i * chunk_duration, total_duration)
    
    # Process chunk
    chunk_result <- trk_dio(
      file_path,
      beginTime = start_time,
      endTime = end_time,
      toFile = TRUE,
      outputDirectory = "chunks/"
    )
    
    results[[i]] <- chunk_result
  }
  
  return(results)
}
```

