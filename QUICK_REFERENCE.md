# Quick Reference: Python DSP Function Migration

## Template Usage

Copy this template to create new optimized Python DSP functions:

```r
#' Your function description
#' @inheritParams swipe_opt
#' @param custom_param Your custom parameter
#' @export
your_function_opt <- function(listOfFiles = NULL,
                             beginTime = 0.0,
                             endTime = 0.0,
                             windowShift = 5.0,
                             custom_param = default_value,
                             explicitExt = "ext",
                             outputDirectory = NULL,
                             toFile = TRUE,
                             assertLossless = NULL,
                             logToFile = FALSE,
                             keepConverted = FALSE,
                             verbose = TRUE) {
  
  ## Setup
  explicitExt <- ifelse(is.null(explicitExt), "ext", explicitExt)
  newTracknames <- c("track1", "track2")  # Customize
  nativeFiletypes <- c("wav")
  
  currCall <- rlang::current_call()
  funName <- rlang::call_name(currCall)
  preferedFiletype <- nativeFiletypes[[1]]
  knownLossless <- c(assertLossless, knownLossless())
  
  beginTime <- if(is.null(beginTime)) 0.0 else beginTime
  endTime <- if(is.null(endTime)) 0.0 else endTime
  n_files <- length(listOfFiles)
  
  # Validate
  if(length(beginTime) > 1 && length(beginTime) != n_files) {
    cli::cli_abort("beginTime must be length 1 or match listOfFiles length")
  }
  if(length(endTime) > 1 && length(endTime) != n_files) {
    cli::cli_abort("endTime must be length 1 or match listOfFiles length")
  }
  
  beginTime <- fast_recycle_times(beginTime, n_files)
  endTime <- fast_recycle_times(endTime, n_files)
  
  makeOutputDirectory(outputDirectory, logToFile, funName)
  
  ## File preparation
  listOfFiles <- fast_strip_file_protocol(listOfFiles)
  listOfFiles <- normalizePath(path.expand(listOfFiles))
  
  file_exts <- fast_file_ext(listOfFiles)
  is_native <- fast_is_native(file_exts, nativeFiletypes)
  needs_timewindow <- (beginTime != 0.0 | endTime != 0.0)
  
  if(all(is_native) && !any(needs_timewindow)) {
    if(verbose) {
      cli::cli_inform("Applying {.fun {funName}} to {cli::no(n_files)} recording{?s}")
    }
    listOfFilesDF <- data.frame(
      audio = listOfFiles, dsp_input = listOfFiles,
      beginTime = beginTime, endTime = endTime, stringsAsFactors = FALSE
    )
    toClear <- character(0)
  } else {
    listOfFiles_toClear <- convertInputMediaFiles(
      listOfFiles, beginTime, endTime, windowShift,
      nativeFiletypes, preferedFiletype, knownLossless,
      funName, keepConverted, verbose
    )
    listOfFilesDF <- listOfFiles_toClear[[1]]
    toClear <- listOfFiles_toClear[[3]]
    
    file_exts <- fast_file_ext(listOfFilesDF$dsp_input)
    if(!all(fast_is_native(file_exts, nativeFiletypes))) {
      cli::cli_abort("File conversion failed")
    }
  }
  
  ## Processing
  if(verbose) {
    cli::cli_inform("Applying {.fun {funName}} to {cli::no(n_files)} recording{?s}")
  }
  
  # Check Python module
  if(!reticulate::py_module_available("your_python_module")) {
    cli::cli_abort(c(
      "Python module {.pkg your_python_module} is not available.",
      "i" = "Install with: pip install your_python_module"
    ))
  }
  
  externalRes <- vector("list", nrow(listOfFilesDF))
  
  for(i in seq_len(nrow(listOfFilesDF))) {
    origSoundFile <- listOfFilesDF$dsp_input[i]
    bt <- listOfFilesDF$beginTime[i]
    et <- listOfFilesDF$endTime[i]
    
    tryCatch({
      externalRes[[i]] <- process_your_function_single(
        origSoundFile, bt, et, windowShift, custom_param,
        explicitExt, outputDirectory, toFile
      )
    }, error = function(e) {
      cli::cli_abort(c(
        "Error processing {.file {basename(origSoundFile)}}",
        "x" = conditionMessage(e)
      ))
    })
  }
  
  if(!toFile && !is.null(newTracknames)) {
    n_tracks <- length(names(externalRes[[1]]))
    if(n_tracks != length(newTracknames)) {
      cli::cli_abort("Wrong number of track names")
    }
    externalRes <- fast_rename_tracks(externalRes, newTracknames)
  }
  
  if(n_files == 1) externalRes <- externalRes[[1]]
  
  ## Cleanup
  cleanupConvertedInputMediaFiles(toClear, keepConverted, verbose)
  
  if(toFile) return(n_files) else return(externalRes)
}

#' Process single file
#' @keywords internal
process_your_function_single <- function(soundFile, beginTime, endTime,
                                        windowShift, custom_param,
                                        explicitExt, outputDirectory, toFile) {
  
  # Initialize Python (once per session)
  if(!exists(".py_your_function_initialized", envir = .GlobalEnv)) {
    reticulate::py_run_string("
import numpy as np
import gc
import librosa as lr
import your_python_module
")
    assign(".py_your_function_initialized", TRUE, envir = .GlobalEnv)
  }
  
  # Set parameters
  py <- reticulate::import_main()
  py$soundFile <- soundFile
  py$ws <- windowShift
  py$custom <- custom_param
  py$bt <- beginTime
  py$et <- endTime
  
  # Run Python DSP code
  reticulate::py_run_string("
if et > 0:
    x, fs = lr.load(soundFile, dtype=np.float64, offset=bt, duration=et - bt)
else:
    x, fs = lr.load(soundFile, dtype=np.float64, offset=bt)

# Your DSP computation here
result1 = your_python_module.compute1(x, fs, ...)
result2 = your_python_module.compute2(x, fs, ...)

del x
gc.collect()
")
  
  # Build result table
  inTable <- data.frame(
    track1 = py$result1,
    track2 = py$result2
  )
  
  sampleRate <- 1 / windowShift * 1000
  startTime <- 1 / sampleRate  # Or from Python: py$times[[1]]
  
  # Create AsspDataObj
  outDataObj <- list()
  attr(outDataObj, "trackFormats") <- c("INT16", "INT16")  # Adjust as needed
  attr(outDataObj, "sampleRate") <- sampleRate
  attr(outDataObj, "origFreq") <- as.numeric(py$fs)
  attr(outDataObj, "startTime") <- as.numeric(startTime)
  attr(outDataObj, "startRecord") <- as.integer(1)
  attr(outDataObj, "endRecord") <- as.integer(nrow(inTable))
  class(outDataObj) <- "AsspDataObj"
  
  AsspFileFormat(outDataObj) <- "SSFF"
  AsspDataFormat(outDataObj) <- as.integer(2)
  
  # Add tracks
  for(track_name in names(inTable)) {
    trackTable <- inTable %>%
      dplyr::select(!!track_name) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(dplyr::across(tidyselect::everything(), as.integer))
    
    names(trackTable) <- NULL
    outDataObj <- addTrack(outDataObj, track_name, 
                          as.matrix(trackTable[,1]), "INT16")
  }
  
  # Fix missing values at start
  if(startTime > (1/sampleRate)) {
    nr_of_missing_samples <- as.integer(floor(startTime / (1/sampleRate)))
    
    for(track_name in names(inTable)) {
      missing_vals <- matrix(0, nrow = nr_of_missing_samples, 
                            ncol = ncol(outDataObj[[track_name]]))
      outDataObj[[track_name]] <- rbind(missing_vals, outDataObj[[track_name]])
    }
    
    attr(outDataObj, "startTime") <- startTime - nr_of_missing_samples * (1/sampleRate)
  }
  
  assertthat::assert_that(is.AsspDataObj(outDataObj),
                          msg = "Invalid AsspDataObj created")
  
  # Write to file if requested
  if(toFile) {
    ssff_file <- sub("wav$", explicitExt, soundFile)
    if(!is.null(outputDirectory)) {
      ssff_file <- file.path(outputDirectory, basename(ssff_file))
    }
    attr(outDataObj, "filePath") <- as.character(ssff_file)
    write.AsspDataObj(dobj = outDataObj, file = ssff_file)
  }
  
  return(outDataObj)
}

# Set function attributes
attr(your_function_opt, "ext") <- "ext"
attr(your_function_opt, "tracks") <- c("track1", "track2")
attr(your_function_opt, "outputType") <- "SSFF"
attr(your_function_opt, "nativeFiletypes") <- c("wav")
attr(your_function_opt, "suggestCaching") <- FALSE
```

## Checklist

When implementing a new function:

- [ ] Copy template above
- [ ] Replace `your_function` with actual name
- [ ] Update parameters and documentation
- [ ] Set correct `explicitExt` and `newTracknames`
- [ ] Update Python module check
- [ ] Implement Python DSP code in string
- [ ] Set correct track formats (INT16, REAL32, etc.)
- [ ] Update function attributes at end
- [ ] Test with single file
- [ ] Test with batch processing
- [ ] Run benchmark comparison
- [ ] Update documentation

## Common Track Formats

```r
# Integer (16-bit)
attr(outDataObj, "trackFormats") <- c("INT16")

# Float (32-bit)
attr(outDataObj, "trackFormats") <- c("REAL32")

# Multiple tracks
attr(outDataObj, "trackFormats") <- c("INT16", "REAL32", "INT16")
```

## Common Python Patterns

### Loading audio with time windowing
```python
if et > 0:
    x, fs = lr.load(soundFile, dtype=np.float64, offset=bt, duration=et - bt)
else:
    x, fs = lr.load(soundFile, dtype=np.float64, offset=bt)
```

### Memory cleanup
```python
del x, intermediate_result
gc.collect()
```

### Convert to specific format
```python
# For REAPER (needs int16)
raw_x = x * 2**15
int_x = raw_x.astype(np.int16)
```

## Testing Commands

```r
# Source your implementation
source("R/python_ssff_optimized.R")

# Test single file
result <- your_function_opt(
  "tests/signalfiles/AVQI/input/sv1.wav",
  toFile = FALSE,
  verbose = TRUE
)

# Check structure
str(result)
names(result)

# Test batch
files <- list.files("tests/signalfiles/AVQI/input/", 
                   pattern = "*.wav", full.names = TRUE)
count <- your_function_opt(files, toFile = TRUE)

# Benchmark
system.time({
  your_function_opt(files, toFile = FALSE, verbose = FALSE)
})
```

## Common Issues

### Python module not found
```r
reticulate::py_install("module_name")
# or
system("pip install module_name")
```

### Wrong number of samples
Check that `windowShift` is correctly converted (ms vs seconds).

### Track name mismatch
Ensure `newTracknames` length matches actual tracks added.

### Memory issues
Add explicit garbage collection:
```python
del x
gc.collect()
```

### File not found
Check that `fast_strip_file_protocol()` is called.

## Performance Tips

1. **Initialize Python once**: Use `.GlobalEnv` flag
2. **Batch file operations**: Use `convertInputMediaFiles()`
3. **Avoid repeated normalizePath**: Store results
4. **Use Rcpp helpers**: Much faster than R equivalents
5. **Preallocate lists**: `vector("list", n)` not `list()`
