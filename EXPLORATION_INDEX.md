# superassp Package Exploration - Complete Index

**Exploration Date:** October 27, 2025
**Package Version:** 0.8.4
**Author:** Fredrik Nylén (Umeå University)

This exploration provides comprehensive documentation of the superassp R package structure, particularly focusing on understanding trk_* functions, AsspDataObj architecture, and DSP integration patterns.

---

## Documentation Files Generated

### 1. SUPERASSP_EXPLORATION_SUMMARY.md (Primary Document)
**Comprehensive structural analysis** of the entire package.

**Contents:**
- Complete package structure and file organization (96 R files, 27,843 lines)
- DESCRIPTION file analysis (key dependencies: reticulate, av, S7, Rcpp)
- NAMESPACE exports (160+ functions across 8 categories)
- AsspDataObj S3 class structure and attributes
- Equal-interval signal processing architecture
- 40+ trk_* functions organized by backend
- Python integration architecture (50+ modules in inst/python/)
- Track naming convention (bracket notation)
- Data structure conversion helpers
- 12 major sections with detailed tables and examples

**Key Sections:**
1. Package structure
2. DESCRIPTION dependencies
3. NAMESPACE & exports
4. AsspDataObj definition
5. trk_* function inventory
6. Python integration
7. Data conversion helpers
8. Track naming & labels
9. Signal processing patterns
10. Key files reference
11. Function metadata
12. Dependencies & ecosystem

**Use This For:** Understanding overall package architecture, finding which files contain which functions, learning the standard patterns used throughout.

---

### 2. SUPERASSP_CODE_EXAMPLES.md (Practical Guide)
**Ready-to-use code patterns and examples** for implementing new functionality.

**Contents:**
- trk_* function wrapper templates (C++ and Python versions)
- AsspDataObj creation patterns (minimal, multi-column, matrix)
- Equal-interval processing in Python (cycle-based → grid-based)
- Audio loading and format conversion
- Track naming and data frame conversion
- Batch processing examples
- Validation and error handling
- Performance optimization tips

**Key Patterns:**
1. Complete trk_dio() C++ wrapper template
2. Complete trk_crepe() Python wrapper template
3. Minimal AsspDataObj creation
4. Multi-column formant AsspDataObj
5. Python equal-interval grid creation
6. Audio loading with universal format support
7. Track name expansion patterns
8. Batch file processing
9. Error handling with tryCatch
10. Parallel processing setup
11. Memory-efficient chunk processing

**Use This For:** Implementing new trk_* functions, creating AsspDataObj correctly, converting between formats, understanding error handling patterns.

---

## Key Findings Summary

### Package Statistics
- **96 R files** with ~27,843 total lines
- **160+ exported functions**
- **60+ trk_* functions** (track processing)
- **50+ Python modules** in inst/python/
- **220+ man pages** (well documented)

### Architecture Highlights

#### 1. Core Data Structure: AsspDataObj
- S3 class (list-based with attributes)
- Stores **equal-interval** analysis tracks
- SSFF file format compatibility
- Tracks are matrices with metadata
- Template naming: `Fi[Hz]` → `F1[Hz]`, `F2[Hz]`, etc.

#### 2. Processing Pattern
```
Audio File (any format)
    ↓ av::read_audio_bin() [in-memory]
    ↓ normalize to float
    ↓ reticulate → Python or Rcpp backend
    ↓ equal-interval analysis
    ↓ create AsspDataObj
    ↓ write.AsspDataObj() → SSFF file
```

#### 3. Function Organization
- **ssff_*.R** files: 85% of trk_* implementations
- **Direct R files:** 3 (trk_creak_union, trk_egg_f0, trk_formants_tvwlp)
- **Naming:** `ssff_<backend>_<method>.R` → `trk_<method>()`

#### 4. Python Integration
- **reticulate** bridges R ↔ Python
- **av package** loads all audio formats in-memory
- **NumPy arrays** for data passing
- **50+ Python modules** for various analyses
- **Equal-interval conversion** in Python (see egg_f0.py pattern)

#### 5. Backend Support
- **C++/Rcpp:** SPTK pitch trackers (DIO, Harvest, RAPT, Swipe, Reaper)
- **C/ASSP:** Forest, RMS, Autocorr, Zero-crossing
- **Python:** Librosa, Torch, SciPy, custom modules
- **Praat/Parselmouth:** Formants, pitch, intensity
- **MATLAB Legacy:** Support via R.matlab

---

## File Location Reference

### Primary Files by Topic

**Core Infrastructure:**
- `/R/assp_dataobj.R` - AsspDataObj methods
- `/R/wrassp_AsspDataObj.R` - SSFF I/O
- `/R/s7_avaudio.R` - AVAudio S7 class
- `/R/track_helpers.R` - Name manipulation
- `/R/track_labels_plotmath.R` - Label generation

**Audio I/O:**
- `/R/av_helpers.R` - Universal audio loading
- `/R/av_python_helpers.R` - NumPy conversion
- `/R/prep_recode.R` - Format conversion

**trk_* Direct Implementations:**
- `/R/trk_creak_union.R` - Creaky voice (Union Method)
- `/R/trk_egg_f0.R` - EGG F0/Oq analysis
- `/R/trk_formants_tvwlp.R` - TVWLP formants

**trk_* Wrappers (Major):**
- `/R/ssff_cpp_sptk_dio.R`, `*_harvest.R`, `*_rapt.R`, etc. - SPTK pitch
- `/R/ssff_c_assp_*.R` - ASSP analysis
- `/R/ssff_python_crepe.R`, `*_yin.R`, `*_pyin.R` - Python pitch trackers
- `/R/ssff_python_pm_*.R` - Praat/Parselmouth methods

**Feature Extraction:**
- `/R/list_*.R` files (30+) - Feature list functions

**Python Modules:**
- `/inst/python/egg_analysis/` - EGG F0 analysis
- `/inst/python/ftrack_tvwlp/` - TVWLP formants
- `/inst/python/union-creak-detection-method/` - Creaky voice
- `/inst/python/covarep_python/` - COVAREP (50+ files)
- `/inst/python/voice_analysis_python/` - Voice analysis (80+ files)
- Other specialized modules for SNACK, Praat, etc.

---

## Standard Workflow for New Functions

### Step 1: Choose Pattern
**If C++/compiled:**
→ Use ssff_cpp_*.R wrapper pattern
→ Call backend function (e.g., `dio_cpp()`)

**If Python:**
→ Use ssff_python_*.R wrapper pattern
→ Use av for audio loading
→ Convert to NumPy array
→ Call Python via reticulate

**If Complex R logic:**
→ Use direct file like trk_egg_f0.R
→ Include helper functions (`.function_name()`)

### Step 2: Audio Loading
```r
audio_obj <- av_to_asspDataObj(file_path, start_time, end_time)
# OR for Python:
audio_data <- av::read_audio_bin(file_path, channels = 1)
audio_float <- as.numeric(audio_data) / 2147483648.0
```

### Step 3: Analysis
- Compute features at **equal time intervals**
- Create time grid: `np.arange(n_frames) * frame_shift_sec`
- Interpolate results to grid

### Step 4: Create AsspDataObj
```r
obj <- list()
attr(obj, "trackFormats") <- c("REAL32", "INT16")
attr(obj, "sampleRate") <- 100  # frames per second
attr(obj, "origFreq") <- 16000  # original audio rate
attr(obj, "startTime") <- 0.0
attr(obj, "startRecord") <- 1
attr(obj, "endRecord") <- n_frames
class(obj) <- "AsspDataObj"
AsspFileFormat(obj) <- "SSFF"
AsspDataFormat(obj) <- as.integer(2)

obj <- addTrack(obj, "fo[Hz]", as.matrix(f0_values), "REAL32")
obj <- addTrack(obj, "voicing", as.matrix(voicing), "INT16")
```

### Step 5: Add Function Metadata
```r
attr(trk_myfunction, "ext") <- "myx"
attr(trk_myfunction, "tracks") <- c("track1", "track2")
attr(trk_myfunction, "outputType") <- "SSFF"
attr(trk_myfunction, "nativeFiletypes") <- c("wav", "mp3")
attr(trk_myfunction, "suggestCaching") <- FALSE
```

### Step 6: Test and Document
- Validate AsspDataObj with `is.AsspDataObj()`
- Test batch processing
- Add roxygen2 documentation
- Create man page via `devtools::document()`

---

## Common Pitfalls & Solutions

### Problem: Wrong samplingRate
**Solution:** samplingRate = 1 / windowShift (in seconds)
- 10ms shift → samplingRate = 100 Hz
- 5ms shift → samplingRate = 200 Hz

### Problem: Unequal-interval data
**Solution:** Always interpolate cycle-based results to equal grid
- See egg_f0.py lines 153-193 for pattern

### Problem: Audio format not supported
**Solution:** Use av package, not direct file reading
- `av::read_audio_bin()` supports 50+ formats
- Automatic fallback to wrassp for legacy formats

### Problem: Multi-column tracks
**Solution:** Use template names with 'i' placeholder
- `Fi[Hz]` automatically expands to F1, F2, F3, etc.
- See `.expand_track_template()` in track_helpers.R

### Problem: Python module not found
**Solution:** Place in inst/python/ and use reticulate::import_from_path()
- Or use reticulate::py_install() for pip packages
- See trk_egg_f0.R lines 223-244 for pattern

---

## Quick Reference: Key Functions

### Data I/O
- `read.AsspDataObj(file)` - Read SSFF file
- `write.AsspDataObj(obj, file)` - Write SSFF file
- `av_to_asspDataObj(file)` - Load audio (universal format)
- `read_avaudio(file)` - Create AVAudio object
- `av_to_python_audio(data, sr)` - Convert to NumPy

### Data Manipulation
- `addTrack(obj, name, data, format)` - Add track
- `delTrack(obj, name)` - Remove track
- `as.data.frame(obj)` - Convert to data frame
- `as_tibble(obj)` - Convert to tibble
- `cut.AsspDataObj(obj, begin, end)` - Time-domain slice

### Validation
- `is.AsspDataObj(obj)` - Check validity
- `validate_assp(obj)` - Comprehensive check

### Track Naming
- `.has_placeholder(name)` - Check for 'i' placeholder
- `.expand_track_template(template, n)` - Expand Fi[Hz] → F1, F2, ...
- `.clean_track_names(names)` - fo[Hz] → fo_Hz
- `get_track_label(name)` - Get display label

### Plotting
- `ggtrack(data, x, y)` - Auto-labeled ggplot2

---

## Example: Building a New trk_* Function

**Goal:** Create `trk_myalgorithm()` that analyzes audio

**File:** `/R/ssff_python_myalgorithm.R`

```r
#' My Signal Processing Algorithm
#' 
#' Analyzes audio to extract features.
#'
#' @param listOfFiles Vector of audio file paths
#' @param beginTime Start time in seconds (default 0)
#' @param endTime End time in seconds (default 0)
#' @param windowShift Frame shift in milliseconds (default 10)
#' @param param1 Algorithm parameter (default 0.5)
#' @param toFile Write to SSFF file? (default TRUE)
#' @param outputDirectory Output directory (default NULL = same as input)
#' @param verbose Print progress? (default TRUE)
#'
#' @return If toFile=TRUE: number of processed files
#'         If toFile=FALSE: AsspDataObj or list of AsspDataObj
#'
#' @export
trk_myalgorithm <- function(listOfFiles,
                            beginTime = 0.0,
                            endTime = 0.0,
                            windowShift = 10.0,
                            param1 = 0.5,
                            toFile = TRUE,
                            outputDirectory = NULL,
                            verbose = TRUE) {
  
  # Validation
  if (is.null(listOfFiles) || length(listOfFiles) == 0) {
    cli::cli_abort("No files specified")
  }
  
  listOfFiles <- normalizePath(listOfFiles, mustWork = FALSE)
  if (!all(file.exists(listOfFiles))) {
    cli::cli_abort("Some files not found")
  }
  
  # Setup
  n_files <- length(listOfFiles)
  results <- vector("list", n_files)
  
  if (verbose) {
    cli::cli_inform("Processing {n_files} file{?s}")
  }
  
  # Main loop
  for (i in seq_len(n_files)) {
    tryCatch({
      # Load audio
      audio_obj <- av_to_asspDataObj(
        listOfFiles[i],
        start_time = beginTime,
        end_time = if (endTime == 0.0) NULL else endTime
      )
      
      # Get audio data for processing
      audio_float <- as.numeric(audio_obj$signal) / 2147483648.0
      fs <- attr(audio_obj, "origFreq")
      
      # Call analysis (Python backend example)
      py_result <- my_algorithm_python(
        audio = audio_float,
        sample_rate = fs,
        window_shift_ms = windowShift,
        param1 = param1
      )
      
      # Create AsspDataObj
      obj <- list()
      attr(obj, "trackFormats") <- "REAL32"
      attr(obj, "sampleRate") <- 1000 / windowShift
      attr(obj, "origFreq") <- fs
      attr(obj, "startTime") <- 0.0
      attr(obj, "startRecord") <- 1
      attr(obj, "endRecord") <- length(py_result)
      class(obj) <- "AsspDataObj"
      AsspFileFormat(obj) <- "SSFF"
      AsspDataFormat(obj) <- as.integer(2)
      
      obj <- addTrack(obj, "result", as.matrix(py_result), "REAL32")
      
      if (toFile) {
        out_file <- sub("\\.[^.]+$", ".myx", listOfFiles[i])
        if (!is.null(outputDirectory)) {
          out_file <- file.path(outputDirectory, basename(out_file))
        }
        write.AsspDataObj(obj, out_file)
        results[[i]] <- TRUE
      } else {
        results[[i]] <- obj
      }
      
    }, error = function(e) {
      cli::cli_warn("Error: {e$message}")
      results[[i]] <- if (toFile) FALSE else NULL
    })
  }
  
  # Return
  if (toFile) {
    n_success <- sum(unlist(results))
    if (verbose) cli::cli_inform("Processed {n_success}/{n_files} files")
    return(invisible(n_success))
  } else {
    results <- Filter(Negate(is.null), results)
    if (length(results) == 1) results[[1]] else results
  }
}

# Add metadata
attr(trk_myalgorithm, "ext") <- "myx"
attr(trk_myalgorithm, "tracks") <- c("result")
attr(trk_myalgorithm, "outputType") <- "SSFF"
attr(trk_myalgorithm, "nativeFiletypes") <- c("wav", "mp3")
attr(trk_myalgorithm, "suggestCaching") <- FALSE
```

---

## Further Reading

### Essential Files to Study
1. **assp_dataobj.R** - Understand AsspDataObj handling
2. **track_helpers.R** - Learn naming patterns
3. **ssff_cpp_sptk_dio.R** - Study C++ wrapper pattern
4. **ssff_python_crepe.R** - Study Python wrapper pattern
5. **trk_egg_f0.R** - Study complex R implementation
6. **egg_analysis/egg_f0.py** - Study Python equal-interval pattern

### Dependencies to Understand
- **wrassp** - Original SSFF format library
- **reticulate** - R ↔ Python bridge
- **av** - ffmpeg wrapper for audio
- **Rcpp** - R ↔ C++ integration
- **S7** - Modern R OOP system

---

## Document Metadata

**Files in this Exploration:**
1. `SUPERASSP_EXPLORATION_SUMMARY.md` - Comprehensive structure (primary)
2. `SUPERASSP_CODE_EXAMPLES.md` - Practical examples and patterns
3. `EXPLORATION_INDEX.md` - This file (guide and index)

**Total Content:** ~10,000 lines of documentation

**Last Updated:** October 27, 2025
**Superassp Version Documented:** 0.8.4
**Git Status:** Clean (master branch)

