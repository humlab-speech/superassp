# Code Quality Analysis - JSTF Phase 2 Integration

## Analysis Date
2025-11-10

## Scope
Analysis of 6 functions integrated with JSTF file output support in Phase 2 Week 3:
- `lst_covarep_vq()` - R/covarep_vq.R
- `lst_avqip()` - R/list_python_pm_pavqi.R
- `lst_dsip()` - R/list_python_pm_pdsi.R
- `lst_voice_reportp()` - R/list_python_pm_pvoice_report.R
- `lst_voice_tremorp()` - R/list_python_pm_pvoice_tremor.R
- `lst_phonet()` - R/ssff_python_phonet.R

## Code Quality Assessment

### ✅ Strengths

1. **Consistent Pattern Application**
   - All 6 functions follow identical 5-step integration pattern
   - Predictable structure makes code easy to maintain
   - Parameter naming is standardized across functions

2. **Backward Compatibility**
   - `toFile = FALSE` default preserves existing behavior
   - No breaking changes to existing API
   - Users can opt-in to JSTF output

3. **Documentation Quality**
   - Comprehensive roxygen2 documentation
   - Usage examples show both modes (in-memory and file output)
   - Clear parameter descriptions

4. **Metadata Capture**
   - All DSP parameters captured in JSTF metadata
   - Audio metadata (sample rate, duration) included
   - Time ranges properly recorded

5. **Error Handling**
   - NULL results handled gracefully (NA_character_ in output paths)
   - File existence validation
   - Input parameter validation

### ⚠️ Areas for Improvement

#### 1. Code Duplication (DRY Violation)

**Issue**: JSTF file writing logic is duplicated across all 6 functions (~40-60 lines each)

**Current Pattern** (repeated 6 times):
```r
if (toFile) {
  output_paths <- character(n_files)
  for (i in seq_len(n_files)) {
    result <- results[[i]]
    if (!is.null(result)) {
      audio_info <- av::av_media_info(file_path)
      json_obj <- create_json_track_obj(...)
      write_json_track(json_obj, output_path)
      output_paths[i] <- output_path
    }
  }
  return(invisible(output_paths))
}
```

**Impact**:
- ~240-360 lines of duplicated code
- Changes require updating 6+ locations
- Higher maintenance burden

**Recommendation**: Extract to helper function

#### 2. Magic Numbers

**Issue**: Hardcoded values without named constants

Examples:
- `16000` (Phonet sample rate requirement)
- Time conversions (ms to seconds: `/1000`)
- Default thresholds scattered across functions

**Recommendation**: Define package-level constants

#### 3. Mixed Concerns

**Issue**: Functions handle both computation AND file I/O

**Example from lst_covarep_vq()**:
```r
lst_covarep_vq <- function(...) {
  # Validation (10 lines)
  # Audio loading (20 lines)
  # DSP computation (80 lines)
  # JSTF file writing (50 lines)  # <-- Mixed concern
  # Result formatting (10 lines)
}
```

**Impact**: Harder to test each concern independently

**Recommendation**: Separate I/O from computation

#### 4. Inconsistent Error Context

**Issue**: Some functions provide detailed error messages, others are generic

**Examples**:
- ✅ Good: `lst_phonet()` - "Invalid phonological class(es): X, Y, Z"
- ⚠️ Generic: `lst_voice_tremorp()` - "Error processing file"

**Recommendation**: Standardize error messaging

#### 5. Parameter Validation Redundancy

**Issue**: Each function validates `toFile`, `explicitExt`, `outputDirectory` independently

**Current**: Each function validates parameters
**Better**: Shared validation function

#### 6. File Path Handling

**Issue**: Inconsistent normalization timing

Some functions normalize early:
```r
origSoundFile <- normalizePath(listOfFiles, mustWork = TRUE)
```

Others normalize late (in loop):
```r
for (i in seq_along(listOfFiles)) {
  file_path <- normalizePath(listOfFiles[i], mustWork = TRUE)
```

**Recommendation**: Normalize consistently at function entry

## Improvement Proposals

### Priority 1: Extract JSTF Writing Helper (HIGH)

**Benefit**: Eliminates ~250 lines of duplication

```r
# New helper function in R/jstf_helpers.R
write_lst_results_to_jstf <- function(results,
                                      file_paths,
                                      beginTime,
                                      endTime,
                                      function_name,
                                      parameters,
                                      explicitExt,
                                      outputDirectory = NULL,
                                      speaker_id = NULL) {
  n_files <- length(results)
  output_paths <- character(n_files)

  for (i in seq_len(n_files)) {
    result <- results[[i]]

    if (!is.null(result) && is.null(result$error)) {
      file_path <- file_paths[i]

      # Get audio metadata
      audio_info <- av::av_media_info(file_path)

      # Calculate time range
      analysis_begin <- beginTime[i]
      analysis_end <- if (endTime[i] > 0) endTime[i] else audio_info$duration

      # Create JSTF object
      json_obj <- create_json_track_obj(
        results = result,
        function_name = function_name,
        file_path = file_path,
        sample_rate = audio_info$audio$sample_rate,
        audio_duration = audio_info$duration,
        beginTime = analysis_begin,
        endTime = analysis_end,
        parameters = parameters
      )

      # Determine output path
      if (!is.null(speaker_id) && !is.null(speaker_id[i])) {
        base_name <- as.character(speaker_id[i])
      } else {
        base_name <- tools::file_path_sans_ext(basename(file_path))
      }

      out_dir <- if (is.null(outputDirectory)) dirname(file_path) else outputDirectory
      output_path <- file.path(out_dir, paste0(base_name, ".", explicitExt))

      # Write to file
      write_json_track(json_obj, output_path)
      output_paths[i] <- output_path
    } else {
      output_paths[i] <- NA_character_
    }
  }

  # Return single path or vector
  if (n_files == 1) output_paths[1] else output_paths
}
```

**Usage in functions**:
```r
# Before (50 lines)
if (toFile) {
  output_paths <- character(n_files)
  for (i in seq_len(n_files)) { ... }
  return(invisible(output_paths))
}

# After (5 lines)
if (toFile) {
  output_paths <- write_lst_results_to_jstf(
    results = results,
    file_paths = listOfFiles,
    beginTime = beginTime,
    endTime = endTime,
    function_name = "lst_covarep_vq",
    parameters = list(f0_provided = !is.null(f0), gci_provided = !is.null(gci)),
    explicitExt = explicitExt,
    outputDirectory = outputDirectory
  )
  return(invisible(output_paths))
}
```

**Impact**: Reduces code from ~300 lines to ~50 lines across 6 functions

### Priority 2: Define Package Constants (MEDIUM)

```r
# New file: R/constants.R

# Audio processing constants
PHONET_SAMPLE_RATE <- 16000L  # Phonet requirement
DEFAULT_CHANNELS <- 1L
MS_TO_SECONDS <- 1000

# JSTF constants
JSTF_FORMAT_VERSION <- "1.0"

# Time conversion helpers
ms_to_sec <- function(ms) ms / MS_TO_SECONDS
sec_to_ms <- function(sec) sec * MS_TO_SECONDS
```

### Priority 3: Standardize Parameter Validation (MEDIUM)

```r
# New file: R/validation_helpers.R

validate_jstf_parameters <- function(toFile, explicitExt, outputDirectory) {
  # Validate toFile
  if (!is.logical(toFile) || length(toFile) != 1) {
    stop("toFile must be a single logical value (TRUE/FALSE)", call. = FALSE)
  }

  # Validate explicitExt
  if (!is.character(explicitExt) || length(explicitExt) != 1) {
    stop("explicitExt must be a single character string", call. = FALSE)
  }

  if (nchar(explicitExt) == 0) {
    stop("explicitExt cannot be empty", call. = FALSE)
  }

  # Validate outputDirectory
  if (!is.null(outputDirectory)) {
    if (!is.character(outputDirectory) || length(outputDirectory) != 1) {
      stop("outputDirectory must be a single character string or NULL", call. = FALSE)
    }

    if (!dir.exists(outputDirectory)) {
      stop("outputDirectory does not exist: ", outputDirectory, call. = FALSE)
    }
  }

  invisible(TRUE)
}

# Usage in functions
validate_jstf_parameters(toFile, explicitExt, outputDirectory)
```

### Priority 4: Improve Error Context (LOW)

**Standardized error helper**:
```r
# Add to R/utils.R
format_processing_error <- function(file_path, error_msg, context = NULL) {
  base_msg <- sprintf("Error processing %s: %s", basename(file_path), error_msg)
  if (!is.null(context)) {
    base_msg <- paste0(base_msg, " (", context, ")")
  }
  base_msg
}

# Usage
tryCatch({
  # processing
}, error = function(e) {
  warning(format_processing_error(file_path, e$message, "IAIF computation"))
})
```

### Priority 5: Refactor for Testability (LOW)

**Separate computation from I/O**:
```r
# Internal computation function
.compute_covarep_vq <- function(audio_data, f0, gci, ...) {
  # Pure computation - no I/O
  # Returns results list
}

# Public interface with I/O
lst_covarep_vq <- function(listOfFiles, ..., toFile = FALSE) {
  # Load audio
  audio_data <- av_load_for_python(...)

  # Compute (testable separately)
  results <- .compute_covarep_vq(audio_data, f0, gci, ...)

  # Handle output
  if (toFile) {
    write_lst_results_to_jstf(...)
  } else {
    return(results)
  }
}
```

## Implementation Roadmap

### Phase 1: Extract JSTF Helper (1-2 hours)
- [ ] Create `R/jstf_helpers.R`
- [ ] Implement `write_lst_results_to_jstf()`
- [ ] Add unit tests for helper
- [ ] Refactor 6 functions to use helper
- [ ] Verify all functions still work
- [ ] Commit: "refactor: Extract JSTF writing logic to helper function"

### Phase 2: Add Constants (30 min)
- [ ] Create `R/constants.R`
- [ ] Define package constants
- [ ] Replace magic numbers in code
- [ ] Update documentation
- [ ] Commit: "refactor: Replace magic numbers with named constants"

### Phase 3: Standardize Validation (1 hour)
- [ ] Create `R/validation_helpers.R`
- [ ] Implement `validate_jstf_parameters()`
- [ ] Add to all 6 functions
- [ ] Add unit tests
- [ ] Commit: "refactor: Standardize JSTF parameter validation"

### Phase 4: Error Messages (30 min)
- [ ] Implement `format_processing_error()`
- [ ] Update error handling in 6 functions
- [ ] Test error scenarios
- [ ] Commit: "refactor: Improve error message consistency"

### Phase 5: Testability (optional, 2-3 hours)
- [ ] Extract pure computation functions
- [ ] Add comprehensive unit tests
- [ ] Refactor functions to use internal helpers
- [ ] Commit: "refactor: Separate computation from I/O for better testability"

## Metrics

### Current State
- **Total LOC**: ~1,100 (6 functions × ~180 avg)
- **Duplicated LOC**: ~300 (JSTF writing blocks)
- **Test Coverage**: Unknown (no tests for JSTF integration yet)
- **Cyclomatic Complexity**: Medium (nested conditionals in file writing)

### After Improvements
- **Total LOC**: ~600 (45% reduction)
- **Duplicated LOC**: ~50 (83% reduction)
- **Test Coverage**: 80%+ (with new unit tests)
- **Cyclomatic Complexity**: Low (logic extracted to helpers)

## Risk Assessment

### Low Risk Improvements
- ✅ Extract JSTF helper (well-tested pattern)
- ✅ Add constants (no logic change)
- ✅ Standardize validation (defensive programming)

### Medium Risk Improvements
- ⚠️ Error message refactoring (test error scenarios)
- ⚠️ Computation extraction (requires thorough testing)

## Conclusion

The JSTF Phase 2 integration is **functionally correct** and follows a **consistent pattern**, but has **significant code duplication** that can be reduced by 80%+ through helper function extraction.

**Recommended Action**: Implement Priority 1 (JSTF helper extraction) immediately to reduce technical debt before additional functions are integrated.

**Overall Grade**: B+ (Good, with room for refactoring)
