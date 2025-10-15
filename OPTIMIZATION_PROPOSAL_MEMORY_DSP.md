# Optimization Proposal: True In-Memory DSP Processing

## Current Problem

The current implementation has unnecessary disk I/O even when processing from memory:

### Current Flow (Inefficient)
```
Media File (any format)
    ↓ av::read_audio_bin()
AsspDataObj (R, in memory)
    ↓ write.AsspDataObj() ← UNNECESSARY DISK WRITE
Temporary WAV file (disk)
    ↓ performAssp / asspFOpen() ← UNNECESSARY DISK READ
DOBJ (C, in memory)
    ↓ computeXXX()
Result DOBJ
    ↓ dobj2AsspDataObj()
Result AsspDataObj (R)
```

### Issues
1. **Unnecessary disk I/O**: Writing to temp file and reading it back
2. **Performance overhead**: File operations are slow
3. **Temporary file management**: Need to track and cleanup temp files
4. **Code complexity**: Extra error handling for file operations

## Solution: Direct Memory Processing

The infrastructure already exists! We just need to connect the pieces:

### Proposed Flow (Efficient)
```
Media File (any format)
    ↓ av::read_audio_bin()
AsspDataObj (R, in memory)
    ↓ sexp2dobj() ← ALREADY EXISTS!
DOBJ (C, in memory)
    ↓ computeXXX()
Result DOBJ
    ↓ dobj2AsspDataObj()
Result AsspDataObj (R)
```

### Existing Infrastructure

**C Functions (src/dataobj.c):**
- ✅ `sexp2dobj()` (line 595): Converts R AsspDataObj → C DOBJ
- ✅ `dobj2AsspDataObj()` (line 155): Converts C DOBJ → R AsspDataObj
- ✅ `addTrackData()` (line 299): Adds data from R matrix to DOBJ

**Current Usage:**
- `dobj2AsspDataObj()` is used in `performAssp.c:1160` to return results
- `sexp2dobj()` is used in `writeDObj_()` (dataobj.c:280) for file writing

## Implementation Plan

### 1. Create New C Function: `performAsspMemory()`

Add to `src/performAssp.c`:

```c
/*
 * Variant of performAssp that accepts a DOBJ directly instead of a file path.
 * This enables true in-memory processing without temporary files.
 */
SEXP performAsspMemory(SEXP args){
  SEXP el, audioDobj, res;
  const char *name;
  AOPTS OPTS;
  AOPTS *opt = &OPTS;
  W_OPT *wrasspOptions;
  A_F_LIST *anaFunc = funclist;
  int toFile = 0;  // Memory mode always returns results
  DOBJ *inPtr, *outPtr;

  args = CDR(args);  // skip function name

  // First element is AsspDataObj (not file path!)
  audioDobj = CAR(args);
  args = CDR(args);

  // Convert R AsspDataObj to C DOBJ
  inPtr = sexp2dobj(audioDobj);
  if (inPtr == NULL)
    error("Failed to convert AsspDataObj to DOBJ");

  // Second element is assp function name (same as performAssp)
  name = isNull(TAG(args)) ? "" : CHAR(PRINTNAME(TAG(args)));
  if (strcmp(name, "fname") != 0)
    error("Second argument must be named 'fname'");

  el = CAR(args);
  while (anaFunc->funcNum != AF_NONE) {
    if (strcmp(CHAR(STRING_ELT(el, 0)), anaFunc->fName) == 0)
      break;
    anaFunc++;
  }
  if (anaFunc->funcNum == AF_NONE)
    error("Invalid analysis function");

  // Generate default settings
  if ((anaFunc->setFunc)(opt) == -1)
    error("%s", getAsspMsg(asspMsgNum));

  args = CDR(args);

  // Parse options (same logic as performAssp)
  for (int i = 0; args != R_NilValue; i++, args = CDR(args)) {
    name = isNull(TAG(args)) ? "" : CHAR(PRINTNAME(TAG(args)));
    wrasspOptions = anaFunc->options;
    while (wrasspOptions->name != NULL) {
      if (strcmp(wrasspOptions->name, name) == 0)
        break;
      wrasspOptions++;
    }
    if (wrasspOptions->name == NULL)
      error("Invalid option %s for ASSP analysis %s.", name, anaFunc->fName);

    el = CAR(args);
    // ... (copy option parsing from performAssp lines 565-1032)
  }

  // Run the analysis function directly on in-memory DOBJ
  outPtr = (anaFunc->compProc)(inPtr, opt, (DOBJ *) NULL);
  if (outPtr == NULL) {
    freeDObj(inPtr);
    error("%s", getAsspMsg(asspMsgNum));
  }

  // Input no longer needed
  freeDObj(inPtr);

  // Convert result to R AsspDataObj
  PROTECT(res = dobj2AsspDataObj(outPtr));
  freeDObj(outPtr);
  UNPROTECT(1);

  return res;
}
```

### 2. Register New Function in `src/superassp_init.c`

```c
{"performAsspMemory", (DL_FUNC) &performAsspMemory, -1},
```

### 3. Update R Function: `processMediaFiles_LoadAndProcess()`

Modify `R/av_helpers.R` lines 254-283:

```r
# OLD (current inefficient version):
if (is_native[i] && bt == 0.0 && et == 0.0) {
  input_file <- file_path
  use_temp <- FALSE
} else {
  use_temp <- TRUE
  audio_obj <- av_to_asspDataObj(...)
  temp_wav <- tempfile(fileext = ".wav")
  write.AsspDataObj(audio_obj, temp_wav)  # ← REMOVE THIS
  temp_files <- c(temp_files, temp_wav)
  input_file <- temp_wav
  bt <- 0.0
  et <- 0.0
}

# Call performAssp
externalRes[[i]] <- .External(
  "performAssp", input_file,  # ← File path
  fname = fname,
  ...
)

# NEW (optimized version):
if (is_native[i] && bt == 0.0 && et == 0.0) {
  # Native file, no time window - process directly from disk
  externalRes[[i]] <- .External(
    "performAssp", file_path,
    fname = fname,
    beginTime = bt,
    endTime = et,
    toFile = toFile,
    progressBar = NULL,
    ...,
    PACKAGE = "superassp"
  )
} else {
  # Non-native or time-windowed - process from memory
  audio_obj <- av_to_asspDataObj(
    file_path,
    start_time = bt,
    end_time = if(et == 0.0) NULL else et,
    target_sample_rate = NULL
  )

  # Call new memory-based function
  externalRes[[i]] <- .External(
    "performAsspMemory", audio_obj,  # ← AsspDataObj directly!
    fname = fname,
    # beginTime/endTime already applied in av_to_asspDataObj
    toFile = toFile,
    progressBar = NULL,
    ...,
    PACKAGE = "superassp"
  )
}
```

### 4. Simplify Helper Functions

Update `rmsana_memory()` and `process_media_file()` in `R/av_helpers.R`:

```r
# OLD:
rmsana_memory <- function(audio_obj, ...) {
  temp_wav <- tempfile(fileext = ".wav")
  on.exit(unlink(temp_wav), add = TRUE)
  write.AsspDataObj(audio_obj, temp_wav)
  result <- rmsana(temp_wav, toFile = FALSE, ...)
  return(result)
}

# NEW:
rmsana_memory <- function(audio_obj, ...) {
  if (!inherits(audio_obj, "AsspDataObj")) {
    stop("audio_obj must be an AsspDataObj")
  }

  # Call directly without temp files!
  result <- .External(
    "performAsspMemory", audio_obj,
    fname = "rmsana",
    toFile = FALSE,
    ...,
    PACKAGE = "superassp"
  )

  return(result)
}
```

## Performance Benefits

### Expected Improvements

1. **I/O Elimination**:
   - Removes 1 write operation per file
   - Removes 1 read operation per file
   - Removes temp file creation/deletion overhead

2. **Speed Improvement**:
   - Conservative estimate: 20-40% faster for small files
   - Large files (>10MB): 50-70% faster
   - Video files with audio extraction: 60-80% faster

3. **Resource Usage**:
   - No temp files → less disk I/O
   - No file system overhead
   - Reduced cleanup code paths

### Benchmarking Plan

```r
library(bench)
library(superassp)

test_file <- "test_video.mp4"

bench::mark(
  current = {
    # Current: av → write → read → process
    audio_obj <- av_to_asspDataObj(test_file)
    temp_wav <- tempfile(fileext = ".wav")
    write.AsspDataObj(audio_obj, temp_wav)
    result <- rmsana(temp_wav, toFile = FALSE)
    unlink(temp_wav)
  },
  optimized = {
    # New: av → process directly
    audio_obj <- av_to_asspDataObj(test_file)
    result <- rmsana_memory(audio_obj)
  },
  check = FALSE,
  iterations = 100
)
```

## Backwards Compatibility

✅ **Fully backwards compatible**:
- Existing `performAssp()` unchanged
- New `performAsspMemory()` is additional function
- All existing R functions work as before
- Only internal `processMediaFiles_LoadAndProcess()` changes

## Testing Plan

1. **Unit tests**: Verify `performAsspMemory()` produces identical results to `performAssp()`
2. **Integration tests**: Test all DSP functions (acfana, rfcana, rmsana, zcrana, forest, etc.)
3. **Format tests**: Test with various input formats (MP4, M4A, OGG, FLAC)
4. **Edge cases**: Empty files, single-sample files, very long files
5. **Memory tests**: Verify no memory leaks with large file processing
6. **Benchmark tests**: Measure performance improvements

## Implementation Steps

1. ✅ **Phase 1**: Create proposal document (this file)
2. **Phase 2**: Implement `performAsspMemory()` in C
3. **Phase 3**: Update `processMediaFiles_LoadAndProcess()` in R
4. **Phase 4**: Create comprehensive tests
5. **Phase 5**: Run benchmarks and validate improvements
6. **Phase 6**: Update documentation

## Questions to Resolve

1. **Error handling**: How should `performAsspMemory()` handle DOBJ conversion errors?
2. **Progress bars**: Can progress bars work with in-memory processing?
3. **Option validation**: Should we validate that beginTime/endTime are not passed (since already applied)?
4. **File info**: Should we preserve the original file path in attributes even for memory processing?

## Alternative Considered

**Direct DOBJ creation**: Instead of AsspDataObj → DOBJ, create DOBJ directly from av data.
- ❌ More complex
- ❌ Duplicates existing code
- ❌ Less maintainable
- ✅ Current approach reuses existing, tested conversion code

## Success Criteria

- ✅ All tests pass
- ✅ No breaking changes to existing API
- ✅ Measurable performance improvement (>20%)
- ✅ No memory leaks
- ✅ Code is well-documented

## References

- `src/dataobj.c`: Conversion functions
- `src/performAssp.c`: DSP analysis dispatcher
- `R/av_helpers.R`: Media file processing
- `tests/test_av_integration.R`: Integration tests

---

**Status**: Proposal
**Date**: 2025-10-14
**Author**: Based on code analysis and optimization opportunity identified by user
