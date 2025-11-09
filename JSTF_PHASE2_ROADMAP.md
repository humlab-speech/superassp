# JSTF Phase 2 Integration Roadmap

**Status**: Ready to begin
**Estimated Duration**: 2-3 weeks
**Prerequisite**: Phase 1 complete ✅

---

## Overview

Phase 2 will integrate JSTF `toFile` support into all 14 existing `lst_*` functions that currently only return in-memory results. This will enable efficient file-based storage for all list-producing DSP functions.

---

## Target Functions (14 total)

### High Priority (Core Voice Quality Functions)

1. **lst_vat()** - Voice Analysis Toolbox (132 measures)
   - File: `R/list_vat.R`
   - Extension: `.vat`
   - Status: No toFile support currently ❌
   - Priority: HIGH (most comprehensive voice quality analysis)

2. **lst_voice_sauce()** - VoiceSauce voice quality (40+ params)
   - File: `R/list_voice_sauce.R`
   - Extension: `.vsj`
   - Status: No toFile support currently ❌
   - Priority: HIGH (widely used voice quality tool)

3. **lst_dysprosody()** - Dysprosody prosodic features (193 features)
   - File: `R/list_dysprosody.R`
   - Extension: `.dyp`
   - Status: No toFile support currently ❌
   - Priority: HIGH (comprehensive prosodic analysis)

4. **lst_voxit()** - Voxit voice/articulation complexity (11 features)
   - File: `R/list_voxit.R`
   - Extension: `.vxt`
   - Status: No toFile support currently ❌
   - Priority: MEDIUM (newer addition)

### Medium Priority (OpenSMILE Feature Sets)

5. **lst_GeMAPS()** - GeMAPS acoustic features (62 features)
   - File: `R/list_cpp_opensmile_gemaps.R` or `R/list_python_opensmile_GeMAPS.R`
   - Extension: `.gem`
   - Status: To be verified
   - Priority: MEDIUM (standard feature set)

6. **lst_eGeMAPS()** - Extended GeMAPS features (88 features)
   - File: `R/list_python_opensmile_eGeMAPS.R`
   - Extension: `.egm`
   - Status: To be verified
   - Priority: MEDIUM (extended feature set)

7. **lst_emobase()** - emobase emotional features (988 features)
   - File: `R/list_cpp_opensmile_emobase.R` or `R/list_python_opensmile_emobase.R`
   - Extension: `.emb`
   - Status: To be verified
   - Priority: MEDIUM (emotion analysis)

8. **lst_ComParE_2016()** - ComParE 2016 features (6373 features)
   - File: `R/list_python_opensmile_ComParE_2016.R`
   - Extension: `.cmp`
   - Status: To be verified
   - Priority: LOW (very large feature set)

### Low Priority (Specialized Functions)

9. **lst_covarep_vq()** - COVAREP voice quality
   - File: To be located
   - Extension: `.cvq`
   - Status: To be verified
   - Priority: LOW

10. **lst_avqip()** - AVQI Acoustic Voice Quality Index
    - File: `R/list_python_pm_pavqi.R`
    - Extension: `.avq`
    - Status: To be verified
    - Priority: LOW (single measure)

11. **lst_dsip()** - Dysphonia Severity Index (7 measures)
    - File: `R/list_python_pm_pdsi.R`
    - Extension: `.dsi`
    - Status: To be verified
    - Priority: LOW

12. **lst_voice_reportp()** - Praat voice report (9 measures)
    - File: `R/list_python_pm_pvoice_report.R`
    - Extension: `.vrp`
    - Status: To be verified
    - Priority: LOW

13. **lst_voice_tremorp()** - Voice tremor analysis (6 measures)
    - File: `R/list_python_pm_pvoice_tremor.R`
    - Extension: `.vtr`
    - Status: To be verified
    - Priority: LOW

14. **lst_phonet()** - Phonological posteriors (18 measures)
    - File: To be located
    - Extension: `.phn`
    - Status: To be verified
    - Priority: LOW

---

## Integration Pattern (Template)

For each function, apply this minimal modification pattern:

### Step 1: Add Parameters

Add three new parameters to the function signature:

```r
lst_function <- function(listOfFiles,
                        # ... existing parameters ...
                        toFile = FALSE,              # NEW
                        explicitExt = "ext",         # NEW (use registered extension)
                        outputDirectory = NULL,      # NEW
                        verbose = TRUE) {
```

### Step 2: Document Parameters

Add roxygen2 documentation for new parameters:

```r
#' @param toFile Logical. If TRUE, write results to JSTF file. Default FALSE.
#' @param explicitExt Character. File extension for output. Default "ext".
#' @param outputDirectory Character. Output directory path. Default NULL (use input directory).
```

### Step 3: Add toFile Logic

Add file writing logic before the return statement:

```r
  # ... existing processing code ...
  results <- your_dsp_processing(audio)

  # NEW: JSTF file writing
  if (toFile) {
    json_obj <- create_json_track_obj(
      results = results,
      function_name = "lst_function",
      file_path = file_path,
      sample_rate = sample_rate,
      audio_duration = audio_duration,
      beginTime = beginTime,
      endTime = endTime,
      parameters = list(...)  # Capture all DSP parameters
    )

    base_name <- tools::file_path_sans_ext(basename(file_path))
    out_dir <- if (is.null(outputDirectory)) dirname(file_path) else outputDirectory
    output_path <- file.path(out_dir, paste0(base_name, ".", explicitExt))

    write_json_track(json_obj, output_path)
    return(invisible(output_path))
  }

  return(results)  # Existing in-memory return
}
```

### Step 4: Set Function Attributes

Add attributes at the end of the file:

```r
attr(lst_function, "ext") <- "ext"
attr(lst_function, "outputType") <- "JSTF"
attr(lst_function, "format") <- "JSON"
```

### Step 5: Update Tests

Add test cases for toFile functionality:

```r
test_that("lst_function works with toFile=TRUE", {
  skip_if_not_installed("superassp")

  test_file <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_file == "", "Test file not found")

  # Test file output
  tmp_dir <- tempdir()
  result <- lst_function(test_file, toFile = TRUE, outputDirectory = tmp_dir)

  expect_true(file.exists(result))

  # Test reading back
  track <- read_track(result)
  expect_s3_class(track, "JsonTrackObj")

  # Test conversion
  df <- as.data.frame(track)
  expect_s3_class(df, "data.frame")
  expect_true("begin_time" %in% names(df))
  expect_true("end_time" %in% names(df))

  # Cleanup
  unlink(result)
})
```

---

## Implementation Schedule

### Week 1: High Priority Functions (4 functions)

**Day 1-2**: `lst_vat()` integration
- Modify function (add toFile support)
- Update tests
- Verify with real data
- Commit and document

**Day 3**: `lst_voice_sauce()` integration
- Same process as above

**Day 4**: `lst_dysprosody()` integration
- Same process as above

**Day 5**: `lst_voxit()` integration
- Same process as above

### Week 2: Medium Priority Functions (4 OpenSMILE functions)

**Day 1**: `lst_GeMAPS()` and `lst_eGeMAPS()`
- Integrate both (similar structure)

**Day 2**: `lst_emobase()` and `lst_ComParE_2016()`
- Integrate both (similar structure)

**Day 3**: Testing and validation
- Run comprehensive tests on all 8 functions
- Fix any issues discovered

**Day 4-5**: Documentation updates
- Update function documentation
- Update CLAUDE.md examples
- Create integration examples

### Week 3: Low Priority Functions + Polish (6 functions)

**Day 1-2**: Specialized functions
- `lst_covarep_vq()`
- `lst_avqip()`
- `lst_dsip()`
- `lst_voice_reportp()`
- `lst_voice_tremorp()`
- `lst_phonet()`

**Day 3**: Comprehensive testing
- Integration tests for all 14 functions
- Performance benchmarks
- File size comparisons

**Day 4**: Documentation finalization
- Update NEWS.md
- Update README.md with examples
- Create Phase 2 completion summary

**Day 5**: Code review and cleanup
- Review all changes
- Ensure consistency across functions
- Prepare for v0.11.0 release

---

## Testing Strategy

### For Each Function

1. **Basic functionality**: toFile=TRUE creates file
2. **Read back**: read_track() successfully reads file
3. **Conversion**: as.data.frame() works correctly
4. **Field preservation**: All measures present in output
5. **Time slicing**: beginTime/endTime parameters work
6. **Batch processing**: Multiple files work correctly
7. **Error handling**: Invalid inputs handled gracefully

### Integration Tests

1. **Cross-function compatibility**: JSTF files from different functions
2. **Merge operations**: Combining outputs from different time slices
3. **Subset operations**: Filtering by time range works
4. **Performance**: File I/O performance acceptable
5. **File size**: Verify space efficiency (99% reduction maintained)

---

## Success Criteria

### Per-Function Checklist

- [ ] toFile parameter added and documented
- [ ] explicitExt parameter added and documented
- [ ] outputDirectory parameter added and documented
- [ ] Function attributes set (ext, outputType, format)
- [ ] JSTF writing logic implemented
- [ ] Tests updated with toFile test cases
- [ ] All tests passing (0 failures)
- [ ] Documentation updated
- [ ] Examples provided
- [ ] Performance verified

### Overall Phase 2 Completion Criteria

- [ ] All 14 functions integrated
- [ ] 100% test pass rate maintained
- [ ] Documentation complete
- [ ] Performance benchmarks completed
- [ ] NEWS.md updated
- [ ] README.md updated with examples
- [ ] Phase 2 summary document created
- [ ] Ready for v0.11.0 release

---

## Risk Mitigation

### Potential Issues

1. **Python function complexity**: Some Python-based functions may have complex return structures
   - **Mitigation**: Use integration example as template, test incrementally

2. **OpenSMILE functions**: May have very large feature sets
   - **Mitigation**: Verify RcppSimdJson can handle large files, benchmark performance

3. **Backward compatibility**: Existing code relies on in-memory returns
   - **Mitigation**: toFile=FALSE is default, no breaking changes

4. **Test coverage**: Some functions may have limited test coverage
   - **Mitigation**: Add comprehensive tests as part of integration

### Rollback Plan

If integration causes issues:
1. Revert changes for problematic function
2. Keep successfully integrated functions
3. Document issues and create separate fix plan
4. Release with partial integration (note in NEWS.md)

---

## Resources Needed

### Development Tools

- R environment with devtools
- All Python dependencies installed
- Test audio files
- Performance benchmarking tools

### Documentation

- Function integration template (available in `R/json_track_integration_example.R`)
- Test template (available in `tests/testthat/test-json-track.R`)
- JSTF specification (available in `JSON_TRACK_FORMAT_SPECIFICATION.md`)

### Review

- Code review after each function integration
- Integration testing after week 1 and week 2
- Final review before Phase 2 completion

---

## Deliverables

### Code

- 14 modified R function files
- 14 updated test files
- Updated NAMESPACE (auto-generated)
- Updated man/*.Rd files (auto-generated)

### Documentation

- Updated CLAUDE.md (integration examples)
- Updated NEWS.md (Phase 2 announcement)
- Updated README.md (usage examples)
- JSTF_PHASE2_COMPLETION_SUMMARY.md (new)

### Testing

- 14 × 10 = 140 new test cases (minimum)
- Performance benchmarks
- Integration test suite

---

## Next Steps (How to Begin)

1. **Verify current status**: Run all existing tests to ensure clean baseline
2. **Choose first function**: Start with `lst_vat()` (most important)
3. **Create feature branch**: `git checkout -b feature/jstf-phase2-integration`
4. **Implement using template**: Follow integration pattern above
5. **Test thoroughly**: Ensure all tests pass before moving to next function
6. **Document progress**: Update this roadmap with completion status

---

**Status**: Ready to begin Phase 2 integration
**Last Updated**: 2025-11-09
**Next Review**: After first 4 functions integrated
