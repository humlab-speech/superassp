# Deprecation Candidates - superassp v0.8.9

**Date**: 2025-10-30
**Status**: Proposal for Review

This document identifies functions that are candidates for deprecation due to:
1. Redundant functionality (better alternatives available)
2. Legacy implementations with modern replacements
3. Low usage or specialized use cases

---

## IMMEDIATE DEPRECATION CANDIDATES

### 1. trk_reaper_pm() - Redundant Pitchmark Function

**Status**: 🔴 **RECOMMEND IMMEDIATE DEPRECATION**

**Issue**:
- Functionality completely duplicated by `trk_pitchmark()`
- `trk_pitchmark()` is more modern, faster, and better maintained
- Creates confusion for users choosing between functions

**Migration Path**:
```r
# Old code
pm <- trk_reaper_pm("audio.wav")

# New code
pm <- trk_pitchmark("audio.wav", method = "reaper")
```

**Timeline**:
- v0.8.10: Add deprecation warning
- v0.9.0: Mark as deprecated in documentation
- v1.0.0: Remove function

**Affected Users**: Low (specialized function, limited usage)

---

## EVALUATE BASED ON USAGE METRICS

### 2. STRAIGHT Library Functions (3 functions)

**Functions**:
- `trk_straight_f0()` - STRAIGHT pitch tracking
- `trk_straight_spec()` - STRAIGHT spectral analysis
- `trk_straight_synth()` - STRAIGHT vocoder synthesis

**Status**: ❓ **EVALUATE USAGE BEFORE DECISION**

**Issues**:
- Legacy Python STRAIGHT library (limited maintenance)
- File-based processing (not in-memory)
- Modern alternatives available:
  - F0: Use `trk_dio()`, `trk_harvest()`, `trk_rapt()`
  - Spectral: Use `trk_cssSpectrum()`, `trk_lpsSpectrum()`
  - Synthesis: Use WORLD (`dio()` + `d4c()`)

**Alternative**: WORLD Algorithm Suite
```r
# STRAIGHT F0
old_f0 <- trk_straight_f0("audio.wav")

# Modern alternative: WORLD DIO or Harvest
new_f0 <- trk_dio("audio.wav")         # Fast, robust
new_f0 <- trk_harvest("audio.wav")     # High quality
```

**Decision Criteria**:
- Survey user base for STRAIGHT usage
- If usage <5%, deprecate
- If usage >5%, migrate to in-memory processing

**Timeline** (if deprecation approved):
- v0.8.10: Collect usage metrics
- v0.9.0: Add deprecation warnings if usage <5%
- v1.0.0: Remove if migration to alternatives is successful

**Affected Users**: Unknown (need metrics)

---

### 3. trk_aperiodicities() - Legacy Aperiodicity Function

**Status**: ❓ **EVALUATE USAGE BEFORE DECISION**

**Issue**:
- File-based processing
- Redundant with modern `trk_d4c()` (WORLD aperiodicity)

**Migration Path**:
```r
# Old code
ap <- trk_aperiodicities("audio.wav")

# New code (WORLD D4C)
ap <- trk_d4c("audio.wav")
```

**Decision Criteria**:
- Check usage in downstream packages
- If <5% usage, deprecate
- If >5% usage, migrate to in-memory

**Timeline** (if deprecation approved):
- v0.9.0: Add deprecation warning
- v1.0.0: Remove function

**Affected Users**: Unknown (need metrics)

---

## NOT RECOMMENDED FOR DEPRECATION (But Document Trade-offs)

### OpenSMILE Python Implementations

**Functions**:
- `lst_GeMAPS_python()`
- `lst_eGeMAPS_python()`
- `lst_emobase_python()`
- `lst_ComParE_2016_python()`

**Status**: ✅ **KEEP AS FALLBACK**

**Rationale**:
- Python versions serve as fallback for systems without C++ build tools
- Some users may prefer Python-only workflows
- Both implementations fully compliant with modern workflow
- Dispatcher functions (`lst_GeMAPS()`, etc.) default to C++ automatically

**Performance Trade-off**:
- C++: ~50ms for 3s audio (**Recommended**)
- Python: ~275ms for 3s audio (5.5x slower, but acceptable fallback)

**Recommendation**:
- ✅ Keep both implementations
- ✅ Default to C++ (`use_cpp = TRUE`)
- 📖 Document performance differences in user guide
- ⚠️ Add note in Python function docs: "Consider using C++ version for better performance"

---

## MIGRATION REQUIRED (Not Deprecation)

The following functions should be **migrated**, not deprecated:

### High Priority Migration (5 functions)

| Function | Issue | Action |
|----------|-------|--------|
| `trk_yaapt` | File-based | Migrate to `av::read_audio_bin()` |
| `trk_snackp` | File-based | Migrate to `av::read_audio_bin()` |
| `trk_snackf` | File-based | Migrate to `av::read_audio_bin()` |
| `trk_formantp` | Temp files | Migrate to `av_load_for_parselmouth()` |
| `trk_formantpathp` | Temp files | Migrate to `av_load_for_parselmouth()` |

These provide unique functionality not available elsewhere and should be modernized, not removed.

---

## USAGE METRICS NEEDED

To make informed deprecation decisions, collect metrics on:

1. **Function call frequency** (via package usage telemetry, if available)
2. **Dependency analysis** (check reverse dependencies on CRAN)
3. **User survey** (query superassp-users mailing list)
4. **GitHub issues** (search for function mentions)

### Suggested Survey Questions

```
Dear superassp users,

We are considering deprecating some legacy functions in favor of modern alternatives.
Please let us know if you actively use any of these:

1. trk_straight_f0, trk_straight_spec, trk_straight_synth (STRAIGHT library)
2. trk_reaper_pm (use trk_pitchmark instead)
3. trk_aperiodicities (use trk_d4c instead)

If you use these functions, please respond with:
- Which functions?
- Use case/purpose?
- Would migration to alternatives be acceptable?

Thank you!
```

---

## DEPRECATION WORKFLOW

### Standard Deprecation Process

1. **v0.X.Y** (Decision version):
   - Collect usage metrics
   - Make deprecation decision
   - Add `.Deprecated()` call with message

2. **v0.X+1.0** (Warning version):
   - Function issues deprecation warning
   - Documentation marked as "Deprecated"
   - Migration guide published
   - Minimum 6 months before removal

3. **v1.0.0** (Removal version):
   - Function removed from package
   - Error message points to alternative
   - NEWS.md documents removal

### Example Deprecation Code

```r
trk_reaper_pm <- function(...) {
  .Deprecated(
    "trk_pitchmark",
    package = "superassp",
    msg = paste(
      "trk_reaper_pm() is deprecated and will be removed in v1.0.0.",
      "Use trk_pitchmark(..., method = 'reaper') instead."
    )
  )

  # Forward to new function
  trk_pitchmark(..., method = "reaper")
}
```

---

## SUMMARY TABLE

| Function | Status | Reason | Alternative | Priority |
|----------|--------|--------|-------------|----------|
| `trk_reaper_pm` | 🔴 Deprecate | Redundant | `trk_pitchmark()` | High |
| `trk_straight_f0` | ❓ Evaluate | Legacy library | `trk_dio()`, `trk_harvest()` | Medium |
| `trk_straight_spec` | ❓ Evaluate | Legacy library | `trk_cssSpectrum()` | Medium |
| `trk_straight_synth` | ❓ Evaluate | Legacy library | WORLD vocoder | Low |
| `trk_aperiodicities` | ❓ Evaluate | Redundant | `trk_d4c()` | Low |
| OpenSMILE Python | ✅ Keep | Fallback needed | Use C++ version | N/A |

---

## RECOMMENDATIONS

### Immediate Actions (v0.8.10)

1. ✅ Add deprecation warning to `trk_reaper_pm()`
2. 📊 Implement usage tracking for STRAIGHT functions
3. 📧 Send survey to user community
4. 📖 Document OpenSMILE C++ vs Python trade-offs

### v0.9.0 Actions

1. ⚠️ Mark `trk_reaper_pm` as officially deprecated
2. ❓ Based on metrics, decide on STRAIGHT/aperiodicities
3. 📝 Publish migration guide for deprecated functions
4. 🧪 Beta test with key users

### v1.0.0 Actions

1. ❌ Remove `trk_reaper_pm`
2. ❌ Remove other functions if deprecation approved
3. 🎉 Celebrate 100% modern workflow compliance!

---

## MIGRATION GUIDES

### From trk_reaper_pm to trk_pitchmark

```r
# Before
pm <- trk_reaper_pm("audio.wav", minF = 50, maxF = 400)

# After
pm <- trk_pitchmark("audio.wav", method = "reaper", minF = 50, maxF = 400)

# Benefits:
# - Unified interface for all pitchmark algorithms
# - Better maintained
# - More options (pda, reaper methods)
```

### From STRAIGHT to WORLD

```r
# STRAIGHT F0 → WORLD DIO/Harvest
# Before
f0 <- trk_straight_f0("audio.wav", minF = 70, maxF = 200)

# After (fast)
f0 <- trk_dio("audio.wav", minF = 70, maxF = 200)

# After (high quality)
f0 <- trk_harvest("audio.wav", minF = 70, maxF = 200)

# Benefits:
# - Actively maintained
# - In-memory processing
# - Part of WORLD vocoder ecosystem
```

### From trk_aperiodicities to trk_d4c

```r
# Before
ap <- trk_aperiodicities("audio.wav")

# After
ap <- trk_d4c("audio.wav")

# Benefits:
# - Modern WORLD algorithm
# - In-memory processing
# - Better integration with WORLD ecosystem
```

---

## CONTACT

For questions or feedback on deprecation decisions:
- **GitHub Issues**: https://github.com/humlab-speech/superassp/issues
- **Maintainer**: fredrik.nylen@umu.se

---

**Last Updated**: 2025-10-30
**Next Review**: After usage metrics collection
