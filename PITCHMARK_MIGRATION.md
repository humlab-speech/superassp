# ESTK Pitchmark Migration to Protoscribe

## Notice

The `trk_pitchmark()` function has been **migrated to protoscribe** as `draft_pitchmark()`.

## Why?

**Pitchmarks are EVENT annotations**, not DSP measurements:
- Pitchmarks mark discrete time points (glottal closure instants)
- They are phonetic event markers, like VOT boundaries or pitch targets
- superassp provides DSP measurements at regular intervals (formants, spectrograms, etc.)
- protoscribe provides draft annotations for corpus transcription

The function conceptually belongs in protoscribe, where it follows the `draft_` function pattern and integrates seamlessly with reindeer workflows.

## Migration Guide

### Old Usage (superassp - deprecated)

```r
library(superassp)

# Write SSFF files
trk_pitchmark("recording.egg")
trk_pitchmark("recording.wav", fill = TRUE)

# Get data frame
result <- trk_pitchmark("recording.egg", toFile = FALSE)
```

### New Usage (protoscribe - recommended)

```r
library(reindeer)
library(protoscribe)

# Load corpus
corp <- corpus("path/to/database_emuDB")
files <- signal_files(corp)

# Generate pitchmark suggestions
pm <- draft_pitchmark(
  audio_files = files$full_path,
  session_names = files$session,
  bundle_names = files$bundle
)

# Validate and apply (reindeer workflow)
pm <- assess(pm)
prepare(pm)
transcribe(pm)

# Query annotations
data <- ask_for(corp, "Pitchmark == *")
```

## Key Differences

| Feature | superassp::trk_pitchmark() | protoscribe::draft_pitchmark() |
|---------|---------------------------|-------------------------------|
| **Output** | SSFF files or data frames | EventSuggestion objects |
| **Workflow** | Standalone file processing | Integrated with reindeer |
| **Interface** | File-based (listOfFiles) | Corpus-based (sessions/bundles) |
| **Pattern** | trk_ (track/measurement) | draft_ (annotation suggestion) |
| **Integration** | Independent | Part of annotation workflow |

## Deprecation Timeline

- **2025-10-25**: Function migrated to protoscribe (protoscribe commit 17c0649)
- **Current**: Deprecated in superassp with warning message
- **Future**: Will be removed from superassp after deprecation period

## Backwards Compatibility

The `trk_pitchmark()` function **remains available in superassp** for backwards compatibility. It will show a deprecation warning but continues to work as before.

## Installation

```r
# Install protoscribe (for new code)
remotes::install_github("humlab-speech/protoscribe")

# superassp still works (for existing code)
library(superassp)
trk_pitchmark("file.egg")  # Shows deprecation warning
```

## Benefits of Migration

The new `protoscribe::draft_pitchmark()` provides:
- ✓ Integration with reindeer annotation workflows
- ✓ EventSuggestion objects for assessment and validation
- ✓ Consistent with other draft_ functions (draft_hgci, draft_periods, draft_vot)
- ✓ Universal audio format support via av package
- ✓ In-memory processing (no temp files)
- ✓ Cleaner corpus-based interface

## Technical Details

Both functions use the same underlying C++ implementation (ESTK pitchmark algorithm) but provide different R interfaces:
- **superassp**: File-oriented interface, writes SSFF files
- **protoscribe**: Corpus-oriented interface, returns annotations

The ESTK source code is shared between both packages.

## References

**Protoscribe repository**: https://github.com/humlab-speech/protoscribe

**Migration commit**: protoscribe@17c0649

**Algorithm source**: Edinburgh Speech Tools
- Macon, M. W., & Taylor, P. (1997). Pitchmark algorithm
- Centre for Speech Technology Research, University of Edinburgh

## Questions?

For issues or questions about the migration, please file an issue at:
- Protoscribe: https://github.com/humlab-speech/protoscribe/issues
- Superassp: https://github.com/humlab-speech/superassp/issues
