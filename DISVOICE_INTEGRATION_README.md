# DisVoice Integration for superassp

## Overview

This package now includes optimized speech analysis functions using DisVoice's in-memory Parselmouth processing. These functions are **1.5-2.5x faster** than traditional file-based Praat methods.

---

## New Functions

### Pitch Tracking

**Function**: `trk_dv_f0()`
**File**: `R/ssff_python_dv_f0.R`
**Backend**: Python (DisVoice + Parselmouth)
**Performance**: ~2x faster than file-based methods

```r
# Basic usage
f0_track <- trk_dv_f0("audio.wav")

# With custom parameters
f0_track <- trk_dv_f0(
  "audio.wav",
  frame_shift = 10,       # 10ms = 100 Hz
  min_f0 = 75,
  max_f0 = 600,
  include_voicing = TRUE  # Binary voicing track
)

# Get as data.frame
f0_df <- trk_dv_f0("audio.wav", output_format = "dataframe")
```

**Output**: AsspDataObj with tracks:
- `f0`: Fundamental frequency (Hz)
- `voicing` (optional): Binary voicing decision (1 = voiced, 0 = unvoiced)

---

### Formant Tracking

**Function**: `trk_dv_formants()`
**File**: `R/ssff_python_dv_formants.R`
**Backend**: Python (DisVoice + Parselmouth)
**Performance**: ~1.5x faster than file-based methods

```r
# Basic usage
formant_track <- trk_dv_formants("audio.wav")

# With custom parameters
formant_track <- trk_dv_formants(
  "audio.wav",
  frame_shift = 5,         # 5ms = 200 Hz
  window_size = 25,         # 25ms window
  max_formants = 5,
  max_formant_freq = 5500   # Max formant frequency
)

# Access individual formants
plot(formant_track$tracks[, "F1"], type = "l", col = "blue")
lines(formant_track$tracks[, "F2"], col = "red")
```

**Output**: AsspDataObj with tracks:
- `F1`: First formant frequency (Hz)
- `F2`: Second formant frequency (Hz)
- `F3`: Third formant frequency (Hz)
- `F4`: Fourth formant frequency (Hz)

**Note**: Returns F1-F4 (4 formants), unlike file-based methods that typically return only F1-F2.

---

## Installation

### Step 1: Install/Update superassp

```r
# From CRAN (when updated)
install.packages("superassp")

# Or from GitHub
# devtools::install_github("humlab-speech/superassp")
```

### Step 2: Check DisVoice Availability

```r
library(superassp)

# Check if Python and DisVoice are available
has_disvoice_support()  # Returns TRUE/FALSE
```

### Step 3: Install Python Dependencies (if needed)

If `has_disvoice_support()` returns `FALSE`:

```r
# Install Python packages (parselmouth + numpy)
install_disvoice_python()

# Verify installation
has_disvoice_support()  # Should now return TRUE
```

**What gets installed:**
- `praat-parselmouth>=0.4.0` - Parselmouth (Praat Python interface)
- `numpy>=1.21.0` - NumPy (numerical computing)

---

## File Naming Convention

DisVoice functions follow superassp's standard naming convention:

**File names**: `ssff_{backend}_{library}_{method}.R`
- `ssff` - Returns SSFF-compatible data structure
- `python` - Python backend via reticulate
- `dv` - DisVoice library
- `{method}` - Specific analysis method (f0, formants, etc.)

**Examples**:
- `ssff_python_dv_f0.R` - F0 tracking with DisVoice
- `ssff_python_dv_formants.R` - Formant tracking with DisVoice

**Function names**: `trk_dv_{method}`
- `trk` - Track function (equal-interval output)
- `dv` - DisVoice library
- `{method}` - Analysis method

**Examples**:
- `trk_dv_f0()` - F0 tracking
- `trk_dv_formants()` - Formant tracking

---

## Technical Details

### In-Memory Processing Pipeline

```
Audio File (.wav, .mp3, .flac, etc.)
    ↓
av::read_audio_bin()  [in-memory loading]
    ↓
INT32 Audio Data (R)
    ↓
as_numpy_audio()  [single conversion]
    ↓
Float64 NumPy Array (Python)
    ↓
DisVoice extract_*()  [no temp files!]
    ↓
F0/Formants + Times (Python)
    ↓
numpy_to_r()  [auto-conversion]
    ↓
AsspDataObj (R)
```

**Key optimization**: Zero temporary files throughout entire pipeline!

### Output Formats

All DisVoice functions support three output formats:

1. **"AsspDataObj"** (default) - Compatible with superassp workflows
2. **"dataframe"** - For dplyr/tidyverse workflows
3. **"list"** - For custom processing

```r
# As AsspDataObj (default)
obj <- trk_dv_f0("audio.wav")

# As data.frame
df <- trk_dv_f0("audio.wav", output_format = "dataframe")

# As list
lst <- trk_dv_f0("audio.wav", output_format = "list")
```

### Lazy Loading

DisVoice Python modules are only loaded when first used:
- No startup overhead if DisVoice functions not called
- Cached after first load for performance
- Graceful fallback if Python unavailable

---

## Performance Comparisons

Based on DisVoice benchmarks (2-second audio file, 16 kHz):

| Operation | File-Based | DisVoice | Speedup |
|-----------|------------|----------|---------|
| F0 extraction | 7.0 ms | 4.7 ms | **1.49x** |
| Formant extraction | 17.5 ms | 11.8 ms | **1.48x** |
| Combined (F0 + voicing) | ~25 ms | ~10 ms | **2.5x** |

**Additional benefits**:
- No temporary files created
- Cleaner `/tmp` directory
- Reduced disk I/O overhead

---

## Helper Functions

### Check Availability

```r
# Check if DisVoice Python support is available
has_disvoice_support()  # TRUE or FALSE
```

### Install Dependencies

```r
# Install Python packages
install_disvoice_python()

# With custom Python environment
install_disvoice_python(method = "conda")
```

### Initialization

```r
# Lazy initialization (automatic)
# DisVoice modules loaded on first trk_dv_*() call

# Manual initialization (rarely needed)
init_disvoice()  # Internal function
```

---

## Comparison with Existing Methods

### F0 Tracking

| Function | Backend | Speed | Voicing | Notes |
|----------|---------|-------|---------|-------|
| `trk_dv_f0()` | Python/Parselmouth | ★★★ | Yes | In-memory, fast |
| `trk_rapt()` | C++/SPTK | ★★ | No | File-based |
| `trk_swipe()` | C++/SPTK | ★★ | No | File-based |
| `trk_ppitch()` | Python/Parselmouth | ★★ | No | File-based |

### Formant Tracking

| Function | Backend | Speed | Formants | Notes |
|----------|---------|-------|----------|-------|
| `trk_dv_formants()` | Python/Parselmouth | ★★★ | F1-F4 | In-memory, fast |
| `trk_pformant()` | Python/Parselmouth | ★★ | F1-F4 | File-based |
| (ASSP formants) | C/ASSP | ★★ | F1-F2 | File-based |

---

## Troubleshooting

### DisVoice Not Available

**Problem**: `has_disvoice_support()` returns `FALSE`

**Solution**:
```r
# Install Python dependencies
install_disvoice_python()
```

### Python Not Found

**Problem**: `Error: Python not found`

**Solution**:
```r
# Install miniconda (recommended)
reticulate::install_miniconda()

# Then install DisVoice
install_disvoice_python()
```

### Import Error

**Problem**: `Error: could not import DisVoice`

**Solution**:
```r
# Check Python path
reticulate::py_config()

# Reinstall
install_disvoice_python()

# Or specify Python version
reticulate::use_python("/path/to/python")
install_disvoice_python()
```

---

## Advanced Usage

### Batch Processing

```r
# Process multiple files
audio_files <- list.files("audio_dir", pattern = "\\.wav$", full.names = TRUE)

f0_results <- lapply(audio_files, function(file) {
  trk_dv_f0(file, output_format = "dataframe")
})

# Combine results
library(data.table)
f0_combined <- rbindlist(f0_results, idcol = "file_id")
```

### Custom Parameters by Speaker

```r
# Male speakers
f0_male <- trk_dv_f0(
  "male_speaker.wav",
  min_f0 = 60,
  max_f0 = 300
)

formants_male <- trk_dv_formants(
  "male_speaker.wav",
  max_formant_freq = 5000
)

# Female speakers
f0_female <- trk_dv_f0(
  "female_speaker.wav",
  min_f0 = 100,
  max_f0 = 500
)

formants_female <- trk_dv_formants(
  "female_speaker.wav",
  max_formant_freq = 5500
)
```

---

## Contributing

To add new DisVoice-based functions:

1. **Create Python module** in `inst/python/DisVoice/`
2. **Create R wrapper** in `R/ssff_python_dv_{method}.R`
3. **Follow naming convention**: `trk_dv_{method}`
4. **Add documentation** with `@export` tag
5. **Add tests** in `tests/testthat/`

---

## References

- **DisVoice**: https://github.com/jcvasquezc/DisVoice
- **Parselmouth**: https://parselmouth.readthedocs.io/
- **superassp**: https://github.com/humlab-speech/superassp

---

## License

DisVoice integration follows DisVoice's MIT License.
superassp package license applies to R wrapper code.
