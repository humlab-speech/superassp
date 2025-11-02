# superassp R Package Structure Analysis

## Overview
superassp is a comprehensive speech signal processing package that wraps multiple DSP frameworks (Python, C++, Praat/Parselmouth, MATLAB) into a unified R interface using the wrassp-compatible AsspDataObj structure.

**Package Version:** 0.8.4 (as of 2025-01-27)
**Author:** Fredrik Nylén (Umeå University)
**License:** GPL (>= 2)

---

## 1. PACKAGE STRUCTURE

### Directory Layout
```
superassp/
├── DESCRIPTION              # Package metadata & dependencies
├── NAMESPACE                # Function exports
├── R/                       # 96 R source files (~27,843 lines)
├── src/                     # C++ source code (Rcpp)
├── inst/
│   ├── python/             # 50+ Python DSP modules
│   ├── opensmile/          # OpenSmile C++ integration
│   ├── matlab/             # MATLAB legacy code
│   ├── praat/              # Praat/Parselmouth integration
│   ├── onnx/               # Neural network models
│   └── samples/            # Test audio files
├── man/                    # 220+ documentation files
├── tests/                  # Test suite
└── benchmarking/           # Performance benchmarks
```

### File Counts
- **R Functions:** 96 files
- **Total Lines:** 27,843 lines of R code
- **Exported Functions:** 160+ (see NAMESPACE)
- **Python Modules:** 50+ integrated modules
- **Man Pages:** 220+ documented functions

---

## 2. DESCRIPTION FILE

### Key Dependencies
```
Imports:
  - tidyr, dplyr, purrr (Data manipulation)
  - reticulate (Python integration)
  - av (Audio/video format support)
  - S7 (Modern OOP)
  - Rcpp (C++ integration)
  - wrassp-compatible interface (AsspDataObj)
  - logger, cli (Messaging)
  - digest, uuid (Hashing, unique IDs)
  - R.matlab (MATLAB legacy support)
  - lifecycle (Function deprecation)

Suggests:
  - pbapply, pbmcapply (Progress bars)

Remotes:
  - github::curso-r/torchaudio
  - github::humlab-speech/av
```

### Special Notes
- **S7 Classes:** Modern S7 object system (not S3/S4)
- **C++11 Required:** System requirements specify C++11
- **Rcpp Modules:** Dynamic compilation via LinkingTo: Rcpp
- **Python Mandatory:** reticulate required for many functions

---

## 3. NAMESPACE & EXPORTS

### Main Export Categories

#### A. Track Processing Functions (trk_*)
**60+ exported trk_* functions** organized by DSP backend:

**Python-based:**
- `trk_crepe` - CREPE neural network pitch
- `trk_yin`, `trk_pyin` - YIN/probabilistic YIN
- `trk_swiftf0` - Swift F0 (requires external Python module)
- `trk_deepformants` - Deep learning formant tracking
- `trk_sacc` - SACC voice quality
- `trk_mfcc` - Mel-frequency cepstral coefficients
- `trk_seenc` - Spectral entropy energy cepstral
- `trk_snackp`, `trk_snackf` - SNACK pitch/formants
- `trk_gfmiaif` - GFMIAIF analysis
- `trk_cepstrum`, `trk_cssSpectrum`, `trk_dftSpectrum`, `trk_lpsSpectrum`
- `trk_covarep_iaif`, `trk_covarep_srh` - COVAREP methods
- `trk_vat_srh` - Voice Activity Toolkit
- `trk_excite` - Excitation analysis
- `trk_praat_sauce`, `trk_intensityp`, `trk_pitchp`, `trk_formantpathp`, `trk_formantp` - Praat/Parselmouth

**C++/SPTK-based:**
- `trk_dio`, `trk_harvest`, `trk_rapt`, `trk_swipe`, `trk_reaper` - Pitch trackers
- `trk_d4c` - Spectral envelope

**C/ASSP-based:**
- `trk_acfana` - Autocorrelation analysis
- `trk_forest` - FOREST formant tracking
- `trk_zcrana` - Zero-crossing analysis
- `trk_rmsana` - RMS amplitude
- `trk_pitchmark`, `trk_pitchp` - Pitch marking

**Custom/Hybrid:**
- `trk_creak_union` - Union method for creaky voice detection
- `trk_egg_f0` - EGG (electroglottographic) F0 and Oq
- `trk_formants_tvwlp` - Time-varying weighted LP formants

#### B. Feature Extraction (lst_*)
**15+ feature lists** from external frameworks:
- `lst_ComParE_2016`, `lst_GeMAPS`, `lst_eGeMAPS` - OpenSmile
- `lst_voice_sauce` - Voice Sauce features
- `lst_covarep_vq` - COVAREP voice quality
- `lst_deepformants` - Deep formant features
- `lst_dysprosody` - Dysprosody indices
- `lst_vat`, `lst_voice_analysis` - Voice analysis toolkit
- `lst_avqip`, `lst_dsip` - Acoustic voice quality indices

#### C. Data I/O & Conversion
- `read.AsspDataObj`, `write.AsspDataObj` - SSFF I/O
- `av_to_asspDataObj`, `av_to_python_audio` - Audio format conversion
- `read_avaudio`, `as_avaudio` - AVAudio object creation
- `process_media_file` - Universal media processing

#### D. Audio Object Classes
- `AsspDataObj` - S3 class (traditional wrassp)
- `AVAudio` - S7 class (modern audio wrapper)

#### E. Unit Conversion
- `hz_to_mel`, `mel_to_hz`, `hz_to_bark`, `bark_to_hz`, `hz_to_erb`, `erb_to_hz`
- `hz_to_semitone`, `semitone_to_hz`
- `db_and_hz_to_phon`, `phon_and_hz_to_db`, `phon_to_sone`, `sone_and_hz_to_db`, `sone_to_phon`

#### F. Installation & Configuration
- `install_covarep`, `install_deepformants`, `install_dysprosody`
- `install_ftrack_tvwlp`, `install_gfmiaif`, `install_sacc`
- `install_swiftf0`, `install_vat`, `install_voice_analysis`, `install_voice_sauce`
- `*_available()`, `*_info()` - Module availability checks

#### G. Plotting & Visualization
- `ggtrack` - Auto-labeled ggplot2 formant/pitch plots
- Track label helpers: `get_track_label`, `get_track_label_expr`

#### H. Data Wrangling
- `addTrack`, `delTrack` - Add/remove SSFF tracks
- `store_slice` - Store time slices
- `differentiate`, `afdiff`, `affilter` - Signal processing
- `cut.AsspDataObj` - Time-domain slicing
- `as.data.frame.AsspDataObj`, `as_tibble.AsspDataObj` - Format conversion

---

## 4. AsspDataObj STRUCTURE

### Definition & Purpose
AsspDataObj is the **core data structure** for all DSP output in superassp. It's compatible with wrassp (R package for ASSP library).

**Class:** S3 (though inherits some S4 attributes)
**Data Storage:** List-based with special attributes
**Output Format:** SSFF (Standard Simple Feature Format) binary files

### Core Attributes
```r
# Audio/sampling metadata
@origFreq           # Original sample rate (Hz)
@samplingRate       # Frame rate for analysis (Hz) - for EQUAL-INTERVAL tracks
@startTime          # Time of first frame (seconds)
@startRecord        # First record index (usually 0 or 1)
@endRecord          # Last record index

# Track metadata
@trackFormats       # Vector of track data types (INT16, REAL32, etc.)
@tracks             # List of track matrices

# File metadata
@filePath          # Source/output file path
@fileFormat        # File format (SSFF, WAV, etc.)

# Data format
@dataFormat        # 1=binary, 2=text
```

### Track Structure
Each track is a **matrix** with dimensions:
- **Rows:** Number of frames (equal-interval samples)
- **Columns:** Number of parameters per frame
  - Single-column: `fo[Hz]`, `RMS[dB]`
  - Multi-column: `Fi[Hz]` (F1, F2, F3...), `LPCi` (LPC1, LPC2...)

### Track Naming Convention
Uses **bracket notation** (Titze 2015 / Nylén 2024 compliant):

```
Pattern: Name[Unit]

Examples:
  fo[Hz]           # Fundamental frequency
  Fi[Hz]           # Formant frequencies (template: expands to F1, F2, F3...)
  Bi[Hz]           # Bandwidths
  RMS[dB]          # Root mean square amplitude
  H1-H2c[dB]       # Harmonic amplitude ratio
  LPCi             # LP coefficients (template)
```

**Template Expansion Rules:**
- Placeholder `i` after uppercase letter indicates multi-column expansion
- Pattern: `[A-Z]i(\\[|$)` (uppercase + i + bracket or end)
- Runtime expansion: `Fi[Hz]` → `F1[Hz]`, `F2[Hz]`, `F3[Hz]`, etc.

### Equal-Interval Signal Processing Pattern

**CRITICAL FEATURE:** All tracks are sampled at **equal time intervals**, NOT at the original signal sample rate.

```
Analysis pattern:
1. Load audio at fs=16000 Hz (e.g., 3.2 seconds = 51,200 samples)
2. Window-based analysis (e.g., 10ms shift = 0.01 second intervals)
3. Extract one value per window → equal-interval track
   - Number of frames: ceil(duration / frame_shift)
   - Frame time: startTime + frame_index * (1/samplingRate)
4. Store in AsspDataObj with:
   - samplingRate = 100 Hz (1/0.01 seconds)
   - 320 frames (3.2s / 0.01s)
   - Track matrix: 320 rows × N columns

Example: F0 with 10ms shift
  Frame 0: time=0.000s, f0=145.2 Hz
  Frame 1: time=0.010s, f0=146.1 Hz
  Frame 2: time=0.020s, f0=145.8 Hz
  ...
```

**Key Metadata:**
```r
# For 10ms shift (common)
samplingRate = 100         # Frames per second
startTime = 0.0            # Usually 0, sometimes window-dependent
endRecord = 319            # For 3.2s audio (320 frames total, 0-indexed)

# Frame time calculation:
frame_time = startTime + frame_index / samplingRate
```

### Data Format Support
```r
trackFormats vector possible values:
  "INT16"       # 16-bit signed integer
  "REAL32"      # 32-bit float
  "REAL64"      # 64-bit double (rare)
  "BYTE"        # 8-bit unsigned
  "CHAR"        # Character string
```

### Example Creation Pattern (from egg_f0.py)
```python
# Equal-interval grid creation
duration = len(audio_array) / sample_rate    # seconds
frame_shift_sec = frame_shift_ms / 1000.0    # e.g., 0.01
n_frames = int(np.ceil(duration / frame_shift_sec))
frame_times = np.arange(n_frames) * frame_shift_sec

# Result: array([0.00, 0.01, 0.02, 0.03, ...])
```

Then in R:
```r
new("AsspDataObj",
    trackFormats = c("REAL32", "REAL32"),
    origFreq = 16000,              # Original audio rate
    samplingRate = 1000/10,        # = 100 Hz for 10ms shift
    startTime = 0.0,
    startRecord = 0,
    endRecord = n_frames - 1,
    filePath = "audio.egg",
    tracks = list(
      EGG_F0 = matrix(f0_values, ncol=1),
      EGG_Oq = matrix(oq_values, ncol=1)
    ))
```

---

## 5. TRK_* FUNCTIONS - DETAILED INVENTORY

### Function Naming Convention
Pattern: `trk_<method_name>`

### Patterns by Implementation Type

#### A. Direct R Implementation Files
Only 3 files contain full R implementations:

1. **trk_creak_union.R** (~400 lines)
   - Creaky voice detection using Union Method
   - Combines AM (antimode) and CD (neural network) methods
   - **Unique Features:**
     - Internal Python module initialization: `.init_creak_python()`
     - Complex AsspDataObj creation: `.create_assp_creak_obj()`
     - Multiple output tracks (AM_creak, CD_creak, union_creak, etc.)
     - File I/O helpers: `.construct_creak_output_path()`, `.write_ssff_creak()`

2. **trk_egg_f0.R** (~380 lines)
   - EGG signal F0 and open quotient analysis
   - Uses egg_python module (imported from external path)
   - **Unique Features:**
     - Dynamic module import: `reticulate::import_from_path()`
     - Equal-interval interpolation at frame level
     - Multi-track output with metadata attributes
     - Comprehensive parameter validation

3. **trk_formants_tvwlp.R** (~415 lines)
   - Time-varying weighted LP formant tracking
   - Complex Python backend with optimization levels
   - **Unique Features:**
     - Module loading: `ensure_ftrack_python()`
     - Helper pipeline: `process_with_av_python()`
     - Numpy array conversion: `numpy_to_assp_dataobj()`
     - Output path construction: `construct_output_path()`

#### B. Wrapper Functions in ssff_*.R Files
Most trk_* functions are thin wrappers defined in `ssff_<backend>_<method>.R` files

**Naming Pattern:**
```
ssff_<backend>_<method>.R  →  trk_<method>() function
  backend: cpp, c, python, pm (Praat/Parselmouth)

Examples:
  ssff_cpp_sptk_dio.R      →  trk_dio()
  ssff_python_crepe.R      →  trk_crepe()
  ssff_c_assp_rmsana.R     →  trk_rmsana()
  ssff_python_pm_ppitch.R  →  trk_pitchp()
```

#### C. Wrapper Function Pattern (Template)
```r
# From ssff_cpp_sptk_dio.R (~120 lines)

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
  
  # 1. Input validation
  if (is.null(listOfFiles) || length(listOfFiles) == 0) {
    cli::cli_abort("No input files specified")
  }
  
  listOfFiles <- normalizePath(path.expand(listOfFiles), mustWork = FALSE)
  files_exist <- file.exists(listOfFiles)
  if (!all(files_exist)) {
    cli::cli_abort(c("!", "x" = "{.file {missing_files}}"))
  }
  
  # 2. Parameter recycling
  n_files <- length(listOfFiles)
  beginTime <- if (length(beginTime) == 1) rep(beginTime, n_files) else beginTime
  endTime <- if (length(endTime) == 1) rep(endTime, n_files) else endTime
  
  # 3. Output directory setup
  makeOutputDirectory(outputDirectory, FALSE, "trk_dio")
  
  # 4. Progress bar
  if (verbose && n_files > 1) {
    cli::cli_progress_bar("Processing files", total = n_files)
  }
  
  # 5. Main loop
  results <- vector("list", n_files)
  for (i in seq_len(n_files)) {
    tryCatch({
      # Load audio: av package (universal format support)
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
      
      # Create output AsspDataObj
      out_obj <- create_f0_asspobj(dio_result, windowShift)
      
      # Write or return
      if (toFile) {
        out_file <- generate_output_path(listOfFiles[i], explicitExt, outputDirectory)
        write.AsspDataObj(out_obj, out_file)
        results[[i]] <- TRUE
      } else {
        results[[i]] <- out_obj
      }
    }, error = function(e) {
      cli::cli_warn("Error processing {basename(listOfFiles[i])}: {e$message}")
      results[[i]] <- if (toFile) FALSE else NULL
    })
    
    if (verbose && n_files > 1) cli::cli_progress_update()
  }
  
  # 6. Return logic
  if (toFile) {
    n_success <- sum(unlist(results), na.rm = TRUE)
    if (verbose) cli::cli_inform("Successfully processed {n_success} of {n_files} files")
    return(invisible(n_success))
  } else {
    results <- results[!sapply(results, is.null)]
    if (length(results) == 1) return(results[[1]]) else return(results)
  }
}

# Function attributes for introspection
attr(trk_dio, "ext") <- "f0"
attr(trk_dio, "tracks") <- c("f0")
attr(trk_dio, "outputType") <- "SSFF"
attr(trk_dio, "nativeFiletypes") <- c("wav", "flac", "mp3", "mp4", "mkv", "avi")
attr(trk_dio, "suggestCaching") <- FALSE
```

#### D. Common Backend Functions
All wrappers delegate to backend-specific functions:

**C++ (Rcpp):**
- `dio_cpp()`, `harvest_cpp()`, `rapt_cpp()`, `reaper_cpp()`, `swipe_cpp()`
- `d4c_cpp()`, `sptk_mfcc_cpp()`
- `estk_pitchmark_cpp()`
- `opensmile_gemaps_cpp()`, `opensmile_extract_cpp()`

**Python (via reticulate):**
- `crepe_pitch()`, `swiftf0_pitch()`, `snack_pitch()`, `yin_pitch()`
- Praat/Parselmouth wrappers (many)

**C (wrassp):**
- `acfana_()`, `forest_()`, `lpcana_()`, `rmsana_()`, `zcrana_()`

---

## 6. PYTHON INTEGRATION

### Architecture
```
R Function (trk_*)
    ↓
Wrapper (ssff_*.R)
    ↓
av audio loading (in-memory)
    ↓
reticulate Python interface
    ↓
inst/python/* modules (NumPy/SciPy/PyTorch)
    ↓
Results → AsspDataObj → SSFF file
```

### Python Module Directory Structure
```
inst/python/
├── egg_analysis/              # EGG F0/Oq analysis
│   ├── __init__.py
│   └── egg_f0.py            # Core analysis function
│
├── ftrack_tvwlp/             # Formant tracking TVWLP
│   ├── ftrack/
│   ├── setup.py
│   └── README.md
│
├── union-creak-detection-method/  # Creaky voice detection
│   ├── __init__.py
│   ├── detect_creak_union_extended()
│   └── ...
│
├── covarep_python/           # COVAREP (50+ files)
├── DeepFormants/             # Deep learning formants
├── dysprosody/               # Dysprosody toolkit
├── voice_analysis_python/    # Voice analysis toolkit (80 files)
├── voicesauce/               # Voice Sauce features
├── gfmiaif/                  # GFMIAIF analysis
├── calc_sbpca/               # PCA components
│
├── snack_pitch.py            # SNACK pitch detection
├── snack_formant.py          # SNACK formant tracking
├── praat_*.py                # Parselmouth wrappers (memory-optimized)
├── praat_sauce_memory.py
├── tremor_analysis.py
│
└── Various utility modules
    ├── yaapt.py
    ├── excite.py
    ├── harvest.py
    └── ...
```

### Equal-Interval Processing in Python
**Pattern from egg_analysis/egg_f0.py:**

```python
def analyze_egg_f0(audio_array, sample_rate, frame_shift_ms=10.0):
    # 1. Cycle-based analysis (irregular timing)
    cycle_times_sec, f0_cycles = egg_analysis(audio_array, sample_rate)
    
    # 2. Create equal-interval grid
    duration = len(audio_array) / sample_rate
    frame_shift_sec = frame_shift_ms / 1000.0
    n_frames = int(np.ceil(duration / frame_shift_sec))
    frame_times = np.arange(n_frames) * frame_shift_sec  # [0.00, 0.01, 0.02, ...]
    
    # 3. Interpolate to equal intervals
    f0_interp = np.zeros(n_frames, dtype=np.float32)
    voicing = np.zeros(n_frames, dtype=np.int16)
    
    for i, frame_time in enumerate(frame_times):
        # Find nearest cycle within ±half frame shift
        nearby = np.where(np.abs(cycle_times_sec - frame_time) <= half_shift)[0]
        if len(nearby) > 0:
            nearest_idx = nearby[np.argmin(np.abs(cycle_times_sec[nearby] - frame_time))]
            f0_interp[i] = f0_cycles[nearest_idx]
            voicing[i] = 1
        else:
            f0_interp[i] = 0.0
            voicing[i] = 0
    
    return {
        'f0': f0_interp,                    # [n_frames] - equal intervals
        'voicing': voicing,                  # [n_frames]
        'times': frame_times,                # [n_frames]
        'raw_f0': f0_cycles,                 # Original cycle-based
        'raw_times': cycle_times_sec,        # Original timing
        'n_cycles': len(f0_cycles)
    }
```

### Audio Loading Pattern
All Python-based trk_* functions use this pattern:

```r
# From ssff_python_crepe.R
# Load audio in-memory (supports all formats via av package)
audio_data <- av::read_audio_bin(
  audio = origSoundFile,
  start_time = if (beginTime > 0) beginTime else NULL,
  end_time = if (endTime > 0) endTime else NULL,
  channels = 1
)

# Get sample rate
fs <- attr(audio_data, "sample_rate")

# Convert to float32 for Python/numpy
audio_float <- as.numeric(audio_data) / 2147483647.0  # INT32_MAX

# Create numpy array
np <- reticulate::import("numpy", convert = FALSE)
audio_array <- np$array(audio_float, dtype = "float32")

# Pass to Python and call analysis function
py$waveform <- audio_array
py$fs <- reticulate::r_to_py(as.integer(fs))
py$other_param <- reticulate::r_to_py(value)

# Execute Python code
reticulate::py_run_string("
import librosa
pitch = librosa.yin(waveform, fmin=fMin, fmax=fMax, hop_length=hop_length, sr=fs)
")

# Extract results
result_data <- py$pitch
```

### Key Python Modules Used
- **numpy** - Numerical arrays
- **scipy** - Signal processing (filters, interpolation)
- **librosa** - Audio DSP (pitch, MFCC, spectrograms)
- **torch, torchaudio** - PyTorch deep learning
- **parselmouth** - Praat interface

---

## 7. DATA STRUCTURE CONVERSION HELPERS

### Audio Format Support via av Package
```
Input formats supported:
  .wav, .mp3, .mp4, .m4a, .flac, .ogg, .aac, .opus, .wma
  + all video formats (extracts audio track)

Function: av::read_audio_bin()
  Returns: Integer vector (32-bit signed, s32le format)
           + attributes: sample_rate, channels

Normalization: [-2147483648, 2147483647] → [-1.0, 1.0]
  float_audio = as.numeric(audio_int) / 2147483648.0
```

### Conversion Pipeline: av → Python → AsspDataObj
```
read_audio_bin() → av_to_python_audio() → reticulate::py_to_r()
↓
process in Python (NumPy arrays)
↓
extract results (float matrices)
↓
create AsspDataObj with equal-interval tracks
↓
write.AsspDataObj() → SSFF binary file
```

### AVAudio S7 Class
Modern object-oriented wrapper with validation:

```r
AVAudio <- S7::new_class(
  "AVAudio",
  properties = list(
    samples = S7::class_integer,     # Audio samples (s32le)
    sample_rate = S7::class_integer, # Hz
    channels = S7::class_integer,    # 1, 2, etc.
    file_path = S7::class_character  # Source file (nullable)
  ),
  validator = function(self) {
    # Validates sample_rate > 0, channels > 0, length consistency
  }
)
```

---

## 8. TRACK NAMING & LABELS

### Three-Layer Naming Strategy

**Layer 1: SSFF/AsspDataObj (Scientific Notation)**
```
fo[Hz]           # Fundamental frequency
Fi[Hz]           # Formant frequencies (template)
H1-H2c[dB]       # Harmonic amplitude ratio
```

**Layer 2: Data Frame (R-Friendly)**
```
fo_Hz            # No brackets, underscores
F1_Hz, F2_Hz, F3_Hz  # Template expanded
H1_H2c_dB        # Clean names
```

**Layer 3: Plotting (Display Labels)**
```
Short:  "fo [Hz]"
Long:   "Frequency of oscillation [Hz]"
```

### Template Expansion Functions
From track_helpers.R:

```r
.has_placeholder(name)              # Detects "Fi[Hz]" pattern
.expand_track_template(template, n) # Fi[Hz], 3 → F1[Hz], F2[Hz], F3[Hz]
.clean_track_names(names)           # fo[Hz] → fo_Hz
```

### Label Generation
```r
get_track_label(track_name)         # fo_Hz → "fo [Hz]"
get_track_label_expr(track_name)    # For plotmath expressions
ggtrack(data, ...)                  # Auto-labeled ggplot2
```

---

## 9. SIGNAL PROCESSING PATTERNS

### Window-Based Analysis (Equal-Interval)
All trk_* functions follow this architecture:

```
1. Load audio at native sample rate (e.g., 16 kHz)
2. Apply analysis window (e.g., 20ms Hamming window)
3. Shift by fixed interval (e.g., 10ms)
4. Extract ONE feature per window
5. Store in AsspDataObj with:
   - samplingRate = 1 / windowShift (in seconds)
   - n_frames = total_duration / windowShift
   - Each track: matrix with n_frames rows
```

**Common Parameters:**
- `windowShift` - Frame shift (ms): 5, 10, 20 typical
- `windowSize` - Analysis window (ms): 20, 25, 40 typical
- `centerTime` - Center analysis window? (logical)

**Example: F0 with 10ms shift**
```
Duration: 3.2 seconds at 16 kHz = 51,200 samples
Window: 20ms Hamming window = 320 samples
Shift: 10ms = 160 samples

Frames:
  Frame 0: samples 0-319, time=0.00s → f0=145.2 Hz
  Frame 1: samples 160-479, time=0.01s → f0=146.1 Hz
  Frame 2: samples 320-639, time=0.02s → f0=145.8 Hz
  ...
  Frame 319: samples 51,040-51,199 (partial), time=3.19s → f0=147.3 Hz

AsspDataObj:
  samplingRate = 100 Hz
  endRecord = 319
  f0 track: 320 × 1 matrix
```

### Multi-Column Track Handling
Formant tracking produces multiple values per frame:

```
F1, F2, F3 at each time point
→ Matrix with 3 columns

Example from trk_forest():
  track name template: Fi[Hz]
  time points: 320 (for 10ms shift)
  columns: 4 (F1, F2, F3, F4)
  
  F1_track = matrix(f1_values, nrow=320, ncol=1)
  F2_track = matrix(f2_values, nrow=320, ncol=1)
  F3_track = matrix(f3_values, nrow=320, ncol=1)
  F4_track = matrix(f4_values, nrow=320, ncol=1)
```

OR combined as single multi-column track (less common):

```
Fi_track = matrix(cbind(F1, F2, F3, F4), nrow=320, ncol=4)
# Converted to 4 tracks in AsspDataObj
```

---

## 10. KEY FILES REFERENCE

### Core Infrastructure
- **assp_dataobj.R** (300+ lines)
  - `as.data.frame.AsspDataObj()` - Data frame conversion
  - `as_tibble.AsspDataObj()` - Tibble conversion
  - `.expand_track_template()`, `.clean_track_names()` - Naming

- **wrassp_AsspDataObj.R** (443 lines)
  - `read.AsspDataObj()` - SSFF file reading
  - `write.AsspDataObj()` - SSFF file writing
  - `is.AsspDataObj()` - Validation
  - `print.AsspDataObj()` - Display

- **s7_avaudio.R** (100+ lines)
  - AVAudio S7 class definition
  - `read_avaudio()` - File loading
  - `as_avaudio()` - Format conversion

- **track_helpers.R**
  - `.has_placeholder()` - Template detection
  - `.expand_track_template()` - Name expansion
  - `.clean_track_names()` - Name normalization
  - Unit assignment helpers

- **track_labels_plotmath.R**
  - `get_track_label()` - Short labels
  - `get_track_label_expr()` - plotmath expressions

### Audio I/O
- **av_helpers.R**
  - `av_to_asspDataObj()` - Universal audio loading
  - Automatic fallback to wrassp for niche formats

- **av_python_helpers.R**
  - `av_to_python_audio()` - NumPy conversion
  - `av_load_for_python()` - Combined loading + conversion

- **prep_recode.R**
  - `prep_recode()` - Audio format conversion & normalization

### Plotting
- **ggtrack.R**
  - `ggtrack()` - Auto-labeled formant/pitch plots
  - Automatic unit assignment
  - Theme customization

### Feature Lists
- **list_*.R** files (30+ files)
  - `lst_ComParE_2016()`, `lst_GeMAPS()` - OpenSmile features
  - `lst_voice_sauce()` - Voice Sauce features
  - `lst_covarep_vq()` - COVAREP voice quality
  - Return: character vector of feature names for extraction

---

## 11. FUNCTION ATTRIBUTE METADATA

All trk_* functions include introspection attributes:

```r
attr(trk_dio, "ext")                  # Output file extension: "f0"
attr(trk_dio, "tracks")               # Track templates: c("f0")
attr(trk_dio, "outputType")           # "SSFF"
attr(trk_dio, "nativeFiletypes")      # Supported formats
attr(trk_dio, "suggestCaching")       # Cache results? (logical)
```

These enable:
- Automatic output file naming
- Track structure discovery
- Format validation
- Performance optimization recommendations

---

## 12. DEPENDENCIES & ECOSYSTEM

### Internal Dependencies
- All functions use same AsspDataObj
- reticulate provides Python bridge
- av provides universal audio loading
- Rcpp/C++ for performance-critical code

### External Python Modules
Specified in various `.py` files and `DESCRIPTION`:
- librosa (audio processing)
- scipy (signal processing)
- torch, torchaudio (deep learning)
- parselmouth (Praat automation)
- Many others (see inst/python/)

### Related R Packages
- **wrassp** - Original ASSP library wrapper
- **emuR** - EMU-SDMS integration (compatible)
- **tuneR** - Audio file I/O
- **av** - ffmpeg wrapper for modern formats

---

## SUMMARY TABLE: TRK_* FUNCTIONS BY BACKEND

| Function | Backend | File | Type | Key Features |
|----------|---------|------|------|--------------|
| trk_dio | C++/SPTK | ssff_cpp_sptk_dio.R | Pitch | WORLD vocoder |
| trk_harvest | C++/SPTK | ssff_cpp_sptk_harvest.R | Pitch | High quality |
| trk_rapt | C++/SPTK | ssff_cpp_sptk_rapt.R | Pitch | Real-time capable |
| trk_swipe | C++/SPTK | ssff_cpp_sptk_swipe.R | Pitch | Sawtooth waveform |
| trk_reaper | C++/SPTK | ssff_cpp_sptk_reaper.R | Pitch | REAPer (external) |
| trk_d4c | C++/SPTK | ssff_cpp_sptk_d4c.R | Spectrum | Aperiodic component |
| trk_sptk_mfcc | C++/SPTK | ssff_cpp_sptk_mfcc.R | MFCC | Mel-frequency |
| trk_acfana | C/ASSP | ssff_c_assp_acfana.R | Autocorr | ACF analysis |
| trk_forest | C/ASSP | ssff_c_assp_forest.R | Formants | ASSP formants |
| trk_rmsana | C/ASSP | ssff_c_assp_rmsana.R | Energy | RMS amplitude |
| trk_zcrana | C/ASSP | ssff_c_assp_zcrana.R | Energy | Zero crossing |
| trk_crepe | Python | ssff_python_crepe.R | Pitch | Neural network CNN |
| trk_yin | Python | ssff_python_yin.R | Pitch | YIN algorithm |
| trk_pyin | Python | ssff_python_pyin.R | Pitch | Probabilistic YIN |
| trk_swiftf0 | Python | ssff_python_swiftf0.R | Pitch | SWIFT method |
| trk_snackp | Python | ssff_python_snack_pitch.R | Pitch | SNACK pitch |
| trk_snackf | Python | ssff_python_snack_formant.R | Formants | SNACK formants |
| trk_deepformants | Python | ssff_python_deepformants.R | Formants | Deep learning |
| trk_mfcc | Python | ssff_python_seenc.R | Spectral | MFCC |
| trk_seenc | Python | ssff_python_seenc.R | Spectral | Spectral entropy |
| trk_formantp | Praat/Pars | ssff_python_pm_pformantb.R | Formants | Parselmouth |
| trk_formantpathp | Praat/Pars | ssff_python_pm_pformantpathb.R | Formants | Path following |
| trk_intensityp | Praat/Pars | ssff_python_pm_pintensity.R | Energy | Intensity |
| trk_pitchp | Praat/Pars | ssff_python_pm_ppitch.R | Pitch | Praat autocorr |
| trk_praat_sauce | Praat/Pars | ssff_python_pm_psauce.R | Voice quality | Voice Sauce |
| trk_spectral_momentsp | Praat/Pars | ssff_python_pm_pspectral_moments.R | Spectral | Spectral moments |
| trk_pitchmark | C/ASSP | ssff_cpp_estk_pitchmark.R | Marking | ESTk pitchmarks |
| trk_reaper_pm | Python | ssff_python_reaper_pm.R | Marking | REAPer marks |
| trk_covarep_iaif | Python | covarep_iaif.R | Voice quality | COVAREP IAIF |
| trk_covarep_srh | Python | covarep_srh.R | Voice quality | COVAREP SRH |
| trk_sacc | Python | ssff_python_sacc.R | Voice quality | SACC toolkit |
| trk_vat_srh | Python | vat_srh.R | Voice quality | VAT SRH |
| trk_gfmiaif | Python | ssff_python_gfmiaif.R | Voice quality | GFMIAIF |
| trk_excite | Python | ssff_python_excite.R | Voice quality | Excitation |
| trk_creak_union | Hybrid | trk_creak_union.R | Voice quality | Union method |
| trk_egg_f0 | Python | trk_egg_f0.R | Pitch | EGG signals |
| trk_formants_tvwlp | Python | trk_formants_tvwlp.R | Formants | TVWLP optimized |
| trk_cepstrum | C/ASSP | ssff_c_assp_cepstrum.R | Spectral | Cepstral |
| trk_cssSpectrum | C/ASSP | ssff_c_assp_cssSpectrum.R | Spectral | CSS |
| trk_dftSpectrum | C/ASSP | ssff_c_assp_dftSpectrum.R | Spectral | DFT |
| trk_lpsSpectrum | C/ASSP | ssff_c_assp_lpsSpectrum.R | Spectral | LP-derived |

---

## BEST PRACTICES FOR NEW TRK_* FUNCTIONS

1. **Input/Output:**
   - Accept `listOfFiles` (vector) or `AVAudio` objects
   - Support `toFile` parameter (default TRUE)
   - Use `outputDirectory` for output placement
   - Return invisible count if toFile=TRUE, AsspDataObj if FALSE

2. **Audio Loading:**
   - Use `av::read_audio_bin()` for universal format support
   - Normalize to float: `as.numeric(x) / 2147483648.0`
   - Get sample rate from audio attributes

3. **Equal-Interval Processing:**
   - Define frame shift (windowShift parameter)
   - Create equal-interval grid: `np.arange(n_frames) * frame_shift_sec`
   - Interpolate results to grid (nearest-neighbor or linear)
   - Store with `samplingRate = 1 / frame_shift_sec`

4. **AsspDataObj Creation:**
   - Use explicit attribute setting or `new("AsspDataObj", ...)`
   - Add track templates with placeholders: `Fi[Hz]`, `LPCi`
   - Include track metadata: format, units, descriptions
   - Set `startTime`, `startRecord`, `endRecord`

5. **Error Handling:**
   - Use `cli::cli_abort()`, `cli::cli_warn()` for messages
   - Validate file existence before processing
   - Try/catch errors per file (don't fail entire batch)
   - Return informative error messages

6. **Performance:**
   - Use `av` package (in-memory) not temp files
   - Leverage Python/C++ for heavy computation
   - Add `verbose` parameter for progress feedback
   - Include progress bars for large batches

