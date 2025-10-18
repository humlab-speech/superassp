# Snack Library Integration Analysis

**Date**: 2025-10-18  
**Objective**: Assess integration of Snack pitch and formant tracking into superassp

## Snack Library Overview

### Source Code Location
- **Path**: `src/tcl-snack/`
- **Key Files**:
  - `generic/jkGetF0.c` (1742 lines) - F0/pitch tracking
  - `generic/jkPitchCmd.c` (1135 lines) - Pitch command interface
  - `generic/jkFormant.c` (1169 lines) - Formant tracking
  - `generic/jkFormant.h` - Formant structures
  - `generic/jkGetF0.h` - F0 structures

### Algorithms

**Pitch Tracking**:
- Dynamic programming-based F0 tracker
- Autocorrelation and cross-correlation methods
- Multiple candidate tracking with cost optimization
- Parameters: min/max F0, frame step, window duration, voicing cost

**Formant Tracking**:
- LPC-based pole estimation
- Dynamic programming for formant mapping
- Tracks F1-F7 (up to 7 formants)
- By David Talkin (AT&T/Entropic Research)

## Integration Challenges

### 1. Tight Tcl/Tk Coupling

**Problem**: Snack's C code is designed for Tcl/Tk:
```c
int Get_f0(Sound *s, Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[])
```

**Issues**:
- Functions expect `Tcl_Interp*` and `Tcl_Obj*` parameters
- Sound data structure is Tcl-specific (`Sound *s`)
- No standalone C API for direct calling
- Heavy dependency on Tcl object system

### 2. Complex Data Structures

**Sound Structure** (Tcl-specific):
```c
typedef struct Sound {
    Tcl_Obj *cmdPtr;      // Tcl command object
    Tcl_Channel rwchan;   // Tcl I/O channel
    // ... many Tcl-specific fields
} Sound;
```

**Formant Structure**:
```c
typedef struct form_latt {
    short ncand;
    short **cand;
    short *prept;
    double *cumerr;
} FORM;

typedef struct pole_array {
    double rms, rms2, f0, pv, change;
    short npoles;
    double *freq;  // pole frequencies
    double *band;  // pole bandwidths
} POLE;
```

### 3. Memory Management

- Extensive use of Tcl's memory allocation (`Tcl_Alloc`, `Tcl_Free`)
- Manual memory management for DP lattices
- Complex cleanup procedures

### 4. Python Interface Exists But...

**tkSnack** (Python binding):
- Requires full Tcl/Tk installation
- Needs Tkinter
- GUI-oriented design
- Not suitable for batch processing

## Options for Integration

### Option 1: Extract and Refactor C Code ❌ NOT RECOMMENDED

**What's needed**:
1. Extract core algorithm from Tcl wrapper
2. Remove all Tcl_* dependencies
3. Create standalone C functions accepting raw audio
4. Write Rcpp wrappers
5. Handle memory management manually

**Effort**: ~40-60 hours  
**Risk**: High - complex refactoring, potential bugs  
**Maintenance**: Difficult - forked code from original

### Option 2: Use Python tkSnack ❌ NOT PRACTICAL

**Approach**: Call tkSnack via reticulate

**Problems**:
- Requires Tcl/Tk + Tkinter setup
- GUI dependencies in headless environment
- Installation complexity
- Not designed for batch processing

### Option 3: Implement Algorithms in Python ✅ RECOMMENDED

**Approach**: Use modern Python libraries with similar algorithms

**For Pitch**:
- **pyworld** - Includes DIO, Harvest (different from Snack's algorithm)
- **parselmouth** - Praat's pitch tracker (already used in package)
- **SPTK** - Already integrated (rapt, swipe, reaper, dio)

**For Formants**:
- **parselmouth** - Praat's formant tracker (can be added)
- **scipy** - LPC-based formant estimation
- Custom implementation using librosa + scipy

### Option 4: Note that Snack Algorithms Already Available ✅ BEST

**Reality Check**:

1. **Snack's F0 tracker** is similar to:
   - `rapt()` - SPTK C++ implementation (already in package)
   - `crepe()` - More accurate modern approach
   - Praat pitch - Via parselmouth

2. **Snack's formant tracker** is based on:
   - David Talkin's work (public domain)
   - Same LPC + DP approach as Praat
   - `praat_formant_burg()` already in package

## Recommendation

### DO NOT integrate Snack directly

**Reasons**:
1. **Already have equivalent/better algorithms**:
   - Pitch: rapt(), swipe(), reaper(), dio(), kaldi_pitch(), crepe()
   - Formants: praat_formant_burg()

2. **Snack's advantages are NOT in algorithms**:
   - Snack's strength is Tcl/Tk integration and GUI
   - Core algorithms are not superior to what we have

3. **Integration cost vs benefit**:
   - High effort to extract C code
   - Low benefit - no new capabilities

4. **Better alternatives exist**:
   - Praat formants (parselmouth) - more features
   - SPTK pitch (C++) - faster, better tested
   - Modern deep learning (crepe) - more accurate

### What TO DO Instead

If users specifically need Snack-style processing:

#### Option A: Document equivalents in superassp

```r
# Snack pitch → Use rapt() or swipe()
f0 <- rapt("audio.wav", minF = 75, maxF = 500)

# Snack formants → Use praat_formant_burg()
formants <- praat_formant_burg("audio.wav", numFormants = 5)
```

#### Option B: Add more Praat formant features

Extend `praat_formant_burg()` with:
- Formant bandwidth tracking
- Formant path optimization
- Additional formant analysis methods

#### Option C: Create comprehensive formant function

Implement new `formants()` function that:
- Uses LPC analysis (via scipy/librosa in Python)
- Provides F1-F5 tracking
- Returns multi-track SSFF output
- Follows same pattern as mfcc()

## Current Package Capabilities

### Pitch Tracking (Excellent Coverage)
| Function | Algorithm | Speed | Accuracy | Implementation |
|----------|-----------|-------|----------|----------------|
| rapt() | NCCF + DP | ★★★★★ | ★★★★☆ | SPTK C++ |
| swipe() | SWIPE | ★★★★★ | ★★★★☆ | SPTK C++ |
| reaper() | REAPER | ★★★★★ | ★★★★☆ | SPTK C++ |
| dio() | DIO | ★★★★☆ | ★★★☆☆ | Python (pyworld) |
| kaldi_pitch() | NCCF | ★★★★☆ | ★★★☆☆ | Python (torch) |
| crepe() | Deep Learning | ★★☆☆☆ | ★★★★★ | Python (keras) |

### Formant Tracking (Limited)
| Function | Algorithm | Implementation |
|----------|-----------|----------------|
| praat_formant_burg() | Burg LPC | Python (parselmouth) |

### Gap: Enhanced Formant Analysis
- Multiple formant tracking algorithms
- Formant bandwidth estimation
- Formant trajectory smoothing
- Dynamic formant tracking

## Proposed: Enhanced Formant Function

Instead of Snack integration, create modern formant analysis:

```r
# New function
formants <- formant_lpc(
  "audio.wav",
  numFormants = 5,
  maxFormantHz = 5500,
  windowShift = 5,
  preEmphasis = 0.97,
  toFile = TRUE
)

# Returns SSFF with tracks:
# - fm_1, fm_2, fm_3, fm_4, fm_5 (formant frequencies)
# - bw_1, bw_2, bw_3, bw_4, bw_5 (bandwidths)
```

**Implementation**: External Python script using:
- `scipy.signal.lfilter()` for pre-emphasis
- `librosa` or `scipy` for LPC analysis
- Root finding for formant extraction
- Following mfcc()/kaldi_pitch() pattern

## Conclusion

**DO NOT** spend time extracting Snack C code. The package already has superior alternatives for pitch tracking, and formant tracking can be enhanced more efficiently using modern Python libraries following the proven external script pattern.

**RECOMMEND**:
1. Document Snack equivalents in package documentation
2. Enhance `praat_formant_burg()` if needed
3. Consider adding `formant_lpc()` function using Python/scipy
4. Focus development on areas where package has gaps, not duplicating existing capabilities

---

**Decision**: Focus on enhancing formant analysis with modern approaches rather than integrating legacy Tcl/Tk code.
