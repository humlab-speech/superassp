# An extended wrassp package

The idea is to make a package that has all the functionality of wrassp, and extend it with analyses made avaiable in Praat or MATLAB. The added functions should behave in a wrassp-like manner, and thereby be callable in a similar way in the `emuR` framwork.

The `praat_formant_burg` provides an illustration of how a Praat script that extracts formant values may be wrapped inside of an R function and produce a SSFF formant track file. 

## Details
By loading this package, you also get all the functions exported by the `wrassp` package into your namespace. This is achieved by the `superassp` package being *Depending*  the `wrassp` package (rather than *Importing*, which is usually the preferred way of creating depmendencies between R packages).

## Installation

The package requires the Praat program to be installed in the user's PATH (or in '/Applications' on Mac OS).

Then simply install the package using
```r
install.packages("devtools") # If not installed already
devtools::install_github("humlab-speech/superassp",dependencies = "Imports")
```

## Performance Benchmarks

The following benchmarks were run on the current version of `superassp` using a 4-second audio file from the package's sample data. Each benchmark shows the distribution of execution times across 100 runs using violin plots.

### Formant Analysis

Multiple formant tracking methods are available with different speed/feature tradeoffs:

![Formant Analysis Benchmark](benchmark_formant.png)

**Performance comparison** (4-second audio file, median of 100 runs):
- **superassp::forest**: ~146 ms - Fastest, optimized with av-based media loading
- **wrassp::forest**: ~166 ms - Fast, native WAV files only
- **praat_formant_burg**: ~893 ms - Slower, Praat algorithm via Parselmouth
- **praat_sauce**: ~947 ms - Slowest, but computes many additional voice quality measures

The `superassp::forest` function provides the best performance while supporting any media format (including video files) via the `av` package. The Praat-based functions offer additional features but with higher computational cost due to Python/Parselmouth overhead.

### Pitch Tracking Algorithms

`superassp` provides a comprehensive suite of pitch tracking algorithms with varying speed/accuracy tradeoffs:

![Pitch Tracking Benchmark](benchmark_pitch.png)

#### Algorithm Categories

**Native C/C++ Implementations** (Fastest):
- **KSV F0** (`ksvfo`): ~18 ms - Fastest, autocorrelation-based from ASSP library
- **MHS Pitch** (`mhspitch`): ~52 ms - Fast, modified harmonic sieve (cepstrum) from ASSP library
- **ESTK PDA** (`estk_pda_cpp`): NEW! Super-resolution pitch detection using cross-correlation
  - C++ implementation from Edinburgh Speech Tools
  - In-memory processing, accepts `AsspDataObj` directly
  - Sub-sample accuracy with optional peak tracking

**SPTK C++ Implementations** (Fast, native performance):
- **RAPT C++** (`rapt_cpp`): NEW! Native C++ RAPT implementation from SPTK
  - Direct C++ implementation, no Python overhead
  - In-memory processing, accepts `AsspDataObj` directly
  - ~2-3x faster than Python version (~40-60 ms estimated)
- **SWIPE C++** (`swipe_cpp`): NEW! Native C++ SWIPE implementation from SPTK
  - Sawtooth waveform inspired pitch estimator
  - In-memory processing, no file I/O required
  - ~2-3x faster than Python version (~35-50 ms estimated)
- **REAPER C++** (`reaper_cpp`): NEW! Native C++ REAPER implementation from SPTK
  - Robust epoch and pitch estimator with glottal closure detection
  - Returns F0, epochs, and polarity information
  - ~2-3x faster than Python version (~150-200 ms estimated)
- **DIO C++** (`dio_cpp`): NEW! DIO algorithm from WORLD vocoder (SPTK)
  - High-quality pitch extraction for speech synthesis
  - In-memory processing, accepts `AsspDataObj` directly

**Python/SPTK Implementations** (Superseded - use C++ versions instead):
- ~~**SWIPE** (`swipe`)~~: Superseded by `swipe_cpp` - Python wrapper is 2-3x slower
- ~~**RAPT** (`rapt`)~~: Superseded by `rapt_cpp` - Python wrapper is 2-3x slower
- ~~**REAPER** (`reaper`)~~: Superseded by `reaper_cpp` - Python wrapper is 2-3x slower
- **Note**: These Python-based functions are no longer exported. Use the C++ versions (`rapt_cpp`, `swipe_cpp`, `reaper_cpp`) for better performance.

**Praat-based Implementation** (Flexible, requires Parselmouth):
- **Praat Pitch** (`praat_pitch`, `praat_pitch_opt`): Uses Praat's autocorrelation method
  - Optimized version available via Parselmouth
  - Compatible with Praat scripts and workflows
  - Extensive parameter control

**Additional ESTK Algorithms**:
- **ESTK Pitchmark** (`estk_pitchmark_cpp`): Glottal closure instant detection
  - Designed for laryngograph (EGG) signals
  - Optional F0 conversion from pitchmark intervals
  - Configurable filtering and period constraints

#### Performance Characteristics

**Speed vs. Accuracy Tradeoff**:
- **Fastest** (< 60 ms): KSV, MHS, SPTK C++ (RAPT, SWIPE) - Best for real-time or batch processing
- **Moderate** (60-200 ms): ESTK PDA, REAPER C++, DIO C++ - Good balance of speed and accuracy
- ~~**Slow** (100-150 ms): Python SWIPE, Python RAPT - Superseded by C++ versions~~
- ~~**Slowest** (> 400 ms): Python REAPER - Superseded by C++ version~~

**Implementation Details**:
- **C/C++ methods** (KSV, MHS, ESTK, SPTK C++): No dependencies, fastest, in-memory processing
- **SPTK C++ methods** (rapt_cpp, swipe_cpp, reaper_cpp, dio_cpp): Native C++ implementations, 2-3x faster than deprecated Python equivalents
- ~~**Python methods** (nonopt_swipe, nonopt_rapt, nonopt_reaper): **DEPRECATED** - No longer exported, use C++ versions instead~~
- **Praat methods**: Require Parselmouth, flexible but slower

All algorithms support:
- Configurable F0 range (minF/maxF)
- Frame shift/window control (windowShift)
- Batch processing with automatic parallelization
- Output to SSFF track format or in-memory objects

### Parallel Processing Performance

As of version 0.5.2, `superassp` automatically uses parallel processing for batch operations:

![Parallel Processing Benchmark](benchmark_parallel.png)

**Speedup: ~3.6x on 9 cores** when processing 20 files (80 seconds of audio total, median of 100 runs).

The violin plots clearly show the performance advantage of parallel processing, with consistently lower execution times and reduced variability.

Parallel processing is:
- **Automatically enabled** for batches (2+ files)
- **Automatically disabled** for single files
- **Platform-aware**: Uses fork-based parallelism on Unix/Mac, socket clusters on Windows
- **Thread-safe**: All DSP functions use independent memory structures

### Running the Benchmarks

You can reproduce these benchmarks by running:

```r
# Install required packages
install.packages(c("microbenchmark", "ggplot2"))

# Run benchmark script (from package root)
source(system.file("benchmarks", "run_benchmarks.R", package = "superassp"))
```

Or manually:

```r
library(superassp)
library(microbenchmark)

# Get sample file
test_file <- system.file("samples", "sustained", "a32b.wav", package = "superassp")

# Benchmark formant analysis methods
microbenchmark(
  "wrassp::forest" = wrassp::forest(test_file, toFile = FALSE),
  "superassp::forest" = forest(test_file, toFile = FALSE, verbose = FALSE),
  "praat_formant_burg" = praat_formant_burg(test_file, toFile = FALSE),
  "praat_sauce" = praat_sauce(test_file, toFile = FALSE),
  times = 100
)

# Benchmark pitch tracking
# First load audio for C++ in-memory algorithms
audio_obj <- av_to_asspDataObj(test_file)

# C/C++ methods (fastest)
microbenchmark(
  "KSV" = ksvfo(test_file, toFile = FALSE, verbose = FALSE),
  "MHS" = mhspitch(test_file, toFile = FALSE, verbose = FALSE),
  "ESTK_PDA" = estk_pda_cpp(audio_obj, minF = 60, maxF = 400),
  times = 100
)

# SPTK C++ methods (fast, native implementations)
microbenchmark(
  "RAPT_CPP" = rapt_cpp(audio_obj, minF = 60, maxF = 400, windowShift = 10),
  "SWIPE_CPP" = swipe_cpp(audio_obj, minF = 60, maxF = 400, windowShift = 10),
  "REAPER_CPP" = reaper_cpp(audio_obj, minF = 60, maxF = 400, windowShift = 10),
  "DIO_CPP" = dio_cpp(audio_obj, minF = 60, maxF = 400, windowShift = 10),
  times = 100
)

# Note: Python/SPTK methods (swipe, rapt, reaper) have been superseded
# by faster C++ implementations (swipe_cpp, rapt_cpp, reaper_cpp) shown above

# Praat method (requires Parselmouth)
microbenchmark(
  "Praat_Pitch" = praat_pitch_opt(test_file, toFile = FALSE,
                                   pitch_floor = 60, pitch_ceiling = 400),
  times = 50
)

# Benchmark parallel processing
test_files <- rep(test_file, 20)
microbenchmark(
  "Sequential" = lapply(test_files, function(f) rmsana(f, toFile = FALSE, verbose = FALSE)),
  "Parallel" = rmsana(test_files, toFile = FALSE, verbose = FALSE),
  times = 100
)
```

**Note**: These timings represent the performance in the current version of `superassp` and include overhead from media file loading via the `av` package. The relative performance between algorithms is indicative, but absolute times will vary based on audio file properties and system specifications.

# Steps to implement a new Praat function

1. Indentify what the output of the function will be
    * A signal track (or tracks) that follows the original sound wave
    * A value (or a limited list of values) that summarises the acoustic properties of a wav file, and can therefore not sensibly be shown alongside the sound wave.
2. Implement the core analysis in a Praat script file, and place it in `inst/praat`.
    * In the case where track(s) that follow the sound wave file are returned, the Praat function should write the output to a CSV table file and return the name of that table. The Praat script should also take the desired output table file name (including full path) as an argument. Please refer to `praat/formant_burg.praat` for some example code that computes formants and bandwidths for them for a (possibly windowed) sound file and writes them to a table.
3. Make a copy of the suitable template function, rename it (please keep the praat_ prefix for clarity) and make modifications to the code to suit the new track computed by Praat. You will need to think about what the tracks should be called in the SSFF file and document your choice.
    * For a function that computes a sound wave following signal track (or tracks), use the code of `praat_formant_burg` as a template. Please refer to a suitable function in wrassp for inspiration on what to call sets of tracks. (The `praat_formant_burg` outputs and "fm" and "bw" set, for formant frequencies and formant bandwidths respectivelly)
    * For single value (or list of values) output, there is currently no template function implemented, but please note that the `tjm.praat::wrap_praat_script()`, which `cs_wrap_praat_script` is a revised version of, has an option to return the "Info window" of Praat, which opens up lots of possibilities.
4. There are many moving parts to this whole package, so make sure to contruct a test file and a test suit for the new function to make sure that it works. 
