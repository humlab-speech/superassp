# Praat to Parselmouth Migration Plan

## Current State Assessment

### Functions Already Migrated to Parselmouth

#### SSFF Functions (in R/praat_python_optimized.R)
- ✅ `praat_formant_burg_opt()` → uses `inst/python/praat_formant_burg.py`
- ✅ `praat_formantpath_burg_opt()` → uses `inst/python/praat_formantpath_burg.py`
- ✅ `praat_intensity_opt()` → uses `inst/python/praat_intensity.py`
- ✅ `praat_spectral_moments_opt()` → uses `inst/python/praat_spectral_moments.py`
- ✅ `praat_pitch_opt()` → uses `inst/python/praat_pitch.py`

#### Slice Functions (in R/praat_slicefunctions.R)
- ✅ `praat_voice_report_opt()` → uses `inst/python/praat_voice_report_memory.py`

### Functions Still Using inst/praat Scripts

#### SSFF Functions (in R/praat_ssff.R)
- ❌ `praat_formant_burg()` → uses `inst/praat/formant_burg.praat`
- ❌ `praat_formantpath_burg()` → uses `inst/praat/formantpath_burg.praat`
- ❌ `praat_praatsauce()` → uses `inst/praat/praatsauce.praat`
- ❌ `praat_intensity()` → uses `inst/praat/intensity.praat`
- ❌ `praat_spectral_moments()` → uses `inst/praat/praat_spectral_moments.praat`
- ❌ `praat_pitch()` → uses `inst/praat/praat_pitch.praat`

#### Slice Functions (in R/praat_slicefunctions.R)
- ❌ `praat_avqi()` → uses `inst/praat/AVQI301.praat`
  - **Python code exists**: `inst/python/avqi_3.01.py` but no R wrapper
- ❌ `praat_voice_report()` → uses `inst/praat/praat_voice_report.praat`
  - **Parselmouth version exists**: `praat_voice_report_opt()`
- ❌ `praat_dsi()` → uses `inst/praat/DSI201.praat`
  - **No Python version exists**
- ❌ `praat_voice_tremor()` → uses `inst/praat/tremor3.05/console_tremor305.praat`
  - **No Python version exists**

## Migration Strategy

### Phase 1: Create Missing Parselmouth Implementations

#### High Priority (Commonly Used)
1. **praat_avqi**: Create R wrapper `praat_avqi_opt()` using existing `avqi_3.01.py`
2. **praat_praatsauce**: Create Python implementation and R wrapper

#### Medium Priority
3. **praat_dsi**: Create Python implementation `dsi.py` and R wrapper `praat_dsi_opt()`
4. **praat_voice_tremor**: Create Python implementation `voice_tremor.py` and R wrapper `praat_voice_tremor_opt()`

### Phase 2: Update Original Functions to Use Parselmouth

Two options:
- **Option A (Conservative)**: Keep both versions, make original functions call `_opt` versions internally if Parselmouth available
- **Option B (Complete Migration)**: Replace original functions entirely with Parselmouth implementations

### Phase 3: Remove inst/praat Dependency

1. Update all functions to use Parselmouth versions
2. Remove all `system.file(..., "praat", ...)` calls
3. Delete `inst/praat/` folder
4. Update package dependencies (remove Praat requirement, add parselmouth)
5. Update documentation to remove Praat references

### Phase 4: Testing & Validation

1. Run all existing tests
2. Create equivalence tests comparing old vs new implementations
3. Benchmark performance improvements
4. Document any behavioral differences

## Implementation Notes

### praatsauce Migration

The `praatsauce.praat` script needs review to understand what it does:
- Formant analysis
- Voice source analysis
- Spectral tilt measurements
- Custom Praat plugin integration

This may require significant Python implementation work.

### DSI Migration

DSI computation requires:
- Maximum phonation time
- Softest intensity measurement
- Highest F0 measurement
- Jitter (ppq5) measurement
- Formula: DSI = 0.13 * MPT + 0.0053 * F0max - 0.26 * I(low) - 1.18 * Jitter(ppq5) + 12.4

Can be implemented using Parselmouth primitives.

### Voice Tremor Migration

Complex analysis requiring:
- Frequency tremor analysis
- Amplitude tremor analysis
- Cyclicality measurements
- HNR calculations

Requires careful porting of tremor3.05 algorithms to Python.

## Benchmarking Infrastructure

### Current State
- `benchmarking/` folder exists with `benchmark_python_ssff.R`
- `vignettes/` folder exists with `benchmark_report.qmd`
- Test benchmarks in `tests/`:
  - `benchmark_opensmile_slice_functions.R`
  - `benchmark_python_memory_improvements.R`
  - `benchmark_suite.R`
  - `benchmark_suite_simple.R`

### Needed Benchmarking Comparisons

1. **Formant Estimation Procedures**:
   - Burg method (Praat) vs optimized
   - FormantPath vs standard
   - Comparison with ASSP formants

2. **Pitch/F0 Estimation Procedures**:
   - Praat autocorrelation
   - Praat cross-correlation
   - RAPT
   - REAPER
   - SWIPE
   - Dio (WORLD)
   - Harvest (WORLD)
   - SPICE
   - Comparison metrics: accuracy, speed, robustness

3. **Voice Quality Analysis**:
   - Voice report metrics
   - AVQI computation time
   - DSI computation time

### Vignette Structure

Create comprehensive benchmarking vignette:
- `vignettes/benchmarking_guide.Rmd` - Main benchmarking documentation
- `vignettes/formant_comparison.Rmd` - Formant estimator comparison
- `vignettes/pitch_comparison.Rmd` - Pitch estimator comparison
- `vignettes/performance_optimization.Rmd` - Memory-based vs file-based performance

## Risks & Considerations

### Breaking Changes
- Removing Praat dependency may break user workflows
- Need clear migration guide
- Consider deprecation warnings before removal

### Behavioral Differences
- Parselmouth may have subtle numerical differences from external Praat
- Need to document acceptable tolerance levels
- Some edge cases may behave differently

### Dependencies
- Adds Python dependency (already present)
- Requires parselmouth package installation
- May need conda environment setup instructions

## Recommended Approach

Given the scope, I recommend:

1. **Immediate**: Create wrappers for existing Python code (avqi_3.01.py)
2. **Short-term**: Implement DSI in Parselmouth (simpler than tremor)
3. **Medium-term**: Implement praatsauce if widely used
4. **Long-term**: Voice tremor (most complex, possibly least used)
5. **Parallel**: Organize benchmarking and create vignettes
6. **Final**: Remove inst/praat only after thorough testing

This allows incremental migration with reduced risk.
