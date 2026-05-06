# Graph Report - .  (2026-05-06)

## Corpus Check
- 77 files · ~150,231 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 273 nodes · 379 edges · 25 communities detected
- Extraction: 96% EXTRACTED · 4% INFERRED · 0% AMBIGUOUS · INFERRED: 15 edges (avg confidence: 0.65)
- Token cost: 0 input · 0 output

## Community Hubs (Navigation)
- [[_COMMUNITY_ASSP C Wrappers & Benchmarks|ASSP C Wrappers & Benchmarks]]
- [[_COMMUNITY_Python Pitch Detection Results|Python Pitch Detection Results]]
- [[_COMMUNITY_Benchmark Results & Visualizations|Benchmark Results & Visualizations]]
- [[_COMMUNITY_Voice Quality Indices (AVQIDSI)|Voice Quality Indices (AVQI/DSI)]]
- [[_COMMUNITY_SACCVAT Voice Analysis Benchmarks|SACC/VAT Voice Analysis Benchmarks]]
- [[_COMMUNITY_CoVaRep F0 Benchmark Suite|CoVaRep F0 Benchmark Suite]]
- [[_COMMUNITY_CoVaRep Benchmark Plot Generation|CoVaRep Benchmark Plot Generation]]
- [[_COMMUNITY_SwiftF0 ONNX Pitch Detector|SwiftF0 ONNX Pitch Detector]]
- [[_COMMUNITY_Structured Array Optimization|Structured Array Optimization]]
- [[_COMMUNITY_VAT Parallelization Benchmarks|VAT Parallelization Benchmarks]]
- [[_COMMUNITY_Voxit Feature Benchmarks|Voxit Feature Benchmarks]]
- [[_COMMUNITY_Parallel Processing Benchmarks|Parallel Processing Benchmarks]]
- [[_COMMUNITY_ParselmouthPython Fixes|Parselmouth/Python Fixes]]
- [[_COMMUNITY_Multi-file Test Runner|Multi-file Test Runner]]
- [[_COMMUNITY_SwiftF0 Dependencies & Changelog|SwiftF0 Dependencies & Changelog]]
- [[_COMMUNITY_AVAudio S7 Class|AVAudio S7 Class]]
- [[_COMMUNITY_GeMAPSeGeMAPS Feature Sets|GeMAPS/eGeMAPS Feature Sets]]
- [[_COMMUNITY_Track Template Placeholders|Track Template Placeholders]]
- [[_COMMUNITY_Comprehensive Benchmark Runner|Comprehensive Benchmark Runner]]
- [[_COMMUNITY_Performance Profiler|Performance Profiler]]
- [[_COMMUNITY_SPTK Pitch Trackers|SPTK Pitch Trackers]]
- [[_COMMUNITY_WORLD DIO Vocoder|WORLD DIO Vocoder]]
- [[_COMMUNITY_SwiftF0 Requirements|SwiftF0 Requirements]]
- [[_COMMUNITY_AVAudio AV Converter|AVAudio AV Converter]]
- [[_COMMUNITY_pladdrr Praat Engine|pladdrr Praat Engine]]

## God Nodes (most connected - your core abstractions)
1. `superassp R package` - 12 edges
2. `Benchmarking Suite README` - 12 edges
3. `PitchResult` - 11 edges
4. `SwiftF0` - 10 edges
5. `AVQI Output Report - Acoustic Voice Quality Index v03.01` - 9 edges
6. `generate_test_data_dicts()` - 8 edges
7. `main()` - 8 edges
8. `Python/Parselmouth Implementation - Fixes Applied` - 8 edges
9. `AVAudio S7 class` - 8 edges
10. `generate_test_data_structured()` - 7 edges

## Surprising Connections (you probably didn't know these)
- `praat_dsi() - Dysphonia Severity Index function` --produces--> `DSI Output Report - Dysphonia Severity Index in Praat v02.01`  [INFERRED]
  inst/benchmarking/README.md → tests/signalfiles/DSI/output/1_2021-12-31.pdf
- `praat_avqi() - Acoustic Voice Quality Index function` --produces--> `AVQI Output Report - Acoustic Voice Quality Index v03.01`  [INFERRED]
  inst/benchmarking/README.md → tests/signalfiles/AVQI/output/1.pdf
- `Parallel Processing Performance - 3.63x speedup on 9 cores (20 files)` --result_of--> `Benchmarking Suite README`  [INFERRED]
  inst/benchmarking/results/benchmark_parallel.png → inst/benchmarking/README.md
- `superassp R package` --exports--> `ucnv_* unit conversion functions (14 total)`  [EXTRACTED]
  tests/TRACK_NAME_DIFFERENCES.md → inst/devdocs/PKGDOWN_FUNCTION_GROUPING.md
- `superassp R package` --exports--> `AVAudio S7 class`  [EXTRACTED]
  tests/TRACK_NAME_DIFFERENCES.md → inst/devdocs/S7_AVAUDIO_IMPLEMENTATION.md

## Communities

### Community 0 - "ASSP C Wrappers & Benchmarks"
Cohesion: 0.12
Nodes (28): acfana() function, AsspDataObj S3 class, ASSP C library, benchmark_suite.R, benchmark_suite_simple.R, JSTF (JSON Track Format), lst_* DSP summary functions (17 total), normalize_track_name() helper function (+20 more)

### Community 1 - "Python Pitch Detection Results"
Cohesion: 0.16
Nodes (17): export_to_csv(), PitchResult, plot_pitch(), Container for pitch detection results containing:     - pitch_hz: Estimated fund, Plot pitch with voicing information, optionally saving and/or showing.      Args, Export pitch detection results to CSV file.      Args:         result: PitchResu, SwiftF0 - A fast and accurate fundamental frequency (F0) detector  SwiftF0 is a, export_to_midi() (+9 more)

### Community 2 - "Benchmark Results & Visualizations"
Cohesion: 0.1
Nodes (22): Formant Analysis Performance Comparison - praat_sauce, praat_formant_burg, superassp::forest, wrassp::forest, benchmark_opensmile_slice_functions.R - OpenSMILE Benchmarks, Parallel Processing Performance - 3.63x speedup on 9 cores (20 files), Pitch Algorithm Microbenchmark Timings (Snack, Praat, pYIN, YIN, CREPE, Harvest, DIO, SWIPE), benchmark_python_memory_improvements.R - Memory Optimization Analysis, benchmark_python_ssff.R - Python SSFF Benchmarks, benchmark_suite.R - Comprehensive Benchmark Script, benchmark_suite_simple.R - Quick Benchmark Script (+14 more)

### Community 3 - "Voice Quality Indices (AVQI/DSI)"
Cohesion: 0.17
Nodes (16): Acoustic Voice Quality Index (AVQI), AVQI Output Report - Acoustic Voice Quality Index v03.01, AVQI Reference Illustrated Report - Acoustic Voice Quality Index v03.01, Smoothed Cepstral Peak Prominence (CPPS), Dysphonia Severity Index (DSI), DSI Output Report - Dysphonia Severity Index in Praat v02.01, Test subject: Fredrik Karlsson (born 1975-12-31), Harmonics-to-Noise Ratio (HNR) (+8 more)

### Community 4 - "SACC/VAT Voice Analysis Benchmarks"
Cohesion: 0.23
Nodes (13): benchmark_dfa(), benchmark_dypsa(), benchmark_full_analysis(), benchmark_rpde(), generate_test_signal(), print_summary(), Benchmark Script for Performance Optimization  Tests the performance improvement, Benchmark DYPSA optimization (+5 more)

### Community 5 - "CoVaRep F0 Benchmark Suite"
Cohesion: 0.24
Nodes (13): benchmark_implementation(), generate_test_signal(), load_test_audio(), plot_results(), F0 Tracking Performance Benchmark  Compares performance between: 1. Original imp, Benchmark a single implementation      Parameters     ----------     func : call, Validate numerical accuracy between implementations      Parameters     --------, Plot benchmark results      Parameters     ----------     times_dict : dict (+5 more)

### Community 6 - "CoVaRep Benchmark Plot Generation"
Cohesion: 0.24
Nodes (13): generate_accuracy_plot(), generate_matlab_comparison_plot(), generate_performance_comparison_plot(), generate_scaling_analysis_plot(), generate_summary_infographic(), main(), Generate Comprehensive Performance Plots  Creates publication-quality visualizat, Generate MATLAB comparison plot (+5 more)

### Community 7 - "SwiftF0 ONNX Pitch Detector"
Cohesion: 0.19
Nodes (8): Run the ONNX model to extract pitch and confidence.          Args:             a, Compute voicing mask based on confidence threshold and frequency limits., Calculate accurate frame timestamps accounting for STFT center positions., Detect pitch from numpy array.          Args:             audio_array: Input aud, SwiftF0 - A fast and accurate fundamental frequency (F0) detector using ONNX mod, Detect pitch from audio file.          Args:             audio_path: Path to aud, Initialize SwiftF0 with the bundled ONNX model.          Args:             confi, SwiftF0

### Community 8 - "Structured Array Optimization"
Cohesion: 0.29
Nodes (12): benchmark_conversion_overhead(), benchmark_scaling(), benchmark_with_without_structured(), generate_test_data_dicts(), generate_test_data_structured(), Compare performance with dict input vs structured array input., Test how performance scales with input size., Generate test data as list of dicts (original format). (+4 more)

### Community 9 - "VAT Parallelization Benchmarks"
Cohesion: 0.32
Nodes (10): benchmark_full_analysis(), benchmark_gne_parallelization(), benchmark_hnr_parallelization(), benchmark_rpde_optimization(), load_test_audio(), main(), Test GNE frame-level parallelization, Test full voice analysis with all optimizations (+2 more)

### Community 10 - "Voxit Feature Benchmarks"
Cohesion: 0.33
Nodes (9): benchmark_function(), benchmark_voxit_features(), generate_test_data(), main(), print_summary(), Print summary of all benchmarks., Generate test data for benchmarking., Benchmark a function with warmup iterations. (+1 more)

### Community 11 - "Parallel Processing Benchmarks"
Cohesion: 0.25
Nodes (9): analyze_scalability(), benchmark_batch_processing(), benchmark_single_file(), Benchmark: Sequential vs Parallel Performance  Compares sequential and parallel, Analyze how speedup scales with number of workers, Benchmark sequential vs parallel on single file, Summarize parallelization opportunities and expected gains, Benchmark batch processing with parallelization (+1 more)

### Community 12 - "Parselmouth/Python Fixes"
Cohesion: 0.27
Nodes (10): Python/Parselmouth Implementation - Fixes Applied, Parselmouth (praat-parselmouth Python library), praat_formant_burg_opt() function, praat_formantpath_burg_opt() function, praat_intensity_opt() function, praat_pitch_opt() function, R/praat_python_optimized.R, praat_spectral_moments_opt() function (+2 more)

### Community 13 - "Multi-file Test Runner"
Cohesion: 0.36
Nodes (7): generate_report(), main(), Test Optimizations with Multiple Audio Files  Tests the optimized implementation, Generate comprehensive test report, Main function to test multiple files, Test both F0 and IAIF on a single audio file      Returns     -------     result, test_audio_file()

### Community 14 - "SwiftF0 Dependencies & Changelog"
Cohesion: 0.22
Nodes (9): NoteSegment dataclass, numpy Python package (>=1.21.0), onnxruntime Python package (>=1.12.0), PitchResult dataclass, SwiftF0 Changelog, SwiftF0 Python class, SwiftF0 pitch detection library, SwiftF0 README (+1 more)

### Community 15 - "AVAudio S7 Class"
Cohesion: 0.29
Nodes (8): as_avaudio() function, AVAudio S7 class, avaudio_to_tempfile() function, prep_recode() function, Rationale: AVAudio S7 class for type safety, memory efficiency, preprocessing pipeline, read_avaudio() function, S7 AVAudio Implementation Summary, R/s7_avaudio.R

### Community 16 - "GeMAPS/eGeMAPS Feature Sets"
Cohesion: 0.38
Nodes (7): eGeMAPS Feature Set Revision History, eGeMAPS feature set (openSMILE), GeMAPS Feature Set Revision History, GeMAPS feature set (openSMILE), lst_eGeMAPS() function, lst_GeMAPS() function, openSMILE C++ library

### Community 17 - "Track Template Placeholders"
Cohesion: 0.47
Nodes (6): .expand_track_template() internal function, .has_placeholder() internal function, Rationale: uniform [A-Z]i pattern (no underscore) for track template placeholders, rfcana() function, Uniform 'i' placeholder strategy for track name templates, Uniform Placeholder Strategy: Consistent 'i' Notation

### Community 18 - "Comprehensive Benchmark Runner"
Cohesion: 0.8
Nodes (3): main(), print_header(), print_section()

### Community 19 - "Performance Profiler"
Cohesion: 0.67
Nodes (1): profile_analysis()

### Community 20 - "SPTK Pitch Trackers"
Cohesion: 0.67
Nodes (3): SPTK C++ library, trk_rapt() - RAPT pitch tracker, trk_swipe() - SWIPE pitch estimator

### Community 21 - "WORLD DIO Vocoder"
Cohesion: 1.0
Nodes (2): trk_dio() - DIO F0 estimator, WORLD vocoder C++ library

### Community 34 - "SwiftF0 Requirements"
Cohesion: 1.0
Nodes (1): SwiftF0 Python Requirements

### Community 35 - "AVAudio AV Converter"
Cohesion: 1.0
Nodes (1): avaudio_to_av() function

### Community 36 - "pladdrr Praat Engine"
Cohesion: 1.0
Nodes (1): pladdrr engine (Praat from R via C++)

## Knowledge Gaps
- **90 isolated node(s):** `Container for pitch detection results containing:     - pitch_hz: Estimated fund`, `SwiftF0 - A fast and accurate fundamental frequency (F0) detector using ONNX mod`, `Initialize SwiftF0 with the bundled ONNX model.          Args:             confi`, `Run the ONNX model to extract pitch and confidence.          Args:             a`, `Compute voicing mask based on confidence threshold and frequency limits.` (+85 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **Thin community `Performance Profiler`** (3 nodes): `profile_performance.py`, `profile_analysis()`, `profile_performance.py`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `WORLD DIO Vocoder`** (2 nodes): `trk_dio() - DIO F0 estimator`, `WORLD vocoder C++ library`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `SwiftF0 Requirements`** (1 nodes): `SwiftF0 Python Requirements`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `AVAudio AV Converter`** (1 nodes): `avaudio_to_av() function`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `pladdrr Praat Engine`** (1 nodes): `pladdrr engine (Praat from R via C++)`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `superassp R package` connect `ASSP C Wrappers & Benchmarks` to `AVAudio S7 Class`?**
  _High betweenness centrality (0.016) - this node is a cross-community bridge._
- **Why does `Benchmarking Suite README` connect `Benchmark Results & Visualizations` to `Voice Quality Indices (AVQI/DSI)`?**
  _High betweenness centrality (0.014) - this node is a cross-community bridge._
- **Why does `Voice Quality Analysis Benchmark Comparison Table` connect `Voice Quality Indices (AVQI/DSI)` to `Benchmark Results & Visualizations`?**
  _High betweenness centrality (0.009) - this node is a cross-community bridge._
- **Are the 7 inferred relationships involving `PitchResult` (e.g. with `SwiftF0 - A fast and accurate fundamental frequency (F0) detector  SwiftF0 is a` and `NoteSegment`) actually correct?**
  _`PitchResult` has 7 INFERRED edges - model-reasoned connections that need verification._
- **What connects `Container for pitch detection results containing:     - pitch_hz: Estimated fund`, `SwiftF0 - A fast and accurate fundamental frequency (F0) detector using ONNX mod`, `Initialize SwiftF0 with the bundled ONNX model.          Args:             confi` to the rest of the system?**
  _90 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `ASSP C Wrappers & Benchmarks` be split into smaller, more focused modules?**
  _Cohesion score 0.12 - nodes in this community are weakly interconnected._
- **Should `Benchmark Results & Visualizations` be split into smaller, more focused modules?**
  _Cohesion score 0.1 - nodes in this community are weakly interconnected._