# superassp Strategic Recommendations

**Date**: 2025-11-08  
**Version**: 0.10.0  
**Audience**: Package maintainers and core contributors

---

## Executive Summary

The superassp package has reached a critical milestone with v0.10.0. This document provides strategic recommendations for continuing package development, focusing on areas that will maximize impact and value for the speech processing community.

---

## Priority Matrix

### Critical (Do First) - High Impact, Urgent

1. **Complete Function Modernization** (46% remaining)
   - Impact: Consistency, performance, universal media support
   - Effort: Medium (11 Python functions)
   - Timeline: 2-3 weeks
   - Dependencies: av package (already integrated)

2. **Formalize Testing Infrastructure**
   - Impact: Reliability, confidence in changes
   - Effort: Medium (write tests for uncovered functions)
   - Timeline: 2 weeks
   - Target: 80%+ coverage

### Important (Plan Soon) - High Impact, Not Urgent

3. **Performance Benchmarking Suite**
   - Impact: Attract users, validate optimizations
   - Effort: Low (infrastructure exists)
   - Timeline: 1 week
   - Deliverable: Published benchmark results

4. **Dependency Optimization**
   - Impact: Faster installation, fewer conflicts
   - Effort: Medium (audit 38 packages)
   - Timeline: 2 weeks
   - Potential: Move 10+ to Suggests

### Strategic (Long-Term) - High Value, Low Urgency

5. **Package Architecture Evolution**
   - Impact: Modularity, maintainability, discoverability
   - Effort: High (breaking changes)
   - Timeline: 6 months
   - Target: v1.0.0

6. **Community Building**
   - Impact: Adoption, contributions, citations
   - Effort: Medium (ongoing)
   - Timeline: Continuous
   - Channels: Website, examples, tutorials

---

## Detailed Recommendations

### 1. Complete Function Modernization ⚡ CRITICAL

**Current Status**: 54% complete (29 of 54 functions)

**Remaining Work**:
- 11 Python functions using librosa.load
- 10 Parselmouth functions needing av integration
- 2 PyTorch functions (consider deprecation)

**Implementation Strategy**:

```r
# Week 1: High-priority pitch functions
- trk_crepe()     # Deep learning pitch
- trk_pyin()      # Probabilistic YIN
- trk_yin()       # Classic YIN

# Week 2: Parselmouth batch (6 functions)
- Follow av_load_for_parselmouth pattern
- Implement Pattern 2 (R creates Sound object)
- Test with multiple formats

# Week 3: Remaining functions + documentation
- Complete librosa migrations
- Update CLAUDE.md
- Comprehensive testing
```

**Success Metrics**:
- 100% of functions accept any media format
- All functions use av package (no librosa dependencies)
- S7 AVAudio dispatch for all trk_* and lst_* functions
- Documentation updated with migration examples

**ROI**: High - Universal media support is a major competitive advantage

---

### 2. Testing Infrastructure ⚡ CRITICAL

**Current Status**: Grade A (417 test cases), but gaps remain

**Target Coverage**:
- Core DSP functions: 100%
- Python integrations: 90% (with availability checks)
- Helper functions: 80%
- Overall package: 85%+

**Testing Strategy**:

```r
# Priority 1: Untested DSP functions
- Identify functions with <50% coverage
- Write basic functionality tests
- Add parameter variation tests
- Test error handling

# Priority 2: Integration tests
- Multi-file batch processing
- Format conversion edge cases
- Parallel processing validation
- Memory leak detection

# Priority 3: Regression tests
- Known bug fixes (RAPT, DIO, etc.)
- Performance benchmarks
- Output format validation
```

**Infrastructure Improvements**:
```r
# Add test utilities in tests/testthat/helper-testing.R
create_test_audio <- function(duration = 1.0, f0 = 220, sr = 16000) {
  # Generate synthetic audio for testing
}

expect_valid_ssff <- function(result, tracks) {
  # Validate SSFF structure
}

expect_avaudio <- function(result) {
  # Validate AVAudio object
}
```

**Success Metrics**:
- covr::package_coverage() > 85%
- All exported functions have at least basic tests
- CI/CD integration (GitHub Actions)
- Automated test reporting

**ROI**: Very High - Prevents regressions, enables confident refactoring

---

### 3. Performance Benchmarking Suite 📊 IMPORTANT

**Objective**: Formalize and publish comprehensive benchmarks

**Current State**: 
- Ad-hoc benchmarking exists
- No published results
- No automated tracking

**Proposed Structure**:

```r
# inst/benchmarks/benchmark_suite.R
benchmark_suite <- function() {
  results <- list(
    pitch_tracking = benchmark_pitch(),
    formant_tracking = benchmark_formants(),
    spectral_analysis = benchmark_spectral(),
    voice_quality = benchmark_vq(),
    feature_extraction = benchmark_features()
  )
  
  # Generate markdown report
  generate_benchmark_report(results)
  
  # Update README.md with results
  update_readme_benchmarks(results)
}

benchmark_pitch <- function() {
  # Test all pitch trackers on common audio
  # Measure: speed, accuracy, memory
  # Compare: C++ vs Python vs Praat
}
```

**Benchmark Metrics**:
1. **Speed**: Real-time factor (processing time / audio duration)
2. **Accuracy**: Compared to ground truth (when available)
3. **Memory**: Peak memory usage
4. **Scalability**: Performance with varying file sizes

**Deliverables**:
- `inst/benchmarks/` directory with scripts
- `docs/reference/BENCHMARKS.md` with results
- README.md updated with performance highlights
- Automated benchmark runs on releases

**Success Metrics**:
- All major functions benchmarked
- Results published in README
- Community can reproduce benchmarks
- Performance regressions detected automatically

**ROI**: Medium-High - Attracts performance-conscious users, validates optimizations

---

### 4. Dependency Optimization 📦 IMPORTANT

**Current State**: 38 imported packages (heavy dependency tree)

**Audit Strategy**:

```r
# Analyze actual usage
library(pkgnet)
report <- CreatePackageReport("superassp")

# For each dependency:
# 1. Count function calls
# 2. Assess if Suggests vs Imports
# 3. Identify potential replacements
```

**Candidates for Suggests** (Preliminary):
- `R.utils` - Used minimally
- `tidyselect` - May be redundant after av migration
- `pbapply`, `pbmcapply` - Already in Suggests ✅
- `lifecycle` - Only for deprecation (could inline)

**Keep in Imports** (Essential):
- `av` - Core media loading ✅
- `Rcpp` - C++ integration ✅
- `reticulate` - Python integration ✅
- `parallel` - Core functionality ✅
- `S7` - AVAudio dispatch ✅

**Potential Eliminations**:
- Can we replace tidyr with base R? (analyze usage)
- Can we replace assertthat with checkmate? (faster, fewer deps)
- Can we replace logger with custom logging? (simpler)

**Success Metrics**:
- Imports reduced from 38 to <30
- Faster installation time
- Fewer CRAN dependency conflicts
- Package size reduced

**ROI**: Medium - Improves installation experience, reduces maintenance burden

---

### 5. Package Architecture Evolution 🏗️ STRATEGIC

**Vision**: Modular architecture for v1.0.0

**Proposed Split**:

```
superassp-core (lightweight)
├── C/C++ implementations (ASSP, SPTK, ESTK)
├── Basic R wrappers
├── Core utilities (av loading, SSFF I/O)
└── No Python dependencies

superassp-ml (deep learning)
├── Neural network pitch trackers (Swift-F0, CREPE, TANDEM)
├── Deep learning formants (DeepFormants)
├── Python required
└── Depends: superassp-core

superassp-features (high-level)
├── OpenSMILE features (GeMAPS, eGeMAPS)
├── Voice quality measures (VoiceSauce, VAT, COVAREP)
├── Prosodic features (Dysprosody)
└── Depends: superassp-core

superassp (meta-package)
├── Imports: all three
└── Convenience wrapper
```

**Benefits**:
- Users can install only what they need
- Lighter core package (<10 MB vs current 50+ MB)
- Clearer dependency management
- Better CRAN compliance

**Challenges**:
- Breaking change for existing users
- More complex release management
- Need migration guide

**Timeline**: v1.0.0 (6+ months)

**ROI**: High (long-term) - Better modularity, wider adoption

---

### 6. Community Building 👥 STRATEGIC

**Current State**:
- Excellent technical package
- Limited visibility
- Few external contributors
- No formal examples/tutorials

**Community Strategy**:

**A. Documentation Website**
```r
# Use pkgdown
usethis::use_pkgdown()

# Customize _pkgdown.yml
navbar:
  - title: Functions
  - title: Articles
    - Pitch Tracking Guide
    - Formant Analysis Tutorial
    - Voice Quality Assessment
    - Batch Processing Workflows
  - title: Benchmarks
```

**B. Example Gallery**
```r
# Create vignettes/
vignettes/
├── pitch-tracking-comparison.Rmd   # Compare 10+ pitch trackers
├── formant-analysis-workflow.Rmd   # End-to-end formant analysis
├── voice-quality-assessment.Rmd    # VQ measures explained
├── batch-processing-guide.Rmd      # Processing large datasets
└── custom-pipelines.Rmd            # Building custom DSP pipelines
```

**C. Publication Strategy**
- Submit to JOSS (Journal of Open Source Software)
- Write methodology paper for journal
- Present at conferences (Interspeech, ICPhS)
- Create citation.cff file

**D. Engagement Channels**
- GitHub Discussions (Q&A, showcases)
- Twitter/Mastodon (updates, tips)
- Speech processing forums (promotion)
- University workshops (training)

**Success Metrics**:
- GitHub stars: 50+ (currently ~10?)
- Citations: 10+ per year
- Contributors: 5+ external
- pkgdown website deployed
- 3+ vignettes published

**ROI**: High - Increases adoption, contributions, academic recognition

---

## Implementation Roadmap

### Quarter 1 (Months 1-3)

**Month 1: v0.10.1 - Function Modernization**
- Week 1-2: Migrate 11 Python functions (librosa → av)
- Week 3-4: Migrate 10 Parselmouth functions
- Release: v0.10.1 (100% av integration)

**Month 2: v0.10.2 - Testing & Quality**
- Week 1-2: Write tests for uncovered functions
- Week 3: CI/CD setup (GitHub Actions)
- Week 4: Coverage reporting, bug fixes
- Release: v0.10.2 (85%+ test coverage)

**Month 3: v0.10.3 - Performance & Docs**
- Week 1-2: Formalize benchmark suite
- Week 3: Publish benchmark results
- Week 4: Dependency audit, optimization
- Release: v0.10.3 (benchmarked, optimized)

### Quarter 2 (Months 4-6)

**Month 4: v0.11.0 - Community Features**
- Week 1-2: pkgdown website setup
- Week 3-4: Write 3 vignettes
- Release: v0.11.0 (documentation website)

**Month 5: v0.11.1 - Refinement**
- Week 1-2: User feedback integration
- Week 3-4: Performance improvements
- Release: v0.11.1 (refined, stable)

**Month 6: v0.12.0 - Preparation for v1.0**
- Week 1-2: Architecture planning (package split)
- Week 3-4: Breaking change documentation
- Release: v0.12.0 (pre-v1.0 candidate)

### Quarter 3-4 (Months 7-12)

**Months 7-9: v1.0.0-rc - Architecture Refactor**
- Package split implementation
- Migration tooling
- Extensive testing

**Months 10-12: v1.0.0 - Production Release**
- CRAN submission
- Publication (JOSS/journal)
- Community outreach
- Celebration! 🎉

---

## Resource Requirements

### Development Time (Estimated)

| Task | Effort | Timeline |
|------|--------|----------|
| Function modernization | 80 hours | 3 weeks |
| Testing infrastructure | 60 hours | 2 weeks |
| Benchmark suite | 20 hours | 1 week |
| Dependency audit | 40 hours | 2 weeks |
| pkgdown + vignettes | 60 hours | 3 weeks |
| Package split (v1.0) | 200 hours | 2-3 months |

**Total**: ~460 hours (~3 months full-time or 6 months part-time)

### Key Skills Needed

- R package development (essential)
- C++ programming (for optimizations)
- Python integration (reticulate)
- Technical writing (vignettes, papers)
- Web development (pkgdown customization)

---

## Risk Assessment

### High Risk

1. **Package split** - Breaking change, user migration complexity
   - Mitigation: Extensive communication, migration tools, grace period

2. **CRAN submission** - Strict requirements, dependency conflicts
   - Mitigation: Early checks, dependency reduction, documentation

### Medium Risk

3. **Test coverage** - Time-consuming, may uncover bugs
   - Mitigation: Incremental approach, prioritize critical functions

4. **Benchmark reproducibility** - Hardware-dependent results
   - Mitigation: Normalize by reference, document environment

### Low Risk

5. **Documentation website** - Tooling mature, low technical risk
   - Mitigation: Use pkgdown defaults, incremental customization

6. **Function modernization** - Pattern established, straightforward
   - Mitigation: Follow existing examples, systematic testing

---

## Success Criteria

### v0.10.x Series (Short-Term)
- ✅ 100% functions use av package
- ✅ 85%+ test coverage
- ✅ Benchmark suite published
- ✅ <30 dependencies in Imports
- ✅ pkgdown website deployed

### v1.0.0 (Long-Term)
- ✅ Modular architecture (core, ml, features)
- ✅ CRAN submission successful
- ✅ 50+ GitHub stars
- ✅ Published paper (JOSS or journal)
- ✅ 3+ external contributors
- ✅ 10+ citations per year

---

## Conclusion

The superassp package has a strong foundation and significant potential. By following this strategic roadmap, the package can:

1. **Achieve consistency** through function modernization
2. **Ensure reliability** through comprehensive testing
3. **Demonstrate performance** through formal benchmarking
4. **Reduce complexity** through dependency optimization
5. **Scale effectively** through modular architecture
6. **Build community** through documentation and outreach

**Recommended Immediate Focus**:
1. Function modernization (46% remaining) - Weeks 1-3
2. Testing infrastructure (85% target) - Weeks 4-6
3. Benchmark suite (formalize + publish) - Week 7

This balanced approach delivers immediate value while building toward long-term sustainability and growth.

---

**Next Review**: After v0.10.3 release (estimated 3 months)

---
