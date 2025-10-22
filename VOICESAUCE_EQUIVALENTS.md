# VoiceSauce Parameter Equivalents in superassp

## Overview

This document maps VoiceSauce parameter names to the proposed Titze 2015 / Nylén 2024 compliant names in superassp.

**VoiceSauce Documentation**: https://www.phonetics.ucla.edu/voicesauce/documentation/parameters.html

## Fundamental Frequency (fo)

| VoiceSauce | superassp Old | superassp New | Notes |
|------------|---------------|---------------|-------|
| `strF0` | `fo` | `fo[Hz]` | Straight F0 method |
| `sF0` | `fo` | `fo[Hz]` | Snack F0 method |
| `pF0` | `fo` | `fo[Hz]` | Praat F0 method |
| `shrF0` | `fo` | `fo[Hz]` | Subharmonic-to-harmonic F0 method |

**Standard**: All fo variants map to `fo[Hz]` (lowercase, with unit)

## Formants

| VoiceSauce | superassp Old | superassp New | Notes |
|------------|---------------|---------------|-------|
| `F1` | `F1` | `F1[Hz]` | First formant frequency |
| `F2` | `F2` | `F2[Hz]` | Second formant frequency |
| `F3` | `F3` | `F3[Hz]` | Third formant frequency |
| `F4` | `F4` | `F4[Hz]` | Fourth formant frequency |

**Standard**: Uppercase `F` + subscript number + `[Hz]` unit (Titze 2015)

## Formant Bandwidths

| VoiceSauce | superassp Old | superassp New | Notes |
|------------|---------------|---------------|-------|
| `B1` | `B1` | `B1[Hz]` | First formant bandwidth |
| `B2` | `B2` | `B2[Hz]` | Second formant bandwidth |
| `B3` | `B3` | `B3[Hz]` | Third formant bandwidth |

**Standard**: Uppercase `B` + subscript number + `[Hz]` unit (Titze 2015)

## Harmonic Amplitudes

| VoiceSauce | superassp Old | superassp New | Notes |
|------------|---------------|---------------|-------|
| `H1` | `H1` | `H1[dB]` | First harmonic (fundamental) |
| `H2` | `H2` | `H2[dB]` | Second harmonic |
| `H4` | `H4` | `H4[dB]` | Fourth harmonic |
| `H2K` | `H2K` | `H2k[dB]` | Harmonic nearest 2000 Hz (lowercase k) |
| `H5K` | `H5K` | `H5k[dB]` | Harmonic nearest 5000 Hz (lowercase k) |
| `A1` | `A1` | `A1[dB]` | Harmonic amplitude nearest F1 |
| `A2` | `A2` | `A2[dB]` | Harmonic amplitude nearest F2 |
| `A3` | `A3` | `A3[dB]` | Harmonic amplitude nearest F3 |

**Standard**:
- `Hi[dB]` where i = harmonic number (Titze 2015)
- `Ai[dB]` where i = formant proximity (Titze 2015)
- Use lowercase `k` for kHz (H2k not H2K)

## Harmonic Difference Measures (Uncorrected)

| VoiceSauce | superassp Old | superassp New | Notes |
|------------|---------------|---------------|-------|
| `H1-H2` or `H1H2u` | `H1H2` | `H1-H2u[dB]` | Open quotient correlate |
| `H2-H4` or `H2H4u` | `H2H4` | `H2-H4u[dB]` | Spectral slope |
| `H1-A1` or `H1A1u` | `H1A1` | `H1-A1u[dB]` | Spectral tilt at F1 |
| `H1-A2` or `H1A2u` | `H1A2` | `H1-A2u[dB]` | Spectral tilt at F2 |
| `H1-A3` or `H1A3u` | `H1A3` | `H1-A3u[dB]` | Spectral tilt at F3 |
| `H4-2K` or `H42Ku` | `H4H2K` | `H4-H2ku[dB]` | Mid-frequency tilt |
| `2K-5K` or `2K5Ku` | `H2KH5K` | `H2k-H5ku[dB]` | High-frequency tilt |

**Standard**:
- Use hyphen separator: `Hi-Hj` (Titze 2015)
- Append `u` suffix for uncorrected (Nylén 2024)
- Add `[dB]` unit after suffix: `H1-H2u[dB]`

## Harmonic Difference Measures (Corrected)

| VoiceSauce | superassp Old | superassp New | Notes |
|------------|---------------|---------------|-------|
| `H1H2c` | `H1H2c` | `H1-H2c[dB]` | Formant-corrected H1-H2 |
| `H2H4c` | `H2H4c` | `H2-H4c[dB]` | Formant-corrected H2-H4 |
| `H1A1c` | `H1A1c` | `H1-A1c[dB]` | Corrected H1-A1 |
| `H1A2c` | `H1A2c` | `H1-A2c[dB]` | Corrected H1-A2 |
| `H1A3c` | `H1A3c` | `H1-A3c[dB]` | Corrected H1-A3 |
| `H42Kc` | `H4H2Kc` | `H4-H2kc[dB]` | Corrected H4-2K |
| `2K5Kc` | `H2KH5Kc` | `H2k-H5kc[dB]` | Corrected 2K-5K |

**Standard** (Nylén 2024 Extension):
- Use hyphen separator: `Hi-Hj`
- Append `c` suffix for corrected: `H1-H2c`
- Add `[dB]` unit after suffix: `H1-H2c[dB]`
- **Rationale**: Corrected measures remove formant influence, used in gender/voice analysis

## Voice Quality Measures

| VoiceSauce | superassp Old | superassp New | Notes |
|------------|---------------|---------------|-------|
| `CPP` | `CPP` | `CPP[dB]` | Cepstral Peak Prominence |
| `HNR05` | `HNR05` | `HNR05[dB]` | Harmonics-to-Noise Ratio (0-500 Hz) |
| `HNR15` | `HNR15` | `HNR15[dB]` | Harmonics-to-Noise Ratio (0-1500 Hz) |
| `HNR25` | `HNR25` | `HNR25[dB]` | Harmonics-to-Noise Ratio (0-2500 Hz) |
| `HNR35` | `HNR35` | `HNR35[dB]` | Harmonics-to-Noise Ratio (0-3500 Hz) |
| `SHR` | `SHR` | `SHR[dB]` | Subharmonic-to-Harmonic Ratio |
| `Energy` | `Energy` | `Energy[dB]` | RMS energy |

**Standard**:
- Add `[dB]` unit to all amplitude-based measures
- Keep VoiceSauce abbreviations (widely recognized)

## Implementation in superassp

### Function: `trk_praat_sauce()`

**File**: `R/ssff_python_pm_psauce.R`

**Current output tracks**:
```r
attr(trk_praat_sauce, "tracks") <- c(
  "t", "fo", "F1", "F2", "F3", "B1", "B2", "B3",
  "H1", "H2", "H4", "A1", "A2", "A3", "H2K", "H5K",
  "H1H2", "H2H4", "H1A1", "H1A2", "H1A3", "H4H2K", "H2KH5K",
  "H1H2c", "H1A1c", "H1A2c", "H1A3c",
  "H1H2u", "H1A1u", "H1A2u", "H1A3u",
  "CPP", "HNR05", "HNR15", "HNR25", "HNR35", "SHR", "Energy"
)
```

**Proposed new tracks**:
```r
attr(trk_praat_sauce, "tracks") <- c(
  "t", "fo[Hz]", "F1[Hz]", "F2[Hz]", "F3[Hz]", "B1[Hz]", "B2[Hz]", "B3[Hz]",
  "H1[dB]", "H2[dB]", "H4[dB]", "A1[dB]", "A2[dB]", "A3[dB]", "H2k[dB]", "H5k[dB]",
  "H1-H2[dB]", "H2-H4[dB]", "H1-A1[dB]", "H1-A2[dB]", "H1-A3[dB]", "H4-H2k[dB]", "H2k-H5k[dB]",
  "H1-H2c[dB]", "H1-A1c[dB]", "H1-A2c[dB]", "H1-A3c[dB]",
  "H1-H2u[dB]", "H1-A1u[dB]", "H1-A2u[dB]", "H1-A3u[dB]",
  "CPP[dB]", "HNR05[dB]", "HNR15[dB]", "HNR25[dB]", "HNR35[dB]", "SHR[dB]", "Energy[dB]"
)
```

**Changes**:
1. `fo` → `fo[Hz]` (lowercase + unit)
2. `F1`, `F2`, `F3` → `F1[Hz]`, `F2[Hz]`, `F3[Hz]` (add units)
3. `B1`, `B2`, `B3` → `B1[Hz]`, `B2[Hz]`, `B3[Hz]` (add units)
4. `H1`, `H2`, etc. → `H1[dB]`, `H2[dB]`, etc. (add units)
5. `H2K`, `H5K` → `H2k[dB]`, `H5k[dB]` (lowercase k + units)
6. `H1H2` → `H1-H2[dB]` (add hyphen + unit)
7. `H1H2c` → `H1-H2c[dB]` (add hyphen + unit)
8. `H1H2u` → `H1-H2u[dB]` (add hyphen + unit)
9. All quality measures get `[dB]` unit

## Backwards Compatibility

During the deprecation period, both old and new names will be supported via the alias system:

```r
# Example: User code with old names still works
result <- trk_praat_sauce("audio.wav", toFile = FALSE)

# Both work:
f0_old <- result$f0        # Deprecated, shows warning
f0_new <- result$`fo[Hz]`  # New standard name

# as.data.frame() handles both automatically
df <- as.data.frame(result)
# Column names will be: fo[Hz], F1[Hz], F2[Hz], ... (new names)
# But df$f0 still works with deprecation warning
```

## Benefits for VoiceSauce Users

1. **Familiar notation**: Core abbreviations (H1, A1, CPP, etc.) preserved
2. **Explicit units**: No ambiguity about dB vs linear, Hz vs kHz
3. **Scientific standard**: Aligns with Titze 2015 JASA publication
4. **Corrected formants**: Nylén 2024 extension for gender/voice research
5. **Automatic conversion**: Integration with R `units` package
6. **Cross-platform**: Matches notation in Praat, phonetics literature

## Migration Guide for VoiceSauce Users

### Old workflow:
```r
library(superassp)

# Extract VoiceSauce parameters
result <- trk_praat_sauce("speech.wav", toFile = FALSE)
df <- as.data.frame(result)

# Access parameters
f0 <- df$f0
h1h2c <- df$H1H2c
```

### New workflow (after migration):
```r
library(superassp)

# Extract VoiceSauce parameters
result <- trk_praat_sauce("speech.wav", toFile = FALSE)
df <- as.data.frame(result)

# Access parameters with new names
f0 <- df$`fo[Hz]`              # Note: brackets require backticks
h1h2c <- df$`H1-H2c[dB]`       # Hyphen separator + unit

# Or use without units package (returns numeric):
f0_numeric <- as.numeric(df$`fo[Hz]`)
```

### During deprecation period (both work):
```r
# Old names show deprecation warning but still work
f0 <- df$f0           # Works, but warns: "Use 'fo[Hz]' instead"
h1h2c <- df$H1H2c     # Works, but warns: "Use 'H1-H2c[dB]' instead"

# New names are the standard
f0 <- df$`fo[Hz]`     # No warning
h1h2c <- df$`H1-H2c[dB]`  # No warning
```

## References

- VoiceSauce Documentation: https://www.phonetics.ucla.edu/voicesauce/documentation/parameters.html
- Titze et al. (2015): "Toward a consensus on symbolic notation..." JASA 137(5), 3005-3007
- Nylén et al. (2024): "Acoustic cues to femininity and masculinity..." JASA 155(2), 1373-1387
- Shue et al. (2011): "VoiceSauce: A program for voice analysis" Proc. ICPhS XVII

## Summary Table: Complete VoiceSauce Mapping

| Category | VoiceSauce | Old | New | Change |
|----------|------------|-----|-----|--------|
| F0 | `strF0`, `sF0`, `pF0`, `shrF0` | `fo` | `fo[Hz]` | Add unit |
| Formants | `F1` - `F4` | `F1` - `F4` | `F1[Hz]` - `F4[Hz]` | Add unit |
| Bandwidths | `B1` - `B3` | `B1` - `B3` | `B1[Hz]` - `B3[Hz]` | Add unit |
| Harmonics | `H1`, `H2`, `H4` | Same | `H1[dB]`, `H2[dB]`, `H4[dB]` | Add unit |
| Harmonics | `H2K`, `H5K` | Same | `H2k[dB]`, `H5k[dB]` | Lowercase k + unit |
| Harmonics | `A1`, `A2`, `A3` | Same | `A1[dB]`, `A2[dB]`, `A3[dB]` | Add unit |
| Differences | `H1-H2`, etc. | `H1H2` | `H1-H2[dB]` | Add hyphen + unit |
| Corrected | `H1H2c`, etc. | Same | `H1-H2c[dB]` | Add hyphen + unit |
| Uncorrected | `H1H2u`, etc. | Same | `H1-H2u[dB]` | Add hyphen + unit |
| Voice Quality | `CPP`, `HNR*`, `SHR` | Same | `CPP[dB]`, `HNR*[dB]`, `SHR[dB]` | Add unit |
| Energy | `Energy` | Same | `Energy[dB]` | Add unit |

**Total**: 40+ VoiceSauce parameters mapped to standardized notation
