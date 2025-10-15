# Units of Measurement in Voice Analysis Functions

This document specifies the units of measurement for all voice analysis parameters returned by Praat/Parselmouth-based functions.

## praat_voice_report / praat_voice_report_opt

Based on Praat Voice Report manual and R documentation.

### Pitch Measurements

| Parameter | Unit | Notes |
|-----------|------|-------|
| Median pitch | Hz | Hertz |
| Mean pitch | Hz | Hertz |
| Standard deviation | Hz | Standard deviation of pitch |
| Minimum pitch | Hz | Hertz |
| Maximum pitch | Hz | Hertz |

### Pulse and Period Measurements

| Parameter | Unit | Notes |
|-----------|------|-------|
| Number of pulses | count | Integer count |
| Number of periods | count | Integer count |
| Mean period | seconds | Average period length |
| Standard deviation of period | seconds | SD of period length |

**Important**: Praat reports periods in **seconds**, not milliseconds.

### Voicing Measurements

| Parameter | Unit | Notes |
|-----------|------|-------|
| Fraction of locally unvoiced frames | fraction | 0.0 to 1.0 (unitless) |
| Number of voice breaks | count | Integer count |
| Degree of voice breaks | fraction | Voice breaks / frames (unitless) |

### Jitter Measurements

| Parameter | Unit | Notes |
|-----------|------|-------|
| Jitter (local) | % | Percentage (e.g., 1.5 means 1.5%) |
| Jitter (local, absolute) | seconds | **NOT microseconds** - Praat returns seconds |
| Jitter (rap) | % | Relative Average Perturbation |
| Jitter (ppq5) | % | Five-point Period Perturbation Quotient |
| Jitter (ddp) | % | Difference of Differences of Periods |

**Critical Note**: While some literature reports jitter (local, absolute) in microseconds
(μs), Praat's Voice Report actually returns this value in **seconds**. The threshold
value of 83.200 μs mentioned in literature is 0.0000832 seconds in Praat output.

Reference: https://www.fon.hum.uva.nl/praat/manual/Voice_2__Jitter.html

### Shimmer Measurements

| Parameter | Unit | Notes |
|-----------|------|-------|
| Shimmer (local) | % | Percentage |
| Shimmer (local, dB) | dB | Decibels |
| Shimmer (apq3) | % | Three-point APQ |
| Shimmer (apq5) | % | Five-point APQ |
| Shimmer (apq11) | % | Eleven-point APQ |
| Shimmer (dda) | % | Difference of Differences of Amplitudes |

Reference: https://www.fon.hum.uva.nl/praat/manual/Voice_3__Shimmer.html

### Harmonicity Measurements

| Parameter | Unit | Notes |
|-----------|------|-------|
| Mean autocorrelation | fraction | 0.0 to 1.0 (unitless) |
| Mean noise-to-harmonics ratio | ratio | Unitless ratio (NHR) |
| Mean harmonics-to-noise ratio | dB | Decibels (HNR) |

**Note**:
- NHR is the inverse of HNR and is reported as a ratio (unitless)
- HNR is reported in decibels (dB)
- Both are derived from the same underlying measurement

## Intensity Measurements

(Not included in standard voice report output, but computed in Praat script)

| Parameter | Unit | Notes |
|-----------|------|-------|
| Mean intensity | dB | Decibels SPL |
| Median intensity | dB | Decibels SPL |
| SD intensity | dB | Standard deviation |

## Common Confusions and Clarifications

### 1. Jitter (local, absolute) - Seconds vs Microseconds

**In Literature**: Often reported in microseconds (μs)
- Example threshold: 83.200 μs (MDVP)

**In Praat Output**: Reported in seconds (s)
- Same threshold: 0.0000832 s

**In this package**: Values are returned in **seconds** (as Praat provides them)

**Conversion**:
```r
jitter_us <- jitter_seconds * 1e6  # Convert seconds to microseconds
```

### 2. Period Measurements - Always in Seconds

All period-related measurements are in seconds:
- Mean period: seconds
- Standard deviation of period: seconds
- Jitter (local, absolute): seconds

### 3. Percentage vs Ratio

**Percentages** (returned as values like 1.5 for 1.5%):
- All jitter measurements except (local, absolute)
- All shimmer measurements except (local, dB)

**Ratios** (returned as values between 0 and 1):
- Fraction of locally unvoiced frames
- Degree of voice breaks
- Mean autocorrelation

**Special case**:
- NHR (noise-to-harmonics ratio) is a ratio but can exceed 1.0

### 4. Decibel Measurements

**Always in dB**:
- Shimmer (local, dB)
- Mean harmonics-to-noise ratio (HNR)
- All intensity measurements

## Verification Against Praat Manual

These units have been verified against:
1. Praat manual: https://www.fon.hum.uva.nl/praat/manual/Manual.html
2. Voice report specific pages
3. Praat script output format
4. R function documentation in superassp

## Implementation Notes

When implementing Parselmouth-based functions:

1. **Do not convert units** - return values exactly as Praat provides them
2. **Document carefully** - especially for jitter (local, absolute)
3. **Match original behavior** - maintain exact compatibility with file-based functions
4. **Test equivalence** - compare Parselmouth output with external Praat output

## References

- Praat Manual: https://www.fon.hum.uva.nl/praat/manual/Manual.html
- Voice Jitter: https://www.fon.hum.uva.nl/praat/manual/Voice_2__Jitter.html
- Voice Shimmer: https://www.fon.hum.uva.nl/praat/manual/Voice_3__Shimmer.html
- MDVP Normative Data (showing μs convention): Various clinical voice literature
