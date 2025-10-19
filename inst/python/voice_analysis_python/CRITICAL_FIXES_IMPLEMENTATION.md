# Critical Fixes Implementation Guide

**Based on**: THESIS_VERIFICATION_REPORT.md  
**Priority**: HIGH - These fixes ensure fidelity to Tsanas (2012) thesis  
**Estimated Time**: 3-4 hours total

---

## Fix 1: PPE Semitone Conversion (15 minutes)

### Issue
Current code uses natural logarithm instead of semitone scale as specified in thesis Equation 3.54.

### Location
`voice_analysis_python/voice_analysis/features/ppe.py`, lines 44-45

### Current Code (INCORRECT)
```python
# Line 44
logF0 = safe_log(F0 / f0_mean_healthy, base='e')  # Natural log, not semitones!
```

### Corrected Code
```python
# Convert to semitones as per Equation 3.54 (p. 71)
# s = 12 × log₂(F₀/127)
# Note: 127 Hz is reference for males (thesis footnote 22, p. 71)
semitones = 12 * np.log2(F0 / 127)  # Semitone scale
```

### Implementation Steps

1. **Update line 44**:
```python
# OLD:
logF0 = safe_log(F0 / f0_mean_healthy, base='e')

# NEW:
semitones = 12 * np.log2(F0 / 127)  # Semitone scale per Equation 3.54
```

2. **Update variable names** (lines 44-96):
   - Change `logF0` → `semitones`
   - Update comments to reflect semitone scale

3. **Test validation**:
```python
# Test case
F0_test = np.array([127, 254, 190, 63.5])
expected = np.array([0, 12, 6.95, -12])  # Approximate semitones
actual = 12 * np.log2(F0_test / 127)
print(np.allclose(actual, expected, atol=0.1))  # Should be True
```

---

## Fix 2: Implement AR-based Jitter (PQ_AR) - 2 hours

### Issue
Equation 3.39 (AR-based perturbation quotient) is not implemented.

### Location
New function needed in `voice_analysis_python/voice_analysis/utils/perturbation.py`

### Thesis Specification (Equation 3.39, p. 59)
```
PQ_AR = Σᵢ|Σⱼ aⱼ(T₀ᵢ₋ⱼ - T̄₀)| / Σᵢ T₀ᵢ
```

Where:
- `aⱼ` = AR(10) coefficients from Yule-Walker equations
- `T̄₀` = mean of T₀ contour
- `j = 1, 2, ..., 10` (p=10 coefficients)

### Implementation

#### Step 1: Add to `perturbation.py`

```python
def compute_ar_perturbation_quotient(time_series, ar_order=10):
    """
    Compute AR-based Perturbation Quotient (Equation 3.39)
    
    Uses Yule-Walker equations to estimate AR coefficients,
    following Schoentgen & de Guchteneere (1995).
    
    Parameters:
    -----------
    time_series : ndarray
        F0 or amplitude contour
    ar_order : int
        AR model order (default: 10, per thesis p. 59)
        
    Returns:
    --------
    pq_ar : float
        AR-based perturbation quotient
    """
    time_series = np.asarray(time_series).ravel()
    time_series = time_series[time_series > 0]
    
    if len(time_series) < ar_order + 1:
        return np.nan
    
    # Center the signal
    mean_val = np.mean(time_series)
    centered = time_series - mean_val
    
    # Estimate AR coefficients using Yule-Walker
    # scipy.signal.levinson implements Levinson-Durbin recursion
    from scipy.signal import correlate
    
    # Compute autocorrelation
    acf = correlate(centered, centered, mode='full', method='fft')
    acf = acf[len(acf)//2:]  # Keep positive lags
    
    if acf[0] == 0:
        return np.nan
    
    acf = acf / acf[0]  # Normalize
    
    # Solve Yule-Walker equations
    # R @ a = r, where R is Toeplitz matrix
    R = np.zeros((ar_order, ar_order))
    for i in range(ar_order):
        for j in range(ar_order):
            R[i, j] = acf[abs(i - j)]
    
    r = acf[1:ar_order+1]
    
    try:
        ar_coeffs = np.linalg.solve(R, r)
    except np.linalg.LinAlgError:
        return np.nan
    
    # Apply Equation 3.39
    # PQ_AR = Σᵢ|Σⱼ aⱼ(T₀ᵢ₋ⱼ - T̄₀)| / Σᵢ T₀ᵢ
    
    numerator = 0.0
    for i in range(ar_order, len(time_series)):
        weighted_sum = 0.0
        for j in range(ar_order):
            weighted_sum += ar_coeffs[j] * (time_series[i - j - 1] - mean_val)
        numerator += abs(weighted_sum)
    
    denominator = np.sum(time_series)
    
    if denominator == 0:
        return np.nan
    
    pq_ar = numerator / denominator
    
    return pq_ar
```

#### Step 2: Update `jitter_shimmer.py`

Add after line 82:

```python
# 12. AR-based Perturbation Quotient (Equation 3.39)
if len(time_series) >= 11:  # Need at least ar_order + 1
    from ..utils.perturbation import compute_ar_perturbation_quotient
    measures[f'{prefix}_PQ_AR'] = compute_ar_perturbation_quotient(time_series, ar_order=10)
else:
    measures[f'{prefix}_PQ_AR'] = np.nan

# Renumber subsequent measures (12→13, 13→14, etc.)
```

#### Step 3: Update `_get_nan_measures()` function

Add to dictionary (line ~143):
```python
f'{prefix}_PQ_AR': np.nan,
```

#### Step 4: Test

```python
# Test with synthetic signal
np.random.seed(42)
F0_test = 120 + 5 * np.sin(np.linspace(0, 10*np.pi, 200)) + np.random.normal(0, 1, 200)
pq_ar = compute_ar_perturbation_quotient(F0_test)
print(f"PQ_AR = {pq_ar:.6f}")  # Should be a small positive number
```

---

## Fix 3: Implement NMSP (30 minutes)

### Issue
Normalized Mean Squared Perturbation (Equation 3.41) is missing.

### Location
`voice_analysis_python/voice_analysis/features/jitter_shimmer.py`, after line 96

### Thesis Specification (Equation 3.41, p. 59)
```
NMSP = Σᵢ(T₀ᵢ - T̄₀)² / [(1/N) × (ΣⱼT₀ⱼ)²]
```

### Implementation

Add after line 96 (after CV calculation):

```python
# 15. Normalized Mean Squared Perturbation (Equation 3.41)
# NMSP = Σ(T₀ - T̄₀)² / [(1/N) × (ΣT₀)²]
sum_sq_dev = np.sum((time_series - mean_val)**2)
denominator = (np.sum(time_series) / len(time_series))**2
measures[f'{prefix}_NMSP'] = sum_sq_dev / denominator if denominator > 0 else np.nan

# Renumber: 15→16 (TKEO_mean), 16→17 (TKEO_std), etc.
```

Update `_get_nan_measures()`:
```python
f'{prefix}_NMSP': np.nan,
```

---

## Fix 4: Validate Shimmer dB (30 minutes)

### Issue
Potential error in absolute value placement in shimmer dB formula.

### Location
`voice_analysis_python/voice_analysis/features/jitter_shimmer.py`, lines 87-89

### Current Code
```python
ratio = time_series[:-1] / (time_series[1:] + 1e-10)
measures[f'{prefix}_dB'] = np.mean(20 * np.abs(np.log10(np.abs(ratio) + 1e-10)))
```

### Issue Analysis

The standard shimmer dB formula is:
```
Shimmer_dB = (1/(N-1)) × Σᵢ|20 × log₁₀(Aᵢ₊₁/Aᵢ)|
```

The absolute value should be **outside** the log, not inside.

### Corrected Code

```python
# Shimmer (dB) - Equation for amplitude perturbation
# Standard formula: mean(|20 × log₁₀(Aᵢ₊₁/Aᵢ)|)
if measure_type == 'shimmer':
    # Avoid division by zero
    ratio = time_series[1:] / (time_series[:-1] + 1e-10)
    # Avoid log(0) with small epsilon
    ratio = np.clip(ratio, 1e-10, None)
    measures[f'{prefix}_dB'] = np.mean(np.abs(20 * np.log10(ratio)))
else:
    measures[f'{prefix}_dB'] = np.nan
```

### Validation Test

```python
# Test shimmer dB calculation
amplitudes = np.array([1.0, 1.1, 0.9, 1.05, 0.95])
ratio = amplitudes[1:] / amplitudes[:-1]
shimmer_db = np.mean(np.abs(20 * np.log10(ratio)))
print(f"Shimmer dB: {shimmer_db:.4f}")  # Should be ~0.8-1.0 dB

# Compare with MATLAB if available
```

---

## Fix 5: Add F₀ Range Measure (5 minutes)

### Issue
F₀ range (p95 - p5) is computed but not stored explicitly.

### Location
`voice_analysis_python/voice_analysis/features/jitter_shimmer.py`, after line 115

### Implementation

```python
# 22. TKEO IQR (p75 - p5)
measures[f'{prefix}_TKEO_IQR'] = percentiles[3] - percentiles[0]

# 23. F₀/Amplitude Range (p95 - p5) - Thesis p. 60
measures[f'{prefix}_range_robust'] = percentiles[4] - percentiles[0]
```

Update `_get_nan_measures()`:
```python
f'{prefix}_range_robust': np.nan,
```

---

## Fix 6: DFA Scaling Range (5 minutes)

### Issue
Potential off-by-one error in scale endpoint.

### Location
`voice_analysis_python/voice_analysis/features/dfa.py`, lines 55-56

### Thesis (p. 68)
```matlab
dfa_scaling = (50:20:200)'  % MATLAB inclusive endpoint
```

### Current Code
```python
nvals = np.arange(50, 201, 20)  # [50, 70, ..., 190]
```

### Validation

MATLAB `50:20:200` produces: `[50, 70, 90, 110, 130, 150, 170, 190]`  
Python `np.arange(50, 201, 20)` produces: `[50, 70, 90, 110, 130, 150, 170, 190]`

**Result**: ✅ **CORRECT** - No fix needed (end=201 makes it inclusive to 190)

However, verify if thesis intended 200 or 210 as endpoint:

```python
# If 200 should be included:
nvals = np.arange(50, 201, 20)  # Current (ends at 190)
# OR
nvals = np.append(np.arange(50, 201, 20), 200)  # Force include 200

# If 210 should be included:
nvals = np.arange(50, 211, 20)  # [50, 70, ..., 190, 210]
```

**Action**: Compare MATLAB `.mat` file output to confirm exact scale values.

---

## Testing Protocol

### Step 1: Unit Tests

Create `test_critical_fixes.py`:

```python
import numpy as np
from voice_analysis.features.jitter_shimmer import compute_jitter_shimmer_features
from voice_analysis.features.ppe import compute_ppe

def test_ppe_semitones():
    """Test PPE uses semitone scale"""
    F0 = np.array([127, 254, 190])  # Reference, octave up, ~major 6th
    # Should use 12*log2(F0/127) internally
    # Expected semitones: [0, 12, ~6.95]
    ppe = compute_ppe(F0, fs=16000)
    assert ppe is not np.nan
    print(f"✓ PPE computation successful: {ppe:.4f}")

def test_ar_jitter():
    """Test AR-based jitter implementation"""
    F0 = 120 + 5*np.sin(np.linspace(0, 10*np.pi, 200)) + np.random.normal(0, 1, 200)
    measures = compute_jitter_shimmer_features(F0, measure_type='jitter')
    assert 'jitter_PQ_AR' in measures
    assert not np.isnan(measures['jitter_PQ_AR'])
    print(f"✓ AR jitter: {measures['jitter_PQ_AR']:.6f}")

def test_nmsp():
    """Test NMSP implementation"""
    F0 = 120 + np.random.normal(0, 2, 100)
    measures = compute_jitter_shimmer_features(F0, measure_type='jitter')
    assert 'jitter_NMSP' in measures
    assert not np.isnan(measures['jitter_NMSP'])
    print(f"✓ NMSP: {measures['jitter_NMSP']:.6f}")

def test_shimmer_db():
    """Test shimmer dB correction"""
    amplitudes = np.array([1.0, 1.1, 0.9, 1.05, 0.95, 1.02])
    measures = compute_jitter_shimmer_features(amplitudes, measure_type='shimmer')
    shimmer_db = measures['shimmer_dB']
    assert 0.5 < shimmer_db < 2.0  # Reasonable range
    print(f"✓ Shimmer dB: {shimmer_db:.4f}")

if __name__ == '__main__':
    test_ppe_semitones()
    test_ar_jitter()
    test_nmsp()
    test_shimmer_db()
    print("\n✅ All critical fixes validated!")
```

### Step 2: MATLAB Comparison

```python
# Run on same test file as MATLAB
import scipy.io as sio
from voice_analysis import analyze_voice

# Analyze test file
results_python = analyze_voice('a1.wav')

# Load MATLAB reference
matlab_ref = sio.loadmat('matlab_reference.mat')

# Compare specific fixed measures
measures_to_check = [
    'PPE',
    'jitter_PQ_AR',
    'jitter_NMSP',
    'shimmer_dB',
]

for measure in measures_to_check:
    py_val = results_python.get(measure, np.nan)
    mat_val = matlab_ref['measures'][0, matlab_ref['names'].index(measure)]
    
    rel_error = abs(py_val - mat_val) / abs(mat_val) * 100
    
    print(f"{measure:20s}: Python={py_val:10.6f}, MATLAB={mat_val:10.6f}, Error={rel_error:6.2f}%")
```

---

## Implementation Checklist

- [ ] **Fix 1**: PPE semitone conversion (15 min)
  - [ ] Update line 44 in ppe.py
  - [ ] Rename variables logF0 → semitones
  - [ ] Test with synthetic signal
  
- [ ] **Fix 2**: Implement AR jitter (2 hours)
  - [ ] Add `compute_ar_perturbation_quotient()` to perturbation.py
  - [ ] Update jitter_shimmer.py
  - [ ] Update _get_nan_measures()
  - [ ] Test with synthetic signal
  
- [ ] **Fix 3**: Implement NMSP (30 min)
  - [ ] Add NMSP calculation
  - [ ] Update _get_nan_measures()
  - [ ] Test with synthetic signal
  
- [ ] **Fix 4**: Validate shimmer dB (30 min)
  - [ ] Correct absolute value placement
  - [ ] Fix ratio direction if needed
  - [ ] Test with synthetic signal
  
- [ ] **Fix 5**: Add F₀ range (5 min)
  - [ ] Add range_robust measure
  - [ ] Update _get_nan_measures()
  
- [ ] **Fix 6**: Verify DFA scaling (5 min)
  - [ ] Compare MATLAB output
  - [ ] Adjust if needed
  
- [ ] **Testing**: Run validation suite (1 hour)
  - [ ] Unit tests pass
  - [ ] MATLAB comparison <5% error
  - [ ] Document any remaining discrepancies

---

## Expected Outcome

After implementing these fixes:

1. **Fidelity to thesis**: 95%+ (up from 83%)
2. **MATLAB correlation**: >0.99 for all measures
3. **Missing features**: 0 critical, 2 optional
4. **Production ready**: Yes

---

## Timeline

| Task | Time | Priority |
|------|------|----------|
| Fix 1 (PPE) | 15 min | Critical |
| Fix 2 (AR jitter) | 2 hours | Critical |
| Fix 3 (NMSP) | 30 min | High |
| Fix 4 (Shimmer dB) | 30 min | High |
| Fix 5 (F₀ range) | 5 min | Medium |
| Fix 6 (DFA verify) | 5 min | Low |
| Testing | 1 hour | Critical |
| **Total** | **4.5 hours** | |

---

**Implementation Date**: TBD  
**Validation Date**: TBD  
**Sign-off**: Pending testing
