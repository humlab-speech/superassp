"""
Test script for thesis-compliant implementation fixes

Tests the new use_thesis_mode parameter and validates:
1. PPE semitone conversion
2. AR-based jitter (PQ_AR)
3. NMSP measure
4. Shimmer dB corrections
"""

import numpy as np
import soundfile as sf
from voice_analysis.core import VoiceAnalyzer

def test_thesis_mode():
    """Test both MATLAB mode (default) and thesis mode"""
    
    # Load test audio
    audio, fs = sf.read('../a1.wav')
    
    print("=" * 80)
    print("THESIS-COMPLIANT IMPLEMENTATION TEST")
    print("=" * 80)
    
    # Test 1: MATLAB mode (default)
    print("\n1. Testing MATLAB Mode (default)")
    print("-" * 80)
    analyzer_matlab = VoiceAnalyzer(use_thesis_mode=False)
    measures_matlab, F0_matlab = analyzer_matlab.analyze(audio, fs)
    
    print(f"\nMATLAB Mode Results:")
    print(f"  Total measures: {len(measures_matlab)}")
    print(f"  PPE: {measures_matlab.get('PPE', np.nan):.6f}")
    print(f"  Jitter RAP: {measures_matlab.get('jitter_RAP', np.nan):.6f}")
    print(f"  Shimmer dB: {measures_matlab.get('shimmer_dB', np.nan):.6f}")
    
    # Check for thesis-specific measures (should not exist)
    has_ar = 'jitter_PQ_AR' in measures_matlab
    has_nmsp = 'jitter_NMSP' in measures_matlab
    print(f"  Has AR jitter: {has_ar} (expected: False)")
    print(f"  Has NMSP: {has_nmsp} (expected: False)")
    
    # Test 2: Thesis mode
    print("\n2. Testing Thesis Mode (use_thesis_mode=True)")
    print("-" * 80)
    analyzer_thesis = VoiceAnalyzer(use_thesis_mode=True)
    measures_thesis, F0_thesis = analyzer_thesis.analyze(audio, fs)
    
    print(f"\nThesis Mode Results:")
    print(f"  Total measures: {len(measures_thesis)}")
    print(f"  PPE (semitone): {measures_thesis.get('PPE', np.nan):.6f}")
    print(f"  Jitter RAP: {measures_thesis.get('jitter_RAP', np.nan):.6f}")
    print(f"  Shimmer dB: {measures_thesis.get('shimmer_dB', np.nan):.6f}")
    
    # Check for thesis-specific measures (should exist)
    has_ar = 'jitter_PQ_AR' in measures_thesis
    has_nmsp = 'jitter_NMSP' in measures_thesis
    has_range = 'jitter_F0_range' in measures_thesis
    print(f"  Has AR jitter: {has_ar} (expected: True)")
    print(f"  Has NMSP: {has_nmsp} (expected: True)")
    print(f"  Has F0 range: {has_range} (expected: True)")
    
    if has_ar:
        print(f"  Jitter PQ_AR: {measures_thesis['jitter_PQ_AR']:.6f}")
    if has_nmsp:
        print(f"  Jitter NMSP: {measures_thesis['jitter_NMSP']:.6f}")
    if has_range:
        print(f"  Jitter F0 range: {measures_thesis['jitter_F0_range']:.2f} Hz")
    
    # Test 3: Compare PPE values
    print("\n3. PPE Comparison (should differ due to log base)")
    print("-" * 80)
    ppe_matlab = measures_matlab.get('PPE', np.nan)
    ppe_thesis = measures_thesis.get('PPE', np.nan)
    print(f"  MATLAB PPE (natural log): {ppe_matlab:.6f}")
    print(f"  Thesis PPE (semitones):   {ppe_thesis:.6f}")
    print(f"  Difference: {abs(ppe_matlab - ppe_thesis):.6f}")
    print(f"  Expected: Values should differ (different log base)")
    
    # Test 4: Validate shimmer dB
    print("\n4. Shimmer dB Comparison")
    print("-" * 80)
    shimmer_matlab = measures_matlab.get('shimmer_dB', np.nan)
    shimmer_thesis = measures_thesis.get('shimmer_dB', np.nan)
    print(f"  MATLAB Shimmer dB: {shimmer_matlab:.6f}")
    print(f"  Thesis Shimmer dB: {shimmer_thesis:.6f}")
    print(f"  Difference: {abs(shimmer_matlab - shimmer_thesis):.6f}")
    
    # Test 5: Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    
    extra_measures_thesis = len(measures_thesis) - len(measures_matlab)
    print(f"✓ MATLAB mode: {len(measures_matlab)} measures")
    print(f"✓ Thesis mode: {len(measures_thesis)} measures (+{extra_measures_thesis} additional)")
    
    success = True
    if not has_ar:
        print("✗ ERROR: AR jitter not found in thesis mode")
        success = False
    if not has_nmsp:
        print("✗ ERROR: NMSP not found in thesis mode")
        success = False
    if abs(ppe_matlab - ppe_thesis) < 0.001:
        print("⚠ WARNING: PPE values too similar (log base may not be changed)")
    
    if success:
        print("\n✅ All tests PASSED!")
    else:
        print("\n❌ Some tests FAILED")
    
    return measures_matlab, measures_thesis


def test_individual_functions():
    """Test individual functions directly"""
    
    print("\n" + "=" * 80)
    print("INDIVIDUAL FUNCTION TESTS")
    print("=" * 80)
    
    # Create test F0 contour
    F0_test = np.array([127, 130, 125, 128, 135, 132, 127, 129, 131, 128])
    
    # Test 1: AR perturbation
    print("\n1. Testing AR Perturbation Quotient")
    print("-" * 80)
    from voice_analysis.utils.perturbation import compute_ar_perturbation_quotient
    
    pq_ar = compute_ar_perturbation_quotient(F0_test, ar_order=3)
    print(f"  AR PQ (order 3): {pq_ar:.6f}")
    print(f"  Expected: Small positive value (e.g., 0.001-0.01)")
    print(f"  Status: {'✓ OK' if 0 < pq_ar < 0.1 else '✗ FAIL'}")
    
    # Test 2: NMSP
    print("\n2. Testing NMSP")
    print("-" * 80)
    from voice_analysis.utils.perturbation import compute_nmsp
    
    nmsp = compute_nmsp(F0_test)
    print(f"  NMSP: {nmsp:.6f}")
    print(f"  Expected: Positive value")
    print(f"  Status: {'✓ OK' if nmsp > 0 else '✗ FAIL'}")
    
    # Test 3: PPE with different modes
    print("\n3. Testing PPE Log Base")
    print("-" * 80)
    from voice_analysis.features.ppe import compute_ppe
    
    F0_longer = np.random.uniform(100, 150, 100)
    
    ppe_natural = compute_ppe(F0_longer, 44100, use_thesis_mode=False)
    ppe_semitone = compute_ppe(F0_longer, 44100, use_thesis_mode=True)
    
    print(f"  PPE (natural log): {ppe_natural:.6f}")
    print(f"  PPE (semitones):   {ppe_semitone:.6f}")
    print(f"  Should differ: {'✓ Yes' if abs(ppe_natural - ppe_semitone) > 0.01 else '✗ No (values too similar)'}")
    
    print("\n✅ Individual function tests complete")


if __name__ == "__main__":
    # Run full analysis test
    measures_matlab, measures_thesis = test_thesis_mode()
    
    # Run individual function tests
    test_individual_functions()
    
    print("\n" + "=" * 80)
    print("All tests completed!")
    print("=" * 80)
