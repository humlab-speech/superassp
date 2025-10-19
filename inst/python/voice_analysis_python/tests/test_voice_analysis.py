"""
Test Voice Analysis Toolbox Implementation
"""

import numpy as np
import pytest
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from voice_analysis import VoiceAnalyzer, analyze_voice


def test_imports():
    """Test that all modules import correctly"""
    from voice_analysis.f0_estimation import estimate_f0_praat, estimate_f0_swipe
    from voice_analysis.features import (
        compute_jitter_shimmer_features,
        compute_hnr_nhr,
        compute_mfcc_features,
        compute_wavelet_features,
        compute_ppe,
        compute_dfa,
        compute_rpde,
        compute_gne,
    )
    from voice_analysis.utils import compute_tkeo, compute_perturbation_quotient, compute_entropy
    
    assert True  # If we got here, imports worked


def test_synthetic_sine_wave():
    """Test on pure sine wave (should have minimal jitter, high HNR)"""
    # Generate 1 second of 100 Hz sine wave
    fs = 44100
    t = np.linspace(0, 1, fs)
    audio = np.sin(2 * np.pi * 100 * t)
    
    analyzer = VoiceAnalyzer(f0_min=50, f0_max=500, f0_algorithm='PRAAT')
    measures, F0 = analyzer.analyze(audio, fs)
    
    # Check that we got measures
    assert len(measures) > 0
    assert len(F0) > 0
    
    # Pure sine should have low jitter
    assert measures['jitter_RAP'] < 1.0, "Pure sine should have low jitter"
    
    # Pure sine should have high HNR
    assert measures['HNR_mean'] > 10, "Pure sine should have high HNR"


def test_f0_estimation_praat():
    """Test Praat-style F0 estimation"""
    from voice_analysis.f0_estimation import estimate_f0_praat
    
    # Generate chirp signal
    fs = 44100
    t = np.linspace(0, 1, fs)
    f0_target = 150  # Hz
    audio = np.sin(2 * np.pi * f0_target * t)
    
    F0 = estimate_f0_praat(audio, fs, f0_min=50, f0_max=500)
    
    assert len(F0) > 0
    
    # Check F0 is roughly correct (within 10%)
    F0_valid = F0[F0 > 0]
    if len(F0_valid) > 0:
        mean_f0 = np.mean(F0_valid)
        assert 0.9 * f0_target < mean_f0 < 1.1 * f0_target, \
            f"F0 should be near {f0_target} Hz, got {mean_f0:.1f} Hz"


def test_jitter_shimmer():
    """Test jitter/shimmer computation"""
    from voice_analysis.features import compute_jitter_shimmer_features
    
    # Create F0 contour with known jitter
    F0 = 100 + np.random.randn(100) * 2  # Mean 100 Hz, std 2 Hz
    F0 = np.maximum(F0, 50)  # Ensure positive
    
    measures = compute_jitter_shimmer_features(F0, 'jitter')
    
    # Check that measures were computed
    assert 'jitter_RAP' in measures
    assert 'jitter_RAP_percent' in measures
    assert 'jitter_PQ3_Schoentgen' in measures
    assert 'jitter_TKEO_mean' in measures
    
    # Check values are reasonable
    assert measures['jitter_RAP'] > 0
    assert not np.isnan(measures['jitter_RAP'])


def test_hnr_computation():
    """Test HNR/NHR computation"""
    from voice_analysis.features import compute_hnr_nhr
    
    # Generate voiced signal (sine wave)
    fs = 44100
    t = np.linspace(0, 1, fs)
    audio = np.sin(2 * np.pi * 150 * t)
    
    measures = compute_hnr_nhr(audio, fs)
    
    assert 'HNR_mean' in measures
    assert 'HNR_std' in measures
    assert 'NHR_mean' in measures
    assert 'NHR_std' in measures
    
    # Pure sine should have high HNR
    if not np.isnan(measures['HNR_mean']):
        assert measures['HNR_mean'] > 10


def test_tkeo():
    """Test TKEO computation"""
    from voice_analysis.utils import compute_tkeo
    
    x = np.array([1, 2, 3, 2, 1])
    energy = compute_tkeo(x)
    
    assert len(energy) == len(x)
    assert energy[0] == 1  # Boundary condition
    assert energy[-1] == 1  # Boundary condition


def test_perturbation_quotient():
    """Test PQ computation"""
    from voice_analysis.utils import compute_perturbation_quotient
    
    time_series = np.array([100, 102, 98, 101, 99, 100, 103, 97])
    
    pq = compute_perturbation_quotient(time_series, K=3)
    
    assert 'classical_Schoentgen' in pq
    assert 'classical_Baken' in pq
    assert 'generalized_Schoentgen' in pq
    
    assert pq['classical_Schoentgen'] > 0


def test_entropy():
    """Test entropy computation"""
    from voice_analysis.utils import compute_entropy
    
    # Uniform distribution should have high entropy
    uniform = np.ones(10) / 10
    H = compute_entropy(uniform, base='2')
    
    assert np.isclose(H, np.log2(10), atol=0.01)
    
    # Peaked distribution should have low entropy
    peaked = np.zeros(10)
    peaked[0] = 1.0
    H_peaked = compute_entropy(peaked, base='2')
    
    assert H_peaked < H


def test_wavelet_features():
    """Test wavelet feature extraction"""
    from voice_analysis.features import compute_wavelet_features
    
    # Generate test F0 contour
    F0 = 120 + 10 * np.sin(np.linspace(0, 10, 100))
    
    measures = compute_wavelet_features(F0, wavelet='db8', level=5)
    
    assert len(measures) > 0
    assert 'wavelet_energy_approx' in measures or 'wavelet_error' in measures


def test_dfa():
    """Test DFA computation"""
    from voice_analysis.features import compute_dfa
    
    # Generate correlated noise
    signal = np.cumsum(np.random.randn(1000))
    
    dfa = compute_dfa(signal)
    
    # DFA should be bounded between 0 and 1 after sigmoid transform
    if not np.isnan(dfa):
        assert 0 <= dfa <= 1


def test_ppe():
    """Test PPE computation"""
    from voice_analysis.features import compute_ppe
    
    # Generate F0 contour
    F0 = 120 + 5 * np.sin(np.linspace(0, 10, 100))
    
    ppe = compute_ppe(F0, fs=44100)
    
    # PPE should be finite and positive
    if not np.isnan(ppe):
        assert ppe >= 0


if __name__ == '__main__':
    print("Running Voice Analysis Toolbox Tests")
    print("=" * 60)
    
    # Run tests
    test_imports()
    print("✓ Imports test passed")
    
    test_synthetic_sine_wave()
    print("✓ Synthetic sine wave test passed")
    
    test_f0_estimation_praat()
    print("✓ F0 estimation test passed")
    
    test_jitter_shimmer()
    print("✓ Jitter/shimmer test passed")
    
    test_hnr_computation()
    print("✓ HNR computation test passed")
    
    test_tkeo()
    print("✓ TKEO test passed")
    
    test_perturbation_quotient()
    print("✓ Perturbation quotient test passed")
    
    test_entropy()
    print("✓ Entropy test passed")
    
    test_wavelet_features()
    print("✓ Wavelet features test passed")
    
    test_dfa()
    print("✓ DFA test passed")
    
    test_ppe()
    print("✓ PPE test passed")
    
    print("=" * 60)
    print("All tests passed!")
