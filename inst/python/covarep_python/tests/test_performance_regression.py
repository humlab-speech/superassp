"""
Performance Regression Tests

Ensures that optimizations maintain their performance over time.
Tests will fail if performance degrades significantly.

Usage:
    pytest tests/test_performance_regression.py
    pytest tests/test_performance_regression.py --benchmark-only
"""

import pytest
import numpy as np
import sys
from pathlib import Path

# Add parent to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from covarep.f0 import pitch_srh
from covarep.f0.f0_optimized import pitch_srh_vectorized
from covarep.glottal import iaif
from covarep.glottal.iaif_optimized import iaif_optimized


# Performance targets (based on benchmark results)
F0_TARGET_MS = 250  # Max 250ms for 3.5s audio
IAIF_TARGET_MS = 1.0  # Max 1ms per frame
F0_SPEEDUP_TARGET = 5.0  # Minimum 5x speedup
IAIF_SPEEDUP_TARGET = 1.3  # Minimum 1.3x speedup


@pytest.fixture
def test_audio():
    """Generate test audio signal"""
    fs = 16000
    duration = 3.5
    t = np.linspace(0, duration, int(duration * fs))

    # Generate signal with varying F0
    f0_base = 150
    f0_var = 20 * np.sin(2 * np.pi * 3 * t)
    phase = np.cumsum(2 * np.pi * (f0_base + f0_var) / fs)

    audio = np.sin(phase)

    # Add some noise
    audio += 0.1 * np.random.randn(len(audio))

    return audio, fs


@pytest.fixture
def test_frame():
    """Generate test frame for IAIF"""
    fs = 16000
    duration = 0.03
    t = np.linspace(0, duration, int(duration * fs))

    # Simple periodic signal
    frame = np.sin(2 * np.pi * 150 * t)

    return frame, fs


# F0 Tracking Tests
class TestF0Performance:
    """F0 tracking performance regression tests"""

    def test_f0_optimized_speed(self, test_audio, benchmark):
        """Test that optimized F0 tracking is fast enough"""
        audio, fs = test_audio

        result = benchmark(pitch_srh_vectorized, audio, fs)

        # Check that we meet performance target
        assert benchmark.stats['mean'] < F0_TARGET_MS / 1000, \
            f"F0 tracking too slow: {benchmark.stats['mean']*1000:.1f}ms > {F0_TARGET_MS}ms"

    def test_f0_speedup(self, test_audio):
        """Test that optimization provides expected speedup"""
        import time

        audio, fs = test_audio

        # Measure original
        times_orig = []
        for _ in range(3):
            start = time.perf_counter()
            _ = pitch_srh(audio, fs)
            times_orig.append(time.perf_counter() - start)

        # Measure optimized
        times_opt = []
        for _ in range(3):
            start = time.perf_counter()
            _ = pitch_srh_vectorized(audio, fs)
            times_opt.append(time.perf_counter() - start)

        speedup = np.mean(times_orig) / np.mean(times_opt)

        assert speedup >= F0_SPEEDUP_TARGET, \
            f"Speedup too low: {speedup:.2f}x < {F0_SPEEDUP_TARGET}x"

    def test_f0_accuracy(self, test_audio):
        """Test that optimization maintains accuracy"""
        audio, fs = test_audio

        f0_orig, vuv_orig, srh_orig, times_orig = pitch_srh(audio, fs)
        f0_opt, vuv_opt, srh_opt, times_opt = pitch_srh_vectorized(audio, fs)

        # Check results match
        assert len(f0_orig) == len(f0_opt), "Frame count mismatch"

        # Check F0 values on voiced frames
        both_voiced = vuv_orig & vuv_opt
        if np.sum(both_voiced) > 0:
            f0_diff = np.abs(f0_orig[both_voiced] - f0_opt[both_voiced])
            median_error = np.median(f0_diff)

            assert median_error < 5.0, \
                f"F0 median error too large: {median_error:.2f} Hz"

        # Check VUV agreement
        vuv_agreement = np.sum(vuv_orig == vuv_opt) / len(vuv_orig)
        assert vuv_agreement > 0.75, \
            f"VUV agreement too low: {vuv_agreement*100:.1f}%"


# IAIF Tests
class TestIAIFPerformance:
    """IAIF performance regression tests"""

    def test_iaif_optimized_speed(self, test_frame, benchmark):
        """Test that optimized IAIF is fast enough"""
        frame, fs = test_frame

        result = benchmark(iaif_optimized, frame, fs)

        # Check that we meet performance target
        assert benchmark.stats['mean'] < IAIF_TARGET_MS / 1000, \
            f"IAIF too slow: {benchmark.stats['mean']*1000:.3f}ms > {IAIF_TARGET_MS}ms"

    def test_iaif_speedup(self, test_frame):
        """Test that optimization provides expected speedup"""
        import time

        frame, fs = test_frame

        # Measure original
        times_orig = []
        for _ in range(10):
            start = time.perf_counter()
            _ = iaif(frame, fs)
            times_orig.append(time.perf_counter() - start)

        # Measure optimized
        times_opt = []
        for _ in range(10):
            start = time.perf_counter()
            _ = iaif_optimized(frame, fs)
            times_opt.append(time.perf_counter() - start)

        speedup = np.mean(times_orig) / np.mean(times_opt)

        assert speedup >= IAIF_SPEEDUP_TARGET, \
            f"Speedup too low: {speedup:.2f}x < {IAIF_SPEEDUP_TARGET}x"

    def test_iaif_accuracy(self, test_frame):
        """Test that optimization maintains accuracy"""
        frame, fs = test_frame

        g_orig, dg_orig, a_orig, ag_orig = iaif(frame, fs)
        g_opt, dg_opt, a_opt, ag_opt = iaif_optimized(frame, fs)

        # Check lengths match
        assert len(g_orig) == len(g_opt), "Glottal flow length mismatch"

        # Check numerical accuracy
        g_diff = np.max(np.abs(g_orig - g_opt))
        dg_diff = np.max(np.abs(dg_orig - dg_opt))

        assert g_diff < 0.1, \
            f"Glottal flow difference too large: {g_diff:.6f}"
        assert dg_diff < 0.1, \
            f"Flow derivative difference too large: {dg_diff:.6f}"


# Memory Tests
class TestMemoryUsage:
    """Memory usage regression tests"""

    @pytest.mark.slow
    def test_f0_memory_reasonable(self, test_audio):
        """Test that F0 tracking doesn't use excessive memory"""
        try:
            import tracemalloc
        except ImportError:
            pytest.skip("tracemalloc not available")

        audio, fs = test_audio

        tracemalloc.start()

        # Run F0 tracking
        _ = pitch_srh_vectorized(audio, fs)

        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()

        # Peak memory should be reasonable (< 50 MB for 3.5s audio)
        peak_mb = peak / 1024 / 1024

        assert peak_mb < 50, \
            f"Memory usage too high: {peak_mb:.1f} MB"


# Correctness Tests
class TestCorrectnessInvariants:
    """Test that optimizations preserve correctness"""

    def test_f0_monotonic_time(self, test_audio):
        """Test that time stamps are monotonically increasing"""
        audio, fs = test_audio

        _, _, _, times = pitch_srh_vectorized(audio, fs)

        # Check monotonicity
        assert np.all(np.diff(times) > 0), "Time stamps not monotonic"

    def test_f0_values_in_range(self, test_audio):
        """Test that F0 values are within specified range"""
        audio, fs = test_audio

        f0, vuv, _, _ = pitch_srh_vectorized(audio, fs, f0min=50, f0max=500)

        # Check that voiced F0 values are in range
        voiced_f0 = f0[vuv]

        if len(voiced_f0) > 0:
            assert np.all(voiced_f0 >= 50), "F0 below minimum"
            assert np.all(voiced_f0 <= 500), "F0 above maximum"

    def test_iaif_output_shapes(self, test_frame):
        """Test that IAIF output shapes are correct"""
        frame, fs = test_frame

        g, dg, a, ag = iaif_optimized(frame, fs)

        # Check that outputs have reasonable lengths
        assert len(g) > 0, "Empty glottal flow"
        assert len(dg) > 0, "Empty flow derivative"
        assert len(a) > 0, "Empty VT coefficients"
        assert len(ag) > 0, "Empty GL coefficients"

        # Check coefficient arrays start with 1.0
        assert a[0] == 1.0, "VT coefficients don't start with 1.0"
        assert ag[0] == 1.0, "GL coefficients don't start with 1.0"


# Parametric Tests
@pytest.mark.parametrize("duration", [1.0, 3.5, 5.0])
def test_f0_scales_linearly(duration):
    """Test that F0 tracking scales linearly with audio length"""
    import time

    fs = 16000
    t = np.linspace(0, duration, int(duration * fs))
    audio = np.sin(2 * np.pi * 150 * t)

    # Measure time
    start = time.perf_counter()
    _ = pitch_srh_vectorized(audio, fs)
    elapsed = time.perf_counter() - start

    # Check RTF (should be < 0.1 for all durations)
    rtf = elapsed / duration
    assert rtf < 0.1, f"RTF too high for {duration}s audio: {rtf:.4f}"


@pytest.mark.parametrize("fs", [8000, 16000, 32000])
def test_f0_handles_sample_rates(fs):
    """Test that F0 tracking handles different sample rates"""
    duration = 1.0
    t = np.linspace(0, duration, int(duration * fs))
    audio = np.sin(2 * np.pi * 150 * t)

    # Should complete without error
    f0, vuv, srh, times = pitch_srh_vectorized(audio, fs)

    assert len(f0) > 0, f"No F0 estimates for fs={fs}"


if __name__ == "__main__":
    # Run tests with pytest
    pytest.main([__file__, "-v", "--tb=short"])
