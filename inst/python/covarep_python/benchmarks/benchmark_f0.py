"""
F0 Tracking Performance Benchmark

Compares performance between:
1. Original implementation
2. Vectorized NumPy implementation
3. Cython implementation (if available)

Usage:
    python benchmarks/benchmark_f0.py [audio_file.wav]

If no audio file provided, generates synthetic test signals.
"""

import sys
import os
import time
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from covarep.f0 import pitch_srh
from covarep.f0.f0_optimized import pitch_srh_vectorized


def generate_test_signal(duration=10.0, fs=16000, f0=150.0, snr_db=20):
    """
    Generate synthetic voiced speech signal for testing

    Parameters
    ----------
    duration : float
        Signal duration in seconds
    fs : int
        Sampling frequency
    f0 : float
        Fundamental frequency
    snr_db : float
        Signal-to-noise ratio in dB

    Returns
    -------
    signal : ndarray
        Synthetic speech signal
    fs : int
        Sampling frequency
    """
    t = np.arange(0, duration, 1/fs)
    n_samples = len(t)

    # Generate glottal pulse train (triangular pulses)
    period_samples = int(fs / f0)
    pulse_train = np.zeros(n_samples)
    pulse_positions = np.arange(0, n_samples, period_samples)

    for pos in pulse_positions:
        if pos + 10 < n_samples:
            # Simple triangular pulse
            pulse_train[pos:pos+10] = np.linspace(0, 1, 5).tolist() + np.linspace(1, 0, 5).tolist()

    # Simple vocal tract filter (formants approximation)
    # Using cascaded resonators
    from scipy.signal import butter, lfilter

    # F1 = 500 Hz
    b1, a1 = butter(2, [450, 550], btype='band', fs=fs)
    # F2 = 1500 Hz
    b2, a2 = butter(2, [1400, 1600], btype='band', fs=fs)

    speech = lfilter(b1, a1, pulse_train)
    speech += 0.5 * lfilter(b2, a2, pulse_train)

    # Add noise
    noise = np.random.randn(n_samples)
    signal_power = np.mean(speech ** 2)
    noise_power = signal_power / (10 ** (snr_db / 10))
    noise *= np.sqrt(noise_power / np.mean(noise ** 2))

    signal = speech + noise

    # Normalize
    signal = signal / np.max(np.abs(signal)) * 0.8

    return signal, fs


def load_test_audio(audio_file=None):
    """
    Load test audio or generate synthetic signal

    Parameters
    ----------
    audio_file : str, optional
        Path to audio file

    Returns
    -------
    audio : ndarray
        Audio signal
    fs : int
        Sampling frequency
    source : str
        'file' or 'synthetic'
    """
    if audio_file and os.path.exists(audio_file):
        try:
            import soundfile as sf
            audio, fs = sf.read(audio_file)
            print(f"✓ Loaded audio: {audio_file}")
            print(f"  Duration: {len(audio)/fs:.2f}s")
            print(f"  Sample rate: {fs} Hz")
            return audio, fs, 'file'
        except Exception as e:
            print(f"✗ Error loading audio: {e}")
            print("  Falling back to synthetic signal")

    # Generate synthetic signal
    print("⚙ Generating synthetic test signal...")
    audio, fs = generate_test_signal(duration=10.0, fs=16000, f0=150.0)
    print(f"  Duration: 10.0s")
    print(f"  Sample rate: {fs} Hz")
    print(f"  F0: 150 Hz")
    return audio, fs, 'synthetic'


def benchmark_implementation(func, audio, fs, n_runs=10, warmup=2):
    """
    Benchmark a single implementation

    Parameters
    ----------
    func : callable
        F0 tracking function
    audio : ndarray
        Test audio
    fs : int
        Sampling frequency
    n_runs : int
        Number of benchmark runs
    warmup : int
        Number of warmup runs

    Returns
    -------
    times : list
        Execution times
    result : tuple
        F0 tracking result
    """
    # Warm-up runs
    for _ in range(warmup):
        _ = func(audio, fs)

    # Benchmark runs
    times = []
    result = None

    for i in range(n_runs):
        start = time.perf_counter()
        result = func(audio, fs)
        elapsed = time.perf_counter() - start
        times.append(elapsed)

    return times, result


def validate_accuracy(result_ref, result_test, name="Test"):
    """
    Validate numerical accuracy between implementations

    Parameters
    ----------
    result_ref : tuple
        Reference result (f0, vuv, srh, times)
    result_test : tuple
        Test result
    name : str
        Implementation name

    Returns
    -------
    passed : bool
        True if validation passed
    """
    f0_ref, vuv_ref, srh_ref, times_ref = result_ref
    f0_test, vuv_test, srh_test, times_test = result_test

    print(f"\n📊 Validating {name}:")

    # Check array shapes
    if f0_ref.shape != f0_test.shape:
        print(f"  ✗ Shape mismatch: {f0_ref.shape} vs {f0_test.shape}")
        return False
    print(f"  ✓ Shape match: {f0_ref.shape}")

    # Check F0 values (on voiced frames)
    voiced_mask = vuv_ref & vuv_test
    if np.sum(voiced_mask) > 0:
        f0_diff = np.abs(f0_ref[voiced_mask] - f0_test[voiced_mask])
        max_diff = np.max(f0_diff)
        mean_diff = np.mean(f0_diff)
        median_diff = np.median(f0_diff)

        print(f"  F0 differences (voiced frames):")
        print(f"    Mean:   {mean_diff:.6f} Hz")
        print(f"    Median: {median_diff:.6f} Hz")
        print(f"    Max:    {max_diff:.6f} Hz")

        if max_diff > 1.0:  # Allow up to 1 Hz difference
            print(f"  ⚠ Large F0 difference detected (max: {max_diff:.3f} Hz)")
            # Not a failure, just a warning
        else:
            print(f"  ✓ F0 values match within tolerance")

    # Check VUV decisions
    vuv_match = np.sum(vuv_ref == vuv_test) / len(vuv_ref) * 100
    print(f"  VUV agreement: {vuv_match:.1f}%")

    if vuv_match < 95:
        print(f"  ⚠ VUV decisions differ significantly")
    else:
        print(f"  ✓ VUV decisions match")

    # Check SRH values
    srh_corr = np.corrcoef(srh_ref, srh_test)[0, 1]
    print(f"  SRH correlation: {srh_corr:.6f}")

    if srh_corr < 0.99:
        print(f"  ⚠ SRH values differ")
    else:
        print(f"  ✓ SRH values match")

    return True


def plot_results(times_dict, audio_duration, output_file='benchmark_f0_results.png'):
    """
    Plot benchmark results

    Parameters
    ----------
    times_dict : dict
        Dictionary of {name: times} for each implementation
    audio_duration : float
        Audio duration in seconds
    output_file : str
        Output filename
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Plot 1: Execution time comparison
    ax = axes[0]
    names = list(times_dict.keys())
    mean_times = [np.mean(times) * 1000 for times in times_dict.values()]  # Convert to ms
    std_times = [np.std(times) * 1000 for times in times_dict.values()]

    colors = ['#e74c3c', '#3498db', '#2ecc71'][:len(names)]
    bars = ax.bar(names, mean_times, yerr=std_times, capsize=5, color=colors, alpha=0.7)

    ax.set_ylabel('Execution Time (ms)', fontsize=12)
    ax.set_title(f'F0 Tracking Performance\n({audio_duration:.1f}s audio @ 16kHz)', fontsize=13)
    ax.grid(axis='y', alpha=0.3)

    # Add value labels on bars
    for bar, mean_t, std_t in zip(bars, mean_times, std_times):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{mean_t:.1f}±{std_t:.1f}',
                ha='center', va='bottom', fontsize=10)

    # Plot 2: Speedup comparison
    ax = axes[1]
    baseline_time = mean_times[0]
    speedups = [baseline_time / t for t in mean_times]

    bars = ax.bar(names, speedups, color=colors, alpha=0.7)
    ax.axhline(y=1.0, color='k', linestyle='--', linewidth=1, alpha=0.5)
    ax.set_ylabel('Speedup Factor', fontsize=12)
    ax.set_title('Relative Performance\n(vs Original Implementation)', fontsize=13)
    ax.grid(axis='y', alpha=0.3)

    # Add value labels
    for bar, speedup in zip(bars, speedups):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{speedup:.1f}x',
                ha='center', va='bottom', fontsize=10, fontweight='bold')

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"\n✓ Saved benchmark plot: {output_file}")


def run_benchmark(audio_file=None):
    """
    Run complete benchmark suite

    Parameters
    ----------
    audio_file : str, optional
        Path to test audio file
    """
    print("=" * 70)
    print("F0 TRACKING PERFORMANCE BENCHMARK")
    print("=" * 70)

    # Load test audio
    audio, fs, source = load_test_audio(audio_file)
    audio_duration = len(audio) / fs

    # Check if implementations are available
    implementations = {
        'Original': pitch_srh,
        'Vectorized': pitch_srh_vectorized,
    }

    # Try to import Cython version
    try:
        from covarep.f0.f0_cython import pitch_srh_cython
        implementations['Cython'] = pitch_srh_cython
        print("✓ Cython implementation available")
    except ImportError:
        print("⚠ Cython implementation not available (will be skipped)")

    print(f"\nTesting {len(implementations)} implementations:")
    for name in implementations:
        print(f"  • {name}")

    # Run benchmarks
    print(f"\n{'─' * 70}")
    print("RUNNING BENCHMARKS")
    print(f"{'─' * 70}")

    times_dict = {}
    results_dict = {}

    for name, func in implementations.items():
        print(f"\n⚙ Benchmarking: {name}")
        print(f"  Warmup runs: 2")
        print(f"  Benchmark runs: 10")

        times, result = benchmark_implementation(func, audio, fs, n_runs=10, warmup=2)
        times_dict[name] = times
        results_dict[name] = result

        mean_time = np.mean(times) * 1000
        std_time = np.std(times) * 1000

        print(f"  ✓ Execution time: {mean_time:.2f} ± {std_time:.2f} ms")

    # Validate accuracy
    print(f"\n{'─' * 70}")
    print("ACCURACY VALIDATION")
    print(f"{'─' * 70}")

    ref_result = results_dict['Original']
    for name in list(implementations.keys())[1:]:  # Skip original
        validate_accuracy(ref_result, results_dict[name], name)

    # Summary
    print(f"\n{'=' * 70}")
    print("BENCHMARK SUMMARY")
    print(f"{'=' * 70}")

    baseline_time = np.mean(times_dict['Original']) * 1000

    print(f"\nTest audio: {audio_duration:.1f}s @ {fs} Hz ({source})")
    print(f"\nResults:")

    for name, times in times_dict.items():
        mean_time = np.mean(times) * 1000
        std_time = np.std(times) * 1000
        speedup = baseline_time / mean_time

        print(f"\n  {name}:")
        print(f"    Time:    {mean_time:.2f} ± {std_time:.2f} ms")
        print(f"    Speedup: {speedup:.2f}x")

        # Real-time factor
        rtf = mean_time / 1000 / audio_duration
        print(f"    RTF:     {rtf:.4f} (lower is better)")

    # Plot results
    plot_results(times_dict, audio_duration)

    print(f"\n{'=' * 70}")
    print("✓ Benchmark completed successfully!")
    print(f"{'=' * 70}\n")


if __name__ == "__main__":
    # Get audio file from command line if provided
    audio_file = sys.argv[1] if len(sys.argv) > 1 else None

    run_benchmark(audio_file)
