"""
Validation Script: Compare Python Implementation with MATLAB COVAREP

This script processes the same audio through both Python and MATLAB versions
to validate numerical accuracy and algorithm correctness.
"""

import numpy as np
import matplotlib.pyplot as plt
import soundfile as sf
import sys
import os
from pathlib import Path

# Add covarep to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from covarep.f0 import F0Tracker
from covarep.glottal import iaif


def test_f0_comparison():
    """
    Compare F0 tracking between Python and MATLAB
    
    This requires running the MATLAB version first to get reference output
    """
    print("=" * 70)
    print("F0 TRACKING VALIDATION")
    print("=" * 70)
    
    # Load test audio from COVAREP howtos
    audio_path = "../../howtos/0011.arctic_bdl1.wav"
    
    if not os.path.exists(audio_path):
        print(f"⚠ Audio file not found: {audio_path}")
        print("Using relative path from covarep_python/validation/")
        return None
    
    # Load audio
    try:
        audio, fs = sf.read(audio_path)
        print(f"✓ Loaded audio: {len(audio)} samples at {fs} Hz ({len(audio)/fs:.2f}s)")
    except Exception as e:
        print(f"✗ Error loading audio: {e}")
        return None
    
    # Python F0 tracking
    print("\nRunning Python SRH F0 tracker...")
    tracker = F0Tracker(method='srh', f0_min=50, f0_max=500, hop_size=10.0)
    
    try:
        f0_py, vuv_py, srh_py, times_py = tracker.estimate(audio, fs)
        print(f"✓ Python F0 tracking complete: {len(f0_py)} frames")
        print(f"  F0 range (voiced): {np.min(f0_py[vuv_py]):.1f} - {np.max(f0_py[vuv_py]):.1f} Hz")
        print(f"  Voicing: {100*np.sum(vuv_py)/len(vuv_py):.1f}%")
        
        # Save Python results
        py_data = np.column_stack([times_py, f0_py, vuv_py.astype(int), srh_py])
        np.savetxt('python_f0_output.txt', py_data, fmt='%.6f', delimiter='\t',
                   header='time(s)\tF0(Hz)\tVUV\tSRH_value')
        print(f"  Saved: python_f0_output.txt")
        
    except Exception as e:
        print(f"✗ Python F0 tracking failed: {e}")
        import traceback
        traceback.print_exc()
        return None
    
    # Check if MATLAB reference exists
    matlab_f0_file = "../../howtos/0011.arctic_bdl1.f0.txt"
    
    if os.path.exists(matlab_f0_file):
        print(f"\n✓ Found MATLAB F0 reference: {matlab_f0_file}")
        
        # Load MATLAB F0 (assuming format: time, f0, vuv)
        try:
            matlab_data = np.loadtxt(matlab_f0_file)
            
            if matlab_data.ndim == 2 and matlab_data.shape[1] >= 2:
                times_mat = matlab_data[:, 0]
                f0_mat = matlab_data[:, 1]
                
                print(f"✓ Loaded MATLAB F0: {len(f0_mat)} frames")
                print(f"  F0 range: {np.min(f0_mat[f0_mat>0]):.1f} - {np.max(f0_mat):.1f} Hz")
                
                # Compare (need to align time scales)
                if len(f0_py) == len(f0_mat):
                    # Direct comparison
                    voiced_mask = (f0_mat > 0) & (f0_py > 0)
                    if np.sum(voiced_mask) > 0:
                        error = np.abs(f0_py[voiced_mask] - f0_mat[voiced_mask])
                        mean_error = np.mean(error)
                        median_error = np.median(error)
                        max_error = np.max(error)
                        
                        print(f"\n📊 Comparison (voiced frames):")
                        print(f"  Mean error: {mean_error:.2f} Hz")
                        print(f"  Median error: {median_error:.2f} Hz")
                        print(f"  Max error: {max_error:.2f} Hz")
                        print(f"  Relative error: {100*mean_error/np.mean(f0_mat[f0_mat>0]):.1f}%")
                        
                        # Plot comparison
                        plot_f0_comparison(times_py, f0_py, vuv_py, times_mat, f0_mat)
                else:
                    print(f"⚠ Length mismatch: Python={len(f0_py)}, MATLAB={len(f0_mat)}")
                    print("  (May need to adjust hop size or alignment)")
                    
                    # Plot anyway for visual inspection
                    plot_f0_comparison(times_py, f0_py, vuv_py, times_mat, f0_mat)
        
        except Exception as e:
            print(f"✗ Error loading MATLAB reference: {e}")
    else:
        print(f"\n⚠ MATLAB F0 reference not found: {matlab_f0_file}")
        print("Run MATLAB COVAREP first or provide reference data")
        
        # Plot Python results only
        plt.figure(figsize=(12, 6))
        
        t_audio = np.arange(len(audio)) / fs
        plt.subplot(2, 1, 1)
        plt.plot(t_audio, audio, 'k', linewidth=0.5, alpha=0.7)
        plt.ylabel('Amplitude')
        plt.title('Speech Signal')
        plt.grid(True, alpha=0.3)
        plt.xlim([0, len(audio)/fs])
        
        plt.subplot(2, 1, 2)
        plt.plot(times_py, f0_py, 'b.-', linewidth=1.5, markersize=3, label='Python F0')
        plt.fill_between(times_py, 0, 500, where=~vuv_py, alpha=0.2, color='red', label='Unvoiced')
        plt.xlabel('Time (s)')
        plt.ylabel('F0 (Hz)')
        plt.title('F0 Contour - Python Implementation')
        plt.ylim([0, 500])
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.xlim([0, len(audio)/fs])
        
        plt.tight_layout()
        plt.savefig('python_f0_only.png', dpi=150)
        print("✓ Saved plot: python_f0_only.png")
    
    return f0_py, vuv_py, times_py


def plot_f0_comparison(times_py, f0_py, vuv_py, times_mat, f0_mat):
    """
    Plot Python vs MATLAB F0 comparison
    """
    fig, axes = plt.subplots(2, 1, figsize=(14, 8))
    
    # Plot both F0 contours
    axes[0].plot(times_mat, f0_mat, 'r.-', linewidth=1.5, markersize=2, 
                 alpha=0.7, label='MATLAB')
    axes[0].plot(times_py, f0_py, 'b.-', linewidth=1.5, markersize=2, 
                 alpha=0.7, label='Python')
    axes[0].fill_between(times_py, 0, 500, where=~vuv_py, alpha=0.15, 
                         color='gray', label='Python Unvoiced')
    axes[0].set_ylabel('F0 (Hz)')
    axes[0].set_title('F0 Comparison: Python vs MATLAB')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    axes[0].set_ylim([0, 500])
    
    # Plot difference (where both are voiced)
    if len(f0_py) == len(f0_mat):
        voiced_mask = (f0_mat > 0) & (f0_py > 0)
        diff = f0_py - f0_mat
        
        axes[1].plot(times_py[voiced_mask], diff[voiced_mask], 'g.-', 
                     linewidth=1, markersize=2, alpha=0.7)
        axes[1].axhline(y=0, color='k', linestyle='--', linewidth=1)
        axes[1].set_xlabel('Time (s)')
        axes[1].set_ylabel('Error (Hz)')
        axes[1].set_title('F0 Error (Python - MATLAB) on Voiced Frames')
        axes[1].grid(True, alpha=0.3)
    else:
        axes[1].text(0.5, 0.5, 'Length mismatch - cannot compute frame-by-frame error', 
                     ha='center', va='center', transform=axes[1].transAxes)
    
    plt.tight_layout()
    plt.savefig('f0_comparison_python_vs_matlab.png', dpi=150)
    print("✓ Saved comparison plot: f0_comparison_python_vs_matlab.png")


def test_iaif_comparison():
    """
    Compare IAIF glottal analysis between Python and MATLAB
    """
    print("\n" + "=" * 70)
    print("IAIF VALIDATION")
    print("=" * 70)
    
    # Load test audio
    audio_path = "../../howtos/0011.arctic_bdl1.wav"
    
    if not os.path.exists(audio_path):
        print(f"⚠ Audio file not found: {audio_path}")
        return None
    
    try:
        audio, fs = sf.read(audio_path)
        print(f"✓ Loaded audio: {len(audio)} samples at {fs} Hz")
    except Exception as e:
        print(f"✗ Error loading audio: {e}")
        return None
    
    # Extract a voiced frame (30 ms around 1 second)
    start_sample = int(1.0 * fs)
    frame_len = int(0.03 * fs)
    
    if start_sample + frame_len > len(audio):
        start_sample = 0
    
    frame = audio[start_sample:start_sample + frame_len]
    print(f"✓ Extracted frame: {len(frame)} samples ({len(frame)/fs*1000:.1f} ms)")
    
    # Python IAIF
    print("\nRunning Python IAIF...")
    try:
        g_py, dg_py, a_py, ag_py = iaif(frame, fs)
        print(f"✓ Python IAIF complete")
        print(f"  Glottal flow: {len(g_py)} samples")
        print(f"  VT order: {len(a_py)-1}")
        print(f"  GL order: {len(ag_py)-1}")
        print(f"  Glottal flow range: [{np.min(g_py):.4f}, {np.max(g_py):.4f}]")
        print(f"  Flow derivative range: [{np.min(dg_py):.4f}, {np.max(dg_py):.4f}]")
        
        # Save results
        np.savetxt('python_iaif_glottal_flow.txt', g_py, fmt='%.8f')
        np.savetxt('python_iaif_flow_derivative.txt', dg_py, fmt='%.8f')
        print(f"  Saved: python_iaif_glottal_flow.txt")
        print(f"  Saved: python_iaif_flow_derivative.txt")
        
    except Exception as e:
        print(f"✗ Python IAIF failed: {e}")
        import traceback
        traceback.print_exc()
        return None
    
    # Plot results
    fig, axes = plt.subplots(3, 1, figsize=(12, 9))
    
    t_ms = np.arange(len(frame)) / fs * 1000
    
    axes[0].plot(t_ms, frame, 'k', linewidth=1)
    axes[0].set_ylabel('Amplitude')
    axes[0].set_title('Original Speech Frame')
    axes[0].grid(True, alpha=0.3)
    
    t_g_ms = np.arange(len(g_py)) / fs * 1000
    axes[1].plot(t_g_ms, g_py, 'b', linewidth=1)
    axes[1].set_ylabel('Amplitude')
    axes[1].set_title('Estimated Glottal Flow (Python IAIF)')
    axes[1].grid(True, alpha=0.3)
    
    axes[2].plot(t_g_ms, dg_py, 'r', linewidth=1)
    axes[2].set_xlabel('Time (ms)')
    axes[2].set_ylabel('Amplitude')
    axes[2].set_title('Glottal Flow Derivative')
    axes[2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('iaif_python_output.png', dpi=150)
    print("✓ Saved IAIF plot: iaif_python_output.png")
    
    print("\n💡 To compare with MATLAB:")
    print("   1. Run HOWTO_glottalsource.m in MATLAB")
    print("   2. Extract glottal flow waveform")
    print("   3. Compare visually or save MATLAB output for numerical comparison")
    
    return g_py, dg_py


def generate_matlab_test_script():
    """
    Generate MATLAB script to produce reference outputs
    """
    matlab_script = """% MATLAB Reference Output Generator for Python Validation
% Run this in MATLAB COVAREP to generate reference data

% Setup COVAREP path
run('../startup.m');

% Load test audio
[wave, fs] = audioread('0011.arctic_bdl1.wav');

%% F0 Tracking with SRH
fprintf('Computing F0 with SRH...\\n');
f0min = 50;
f0max = 500;
hopsize = 10; % ms

[F0s, VUVDecisions, SRHVal, time] = pitch_srh(wave, fs, f0min, f0max, hopsize);

% Save F0 results
f0_data = [time(:), F0s(:), double(VUVDecisions(:)), SRHVal(:)];
save('matlab_f0_reference.mat', 'F0s', 'VUVDecisions', 'SRHVal', 'time');
dlmwrite('matlab_f0_reference.txt', f0_data, 'delimiter', '\\t', 'precision', 6);
fprintf('✓ Saved MATLAB F0 reference\\n');

%% IAIF on sample frame
fprintf('\\nComputing IAIF on sample frame...\\n');
% Extract 30ms frame at 1 second
start_sample = round(1.0 * fs);
frame_len = round(0.03 * fs);
frame = wave(start_sample:start_sample+frame_len-1);

% Run IAIF
[g, dg, a, ag] = iaif(frame, fs);

% Save IAIF results
iaif_data = struct('g', g, 'dg', dg, 'a', a, 'ag', ag, 'fs', fs);
save('matlab_iaif_reference.mat', '-struct', 'iaif_data');

% Save as text too
dlmwrite('matlab_iaif_glottal_flow.txt', g, 'delimiter', '\\t', 'precision', 8);
dlmwrite('matlab_iaif_flow_derivative.txt', dg, 'delimiter', '\\t', 'precision', 8);

fprintf('✓ Saved MATLAB IAIF reference\\n');
fprintf('\\n');
fprintf('Reference files created:\\n');
fprintf('  - matlab_f0_reference.mat/.txt\\n');
fprintf('  - matlab_iaif_reference.mat\\n');
fprintf('  - matlab_iaif_glottal_flow.txt\\n');
fprintf('  - matlab_iaif_flow_derivative.txt\\n');
"""
    
    output_path = "generate_matlab_reference.m"
    with open(output_path, 'w') as f:
        f.write(matlab_script)
    
    print(f"\n✓ Generated MATLAB script: {output_path}")
    print("\nTo generate reference data:")
    print("  1. Copy validation/generate_matlab_reference.m to covarep/howtos/")
    print("  2. Run in MATLAB: generate_matlab_reference")
    print("  3. Run this Python script again to compare")


def run_validation():
    """
    Run complete validation suite
    """
    print("╔" + "═" * 68 + "╗")
    print("║" + " " * 15 + "COVAREP PYTHON - VALIDATION PHASE" + " " * 20 + "║")
    print("╚" + "═" * 68 + "╝")
    
    print("\nThis validation compares Python implementation with MATLAB COVAREP")
    print("to ensure numerical accuracy and algorithm correctness.\n")
    
    # Check if we're in the right directory
    if not os.path.exists("../../howtos"):
        print("⚠ Warning: Not in covarep_python/validation/ directory")
        print("Please run from: covarep_python/validation/")
        # Try anyway
    
    # Test F0 tracking
    result_f0 = test_f0_comparison()
    
    # Test IAIF
    result_iaif = test_iaif_comparison()
    
    # Generate MATLAB reference script
    generate_matlab_test_script()
    
    print("\n" + "=" * 70)
    print("VALIDATION SUMMARY")
    print("=" * 70)
    
    if result_f0 is not None:
        print("✓ F0 tracking: Python implementation runs successfully")
    else:
        print("⚠ F0 tracking: Issues encountered")
    
    if result_iaif is not None:
        print("✓ IAIF: Python implementation runs successfully")
    else:
        print("⚠ IAIF: Issues encountered")
    
    print("\n📋 Next Steps:")
    print("  1. Run generate_matlab_reference.m in MATLAB COVAREP")
    print("  2. Re-run this validation script to compare outputs")
    print("  3. Analyze differences and tune parameters")
    print("  4. Iterate until error < 10%")
    
    print("\n" + "=" * 70)


if __name__ == "__main__":
    run_validation()
