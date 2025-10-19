#!/usr/bin/env python3
"""
Quick test of IAIF fixes - validates improved implementation
"""

import numpy as np
import soundfile as sf
import matplotlib.pyplot as plt
from covarep.glottal import iaif
from covarep.f0 import F0Tracker

# Load test audio
audio_path = '../howtos/arctic_a0007.wav'
audio, fs = sf.read(audio_path)

print("="*70)
print("IAIF FIX VALIDATION TEST")
print("="*70)
print(f"Audio: {len(audio)} samples @ {fs} Hz ({len(audio)/fs:.2f}s)")

# Get F0 for a voiced frame
tracker = F0Tracker(method='srh', f0_min=50, f0_max=500)
f0, vuv, srh, times = tracker.estimate(audio, fs)
print(f"F0 estimated: {np.median(f0[vuv>0]):.1f} Hz (median)")

# Take a voiced frame (where F0 is detected)
voiced_frames = np.where(vuv > 0)[0]
if len(voiced_frames) > 10:
    # Take a frame from middle
    frame_idx = voiced_frames[len(voiced_frames)//2]
    frame_time = times[frame_idx]
    
    # Extract 30ms frame
    frame_start = int((frame_time - 0.015) * fs)
    frame_end = int((frame_time + 0.015) * fs)
    frame_start = max(0, frame_start)
    frame_end = min(len(audio), frame_end)
    
    frame = audio[frame_start:frame_end]
    
    print(f"\nTesting on voiced frame:")
    print(f"  Time: {frame_time:.3f}s")
    print(f"  F0: {f0[frame_idx]:.1f} Hz")
    print(f"  Frame length: {len(frame)} samples")
    
    # Run IAIF
    g, dg, a, ag = iaif(frame, fs)
    
    print(f"\nIAIF Results:")
    print(f"  Glottal flow (g):")
    print(f"    - Length: {len(g)} samples")
    print(f"    - Range: [{g.min():.6f}, {g.max():.6f}]")
    print(f"    - Mean: {g.mean():.6f}")
    print(f"    - Std: {g.std():.6f}")
    
    print(f"  Glottal derivative (dg):")
    print(f"    - Length: {len(dg)} samples")
    print(f"    - Range: [{dg.min():.6f}, {dg.max():.6f}]")
    
    print(f"  Vocal tract filter (a):")
    print(f"    - Order: {len(a)-1}")
    print(f"    - First 5 coeffs: {a[:5]}")
    
    print(f"  Glottal source filter (ag):")
    print(f"    - Order: {len(ag)-1}")
    print(f"    - Coeffs: {ag}")
    
    # Visual test - plot results
    fig, axes = plt.subplots(3, 1, figsize=(12, 8))
    
    # Plot 1: Original speech frame
    t_frame = np.arange(len(frame)) / fs * 1000  # ms
    axes[0].plot(t_frame, frame, 'b-', linewidth=0.8)
    axes[0].set_title(f'Original Speech Frame (F0={f0[frame_idx]:.1f} Hz)', fontsize=12)
    axes[0].set_ylabel('Amplitude')
    axes[0].grid(True, alpha=0.3)
    axes[0].set_xlim([0, len(frame)/fs*1000])
    
    # Plot 2: Glottal flow
    t_g = np.arange(len(g)) / fs * 1000
    axes[1].plot(t_g, g, 'r-', linewidth=0.8)
    axes[1].set_title('Estimated Glottal Flow', fontsize=12)
    axes[1].set_ylabel('Flow')
    axes[1].grid(True, alpha=0.3)
    axes[1].set_xlim([0, len(g)/fs*1000])
    
    # Plot 3: Glottal flow derivative
    t_dg = np.arange(len(dg)) / fs * 1000
    axes[2].plot(t_dg, dg, 'g-', linewidth=0.8)
    axes[2].set_title('Glottal Flow Derivative', fontsize=12)
    axes[2].set_xlabel('Time (ms)')
    axes[2].set_ylabel('Derivative')
    axes[2].grid(True, alpha=0.3)
    axes[2].set_xlim([0, len(dg)/fs*1000])
    
    plt.tight_layout()
    plt.savefig('validation/iaif_test_output.png', dpi=150, bbox_inches='tight')
    print(f"\n✓ Plot saved to: validation/iaif_test_output.png")
    
    # Save output for MATLAB comparison
    np.savetxt('validation/python_iaif_glottal_flow_fixed.txt', g)
    np.savetxt('validation/python_iaif_flow_derivative_fixed.txt', dg)
    np.savetxt('validation/python_iaif_vt_coeffs.txt', a)
    np.savetxt('validation/python_iaif_gl_coeffs.txt', ag)
    
    print(f"✓ Outputs saved for MATLAB comparison")
    
    print("\n" + "="*70)
    print("✓ IAIF TEST COMPLETE")
    print("="*70)
    print("\nNext steps:")
    print("  1. Compare with MATLAB reference (if available)")
    print("  2. Check glottal flow shape - should show periodic pulses")
    print("  3. Verify VT order matches expected: 2*round(fs/2000)+4")
    print(f"     Expected for {fs}Hz: {2*round(fs/2000)+4}, Got: {len(a)-1}")
    print("  4. Verify GL order: 2*round(fs/4000)")
    print(f"     Expected for {fs}Hz: {2*round(fs/4000)}, Got: {len(ag)-1}")
else:
    print("\n⚠ No voiced frames found in audio!")
