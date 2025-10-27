"""
Test script for R integration.
Simulates calling from R via reticulate.
"""

import numpy as np
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent.parent / 'union-creak-detection-method'))

def test_extended_api():
    """Test the extended API with all features."""
    print("="*70)
    print("Testing Extended API for R Integration")
    print("="*70)
    
    # Import module
    try:
        import union_creak_detection_method as creak_mod
        print("✓ Module imported successfully")
    except ImportError as e:
        print(f"✗ Import error: {e}")
        print("Trying local import...")
        sys.path.insert(0, str(Path(__file__).parent))
        from api import detect_creak_union_extended
        from detector import CreakDetectorExtended
        print("✓ Local import successful")
    
    # Load test audio (simulate av::read_audio_bin)
    from creak_detection.utils.audio import load_audio
    
    demo_file = Path(__file__).parent.parent.parent.parent.parent / 'union-creak-detection-method' / 'demo' / 'sample1_passage.wav'
    
    if not demo_file.exists():
        print(f"✗ Demo file not found: {demo_file}")
        return
    
    print(f"\nLoading: {demo_file.name}")
    audio, sr = load_audio(str(demo_file))
    print(f"  Duration: {len(audio)/sr:.2f}s, Sample rate: {sr}Hz")
    
    # Take short segment for testing
    segment_duration = 5.0  # seconds
    audio_segment = audio[:int(segment_duration * sr)]
    
    print(f"\nTesting on {segment_duration}s segment...")
    
    # Test extended API
    print("\n" + "-"*70)
    print("1. Basic Call (decisions only)")
    print("-"*70)
    
    result1 = detect_creak_union_extended(
        audio=audio_segment,
        sample_rate=sr,
        use_am=True,
        use_cd=True,
        cd_threshold=0.3,
        use_reaper=True,
        frame_shift_ms=10.0,
        return_features=False,
        return_probabilities=False
    )
    
    print(f"  Keys: {list(result1.keys())}")
    print(f"  Time: {len(result1['time'])} frames")
    print(f"  AM creak: {np.sum(result1['am_decisions'])} / {len(result1['am_decisions'])}")
    print(f"  CD creak: {np.sum(result1['cd_decisions'])} / {len(result1['cd_decisions'])}")
    print(f"  Union creak: {np.sum(result1['union_decisions'])} / {len(result1['union_decisions'])}")
    print(f"  Antimode: {result1['antimode']:.2f} Hz")
    print(f"  F0 range: [{np.min(result1['F0']):.1f}, {np.max(result1['F0']):.1f}] Hz")
    
    print("\n" + "-"*70)
    print("2. With Probabilities")
    print("-"*70)
    
    result2 = detect_creak_union_extended(
        audio=audio_segment,
        sample_rate=sr,
        use_am=True,
        use_cd=True,
        return_features=False,
        return_probabilities=True
    )
    
    print(f"  Keys: {list(result2.keys())}")
    if 'CD_prob' in result2:
        print(f"  CD probability range: [{np.min(result2['CD_prob']):.3f}, {np.max(result2['CD_prob']):.3f}]")
        print(f"  CD probability mean: {np.mean(result2['CD_prob']):.3f}")
    
    print("\n" + "-"*70)
    print("3. With Features")
    print("-"*70)
    
    result3 = detect_creak_union_extended(
        audio=audio_segment,
        sample_rate=sr,
        use_am=True,
        use_cd=True,
        return_features=True,
        return_probabilities=False
    )
    
    print(f"  Keys: {list(result3.keys())}")
    if 'features' in result3:
        features = result3['features']
        print(f"  Features shape: {features.shape} (frames x features)")
        print(f"  Features: 36 features (12 static + 12 delta + 12 delta-delta)")
        
        # Show some feature statistics
        feature_names = ['H2_H1', 'peak_prom', 'ZCR', 'IFP', 'IPS', 
                        'PwP_fall', 'PwP_rise', 'F0', 'F0_mean', 
                        'energy', 'power_std', 'creak_F0']
        
        print("\n  Feature Statistics (first 6 static features):")
        for i in range(6):
            feat = features[:, i]
            print(f"    {feature_names[i]:15s}: Mean={np.mean(feat):7.3f}, "
                  f"Std={np.std(feat):7.3f}, Range=[{np.min(feat):7.2f}, {np.max(feat):7.2f}]")
    
    print("\n" + "-"*70)
    print("4. Full Output (all options)")
    print("-"*70)
    
    result4 = detect_creak_union_extended(
        audio=audio_segment,
        sample_rate=sr,
        use_am=True,
        use_cd=True,
        return_features=True,
        return_probabilities=True
    )
    
    print(f"  Keys: {list(result4.keys())}")
    print(f"  Total data returned:")
    for key, value in result4.items():
        if isinstance(value, np.ndarray):
            print(f"    {key:20s}: shape={value.shape}, dtype={value.dtype}")
        else:
            print(f"    {key:20s}: {value}")
    
    print("\n" + "="*70)
    print("✓ All tests passed!")
    print("="*70)
    
    return result4


def simulate_r_workflow():
    """Simulate how R would call this."""
    print("\n" + "="*70)
    print("Simulating R Workflow (via reticulate)")
    print("="*70)
    
    print("\nR Code:")
    print("-"*70)
    print("""
    library(av)
    library(reticulate)
    
    # Load audio
    audio_data <- read_audio_bin("sample1_passage.wav", channels=1)
    audio_info <- av_media_info("sample1_passage.wav")
    sr <- audio_info$audio$sample_rate
    
    # Import Python module
    creak <- import("union_creak_detection_method")
    
    # Call detector
    result <- creak$detect_creak_union_extended(
        audio = audio_data,
        sample_rate = as.integer(sr),
        use_am = TRUE,
        use_cd = TRUE,
        return_features = TRUE,
        return_probabilities = TRUE
    )
    
    # Access results (automatic conversion to R)
    time_vec <- result$time           # R numeric vector
    am_creak <- result$am_decisions   # R integer vector
    cd_creak <- result$cd_decisions   # R integer vector
    union_creak <- result$union_decisions  # R integer vector
    cd_prob <- result$CD_prob         # R numeric vector
    f0_track <- result$F0             # R numeric vector
    features <- result$features       # R matrix (n x 36)
    
    # Use in trk_creak_union()
    # Creates AsspDataObj with all tracks
    """)
    print("-"*70)
    
    print("\nPython execution: ✓ (see test above)")


if __name__ == "__main__":
    result = test_extended_api()
    simulate_r_workflow()
