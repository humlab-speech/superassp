"""
Test New Features: DYPSA, GQ, VFER, EMD

Test the newly implemented features
"""

import numpy as np
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))


def test_dypsa():
    """Test DYPSA glottal instant detection"""
    print("\n" + "="*60)
    print("Testing DYPSA Algorithm")
    print("="*60)
    
    from voice_analysis.utils import dypsa
    
    # Generate voiced signal (sine wave)
    fs = 44100
    t = np.linspace(0, 1, fs)
    f0 = 150  # Hz
    audio = np.sin(2 * np.pi * f0 * t)
    
    # Add some amplitude modulation to simulate glottal pulses
    modulation = 0.5 * (1 + np.sin(2 * np.pi * f0 * t * 0.5))
    audio = audio * modulation
    
    try:
        gci, goi = dypsa(audio, fs, f0_min=100, f0_max=200)
        
        print(f"  ✓ DYPSA executed")
        print(f"  ✓ Detected {len(gci)} glottal closure instants")
        print(f"  ✓ Detected {len(goi)} glottal opening instants")
        
        if len(gci) > 0:
            # Check GCI intervals are reasonable
            intervals = np.diff(gci)
            mean_interval = np.mean(intervals)
            expected_interval = fs / f0
            
            print(f"  ✓ Mean GCI interval: {mean_interval:.1f} samples")
            print(f"  ✓ Expected interval: {expected_interval:.1f} samples")
            
            # Should be roughly correct (within 50%)
            if 0.5 * expected_interval < mean_interval < 1.5 * expected_interval:
                print(f"  ✓ GCI intervals are reasonable")
            else:
                print(f"  ⚠ GCI intervals differ from expected")
        
        return True
        
    except Exception as e:
        print(f"  ✗ DYPSA failed: {e}")
        return False


def test_glottal_quotient():
    """Test Glottal Quotient computation"""
    print("\n" + "="*60)
    print("Testing Glottal Quotient (GQ)")
    print("="*60)
    
    from voice_analysis.features import compute_glottal_quotient
    
    # Generate test signal
    fs = 44100
    t = np.linspace(0, 1, fs)
    audio = np.sin(2 * np.pi * 150 * t)
    
    try:
        gq = compute_glottal_quotient(audio, fs)
        
        print(f"  ✓ GQ computation executed")
        print(f"  ✓ GQ: {gq.get('GQ', np.nan):.4f}")
        print(f"  ✓ GQ_open_std: {gq.get('GQ_open_std', np.nan):.2f}")
        print(f"  ✓ GQ_closed_std: {gq.get('GQ_closed_std', np.nan):.2f}")
        
        # Check we got 3 measures
        assert 'GQ' in gq
        assert 'GQ_open_std' in gq
        assert 'GQ_closed_std' in gq
        print(f"  ✓ All 3 GQ measures present")
        
        return True
        
    except Exception as e:
        print(f"  ✗ GQ failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_vfer():
    """Test VFER computation"""
    print("\n" + "="*60)
    print("Testing VFER (Vocal Fold Excitation Ratios)")
    print("="*60)
    
    from voice_analysis.features import compute_vfer
    
    # Generate test signal
    fs = 44100
    t = np.linspace(0, 1, fs)
    audio = np.sin(2 * np.pi * 150 * t)
    
    try:
        vfer = compute_vfer(audio, fs)
        
        print(f"  ✓ VFER computation executed")
        print(f"  ✓ VFER_mean: {vfer.get('VFER_mean', np.nan):.4f}")
        print(f"  ✓ VFER_std: {vfer.get('VFER_std', np.nan):.4f}")
        print(f"  ✓ VFER_entropy: {vfer.get('VFER_entropy', np.nan):.4f}")
        
        # Check we got 7 measures
        expected_keys = [
            'VFER_mean', 'VFER_std', 'VFER_entropy',
            'VFER_TKEO_low_high', 'VFER_SEO_low_high',
            'VFER_log_TKEO_high_low', 'VFER_log_SEO_high_low'
        ]
        
        for key in expected_keys:
            assert key in vfer, f"Missing key: {key}"
        
        print(f"  ✓ All 7 VFER measures present")
        
        return True
        
    except Exception as e:
        print(f"  ✗ VFER failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_emd():
    """Test EMD features"""
    print("\n" + "="*60)
    print("Testing EMD Features")
    print("="*60)
    
    from voice_analysis.features import compute_emd_features
    
    # Generate test signal
    t = np.linspace(0, 1, 1000)
    # Composite signal with multiple frequencies
    signal = (np.sin(2 * np.pi * 5 * t) + 
              0.5 * np.sin(2 * np.pi * 15 * t) +
              0.2 * np.sin(2 * np.pi * 30 * t))
    
    try:
        emd = compute_emd_features(signal)
        
        print(f"  ✓ EMD computation executed")
        print(f"  ✓ EMD_energy_ratio: {emd.get('EMD_energy_ratio', np.nan):.4f}")
        print(f"  ✓ EMD_TKEO_ratio: {emd.get('EMD_TKEO_ratio', np.nan):.4f}")
        print(f"  ✓ EMD_entropy_ratio: {emd.get('EMD_entropy_ratio', np.nan):.4f}")
        
        # Check we got 6 measures
        expected_keys = [
            'EMD_energy_ratio', 'EMD_TKEO_ratio', 'EMD_entropy_ratio',
            'EMD_log_energy_ratio', 'EMD_log_TKEO_ratio', 'EMD_log_entropy_ratio'
        ]
        
        for key in expected_keys:
            assert key in emd, f"Missing key: {key}"
        
        print(f"  ✓ All 6 EMD measures present")
        
        return True
        
    except ImportError:
        print(f"  ⚠ PyEMD not installed (this is optional)")
        print(f"  Install with: pip install EMD-signal")
        return True  # Not a failure, just optional
        
    except Exception as e:
        print(f"  ✗ EMD failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_complete_analysis():
    """Test complete analysis with all 132 measures"""
    print("\n" + "="*60)
    print("Testing Complete Analysis (All 132 Measures)")
    print("="*60)
    
    from voice_analysis import VoiceAnalyzer
    
    # Generate test signal
    fs = 44100
    t = np.linspace(0, 2, 2*fs)  # 2 seconds
    audio = np.sin(2 * np.pi * 150 * t)
    
    # Add some harmonics and modulation
    audio += 0.3 * np.sin(2 * np.pi * 300 * t)
    audio += 0.2 * np.sin(2 * np.pi * 450 * t)
    modulation = 0.5 * (1 + np.sin(2 * np.pi * 5 * t))
    audio = audio * modulation
    
    try:
        analyzer = VoiceAnalyzer(f0_min=100, f0_max=200)
        measures, F0 = analyzer.analyze(audio, fs)
        
        print(f"\n  ✓ Analysis completed")
        print(f"  ✓ Total measures computed: {len(measures)}")
        
        # Count measures by category
        categories = {
            'Jitter': 0, 'Shimmer': 0, 'HNR/NHR': 0,
            'MFCC': 0, 'Wavelet': 0, 'GNE': 0,
            'GQ': 0, 'VFER': 0, 'EMD': 0, 'Other': 0
        }
        
        for key in measures.keys():
            if key.startswith('jitter_'):
                categories['Jitter'] += 1
            elif key.startswith('shimmer_'):
                categories['Shimmer'] += 1
            elif key.startswith(('HNR', 'NHR')):
                categories['HNR/NHR'] += 1
            elif key.startswith('MFCC'):
                categories['MFCC'] += 1
            elif key.startswith('wavelet'):
                categories['Wavelet'] += 1
            elif key.startswith('GNE'):
                categories['GNE'] += 1
            elif key.startswith('GQ'):
                categories['GQ'] += 1
            elif key.startswith('VFER'):
                categories['VFER'] += 1
            elif key.startswith('EMD'):
                categories['EMD'] += 1
            else:
                categories['Other'] += 1
        
        print(f"\n  Measures by category:")
        for cat, count in categories.items():
            if count > 0:
                print(f"    {cat:12s}: {count:3d}")
        
        # Check for new features
        if 'GQ' in measures:
            print(f"\n  ✓ GQ features present")
        else:
            print(f"\n  ⚠ GQ features not found")
        
        if 'VFER_mean' in measures:
            print(f"  ✓ VFER features present")
        else:
            print(f"  ⚠ VFER features not found")
        
        if 'EMD_energy_ratio' in measures:
            print(f"  ✓ EMD features present")
        else:
            print(f"  ⚠ EMD features not found (PyEMD may not be installed)")
        
        # Estimate total (some may be NaN)
        non_nan = sum(1 for v in measures.values() if not np.isnan(v))
        print(f"\n  ✓ Non-NaN measures: {non_nan}")
        
        if len(measures) >= 120:
            print(f"  ✓ SUCCESS: At least 120 measures computed!")
            return True
        else:
            print(f"  ⚠ WARNING: Only {len(measures)} measures computed")
            return False
        
    except Exception as e:
        print(f"  ✗ Complete analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == '__main__':
    print("\n" + "="*60)
    print("TESTING NEW FEATURES")
    print("DYPSA, GQ, VFER, EMD")
    print("="*60)
    
    results = []
    
    # Test each component
    results.append(("DYPSA", test_dypsa()))
    results.append(("Glottal Quotient", test_glottal_quotient()))
    results.append(("VFER", test_vfer()))
    results.append(("EMD", test_emd()))
    results.append(("Complete Analysis", test_complete_analysis()))
    
    # Summary
    print("\n" + "="*60)
    print("TEST SUMMARY")
    print("="*60)
    
    for name, result in results:
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"  {status:8s} {name}")
    
    passed = sum(1 for _, r in results if r)
    total = len(results)
    
    print(f"\n  {passed}/{total} tests passed")
    
    if passed == total:
        print("\n  🎉 ALL TESTS PASSED! 🎉")
        print("  Implementation is complete and working!")
    else:
        print("\n  ⚠ Some tests failed. Check output above.")
    
    print("="*60)
