"""
Main CreakDetector class - optimized for memory efficiency.
Extended version for R integration with all tracks.
"""

import numpy as np
from typing import Dict, Optional, Tuple
from pathlib import Path

class CreakDetectorExtended:
    """
    Creak detection using Union Method with extended output.
    Optimized for in-memory processing with R.
    Returns all intermediate features and probabilities.
    """
    
    def __init__(self, use_trained_ann: bool = True, use_reaper: bool = True, 
                 cd_threshold: float = 0.3, models_dir: Optional[str] = None):
        """
        Initialize detector.
        
        Parameters
        ----------
        use_trained_ann : bool
            Use trained MATLAB ANN weights
        use_reaper : bool
            Use REAPER for F0 (requires pyreaper)
        cd_threshold : float
            CD method threshold
        models_dir : str, optional
            Directory with .mat files (auto-detected if None)
        """
        self.use_trained_ann = use_trained_ann
        self.use_reaper = use_reaper
        self.cd_threshold = cd_threshold
        
        # Auto-detect models directory
        if models_dir is None:
            # Look for models relative to this file
            pkg_dir = Path(__file__).parent
            
            # Try local models directory
            models_dir = pkg_dir / 'models'
            if not models_dir.exists():
                # Try original location (union-creak-detection-method/CD_method)
                orig_location = pkg_dir.parent.parent.parent / 'union-creak-detection-method' / 'CD_method'
                if orig_location.exists():
                    models_dir = orig_location
        
        self.models_dir = Path(models_dir) if models_dir else None
        self.classifier = None
        self._f0_cache = None
        self._features_cache = None
        
        # Lazy loading
        if use_trained_ann and self.models_dir and self.models_dir.exists():
            self._load_classifier()
    
    def _load_classifier(self):
        """Lazy load classifier."""
        if self.classifier is None:
            from .classification import TrainedANN
            try:
                self.classifier = TrainedANN(self.models_dir)
            except Exception as e:
                print(f"Warning: Could not load trained ANN: {e}")
                print("Using simple classifier")
                from .classification import SimpleClassifier
                self.classifier = SimpleClassifier()
    
    def process_extended(self, audio: np.ndarray, sample_rate: int = 16000,
                        use_am: bool = True, use_cd: bool = True,
                        frame_shift_ms: float = 10.0,
                        return_features: bool = False,
                        return_probabilities: bool = False) -> Dict[str, np.ndarray]:
        """
        Process audio with extended output including all intermediate tracks.
        
        Parameters
        ----------
        audio : np.ndarray
            Audio signal (1D array)
        sample_rate : int
            Sampling rate
        use_am : bool
            Use AM method
        use_cd : bool
            Use CD method
        frame_shift_ms : float
            Frame shift in ms
        return_features : bool
            Include feature matrix in output
        return_probabilities : bool
            Include classification probabilities
            
        Returns
        -------
        dict
            Results with all tracks as NumPy arrays
        """
        audio = np.asarray(audio, dtype=np.float64).flatten()
        duration = len(audio) / sample_rate
        
        # Time vector
        frame_shift_sec = frame_shift_ms / 1000.0
        n_frames = int(duration / frame_shift_sec) + 1
        time = np.arange(n_frames) * frame_shift_sec
        
        results = {
            'time': time,
            'sample_rate': sample_rate,
            'duration': duration
        }
        
        # AM Method with F0 extraction
        if use_am:
            am_decisions, antimode, f0_contour = self._am_method_extended(
                audio, sample_rate, time
            )
            results['am_decisions'] = am_decisions
            results['antimode'] = antimode
            results['F0'] = f0_contour
            self._f0_cache = (f0_contour, time)
        
        # CD Method with probabilities and features
        if use_cd:
            cd_decisions, cd_prob, features = self._cd_method_extended(
                audio, sample_rate, time, return_features, return_probabilities
            )
            results['cd_decisions'] = cd_decisions
            
            if return_probabilities:
                results['CD_prob'] = cd_prob
            
            if return_features:
                results['features'] = features
            
            self._features_cache = features
        
        # Union
        if use_am and use_cd:
            results['union_decisions'] = np.logical_or(
                results['am_decisions'],
                results['cd_decisions']
            ).astype(np.int32)
        elif use_cd:
            results['union_decisions'] = results['cd_decisions']
        elif use_am:
            results['union_decisions'] = results['am_decisions']
        
        return results
    
    def _am_method_extended(self, audio: np.ndarray, sr: int, 
                           time: np.ndarray) -> Tuple[np.ndarray, float, np.ndarray]:
        """
        AM method with F0 contour output.
        
        Returns
        -------
        decisions : ndarray
            Binary creak decisions
        antimode : float
            Antimode value in Hz
        f0_contour : ndarray
            F0 contour interpolated to time vector
        """
        from .am_method import extract_f0, detect_antimode
        
        # Get F0
        if self.use_reaper:
            try:
                f0, voicing, f0_times = extract_f0(audio, sr, method='reaper')
            except:
                f0, voicing, f0_times = extract_f0(audio, sr, method='simple')
        else:
            f0, voicing, f0_times = extract_f0(audio, sr, method='simple')
        
        # Interpolate to frame times
        f0_inter = np.interp(time, f0_times, f0)
        vuv_inter = np.interp(time, f0_times, voicing)
        vuv_inter = (vuv_inter > 0.5).astype(np.int32)
        
        # Detect antimode
        voiced_f0 = f0_inter[vuv_inter == 1]
        antimode = detect_antimode(voiced_f0) if len(voiced_f0) > 10 else 0.0
        
        # Mark frames below antimode as creaky
        decisions = np.zeros(len(time), dtype=np.int32)
        if antimode > 0:
            decisions[(f0_inter < antimode) & (f0_inter > 0)] = 1
        
        return decisions, antimode, f0_inter
    
    def _cd_method_extended(self, audio: np.ndarray, sr: int, time: np.ndarray,
                           return_features: bool, 
                           return_probabilities: bool) -> Tuple[np.ndarray, Optional[np.ndarray], Optional[np.ndarray]]:
        """
        CD method with probabilities and features.
        
        Returns
        -------
        decisions : ndarray
            Binary creak decisions
        probabilities : ndarray or None
            Classification probabilities (if return_probabilities=True)
        features : ndarray or None
            Feature matrix (n_frames, 36) (if return_features=True)
        """
        from .features import get_all_features
        
        # Extract features
        features, feat_time = get_all_features(
            audio,
            sample_rate=sr,
            frame_shift_ms=(time[1] - time[0]) * 1000 if len(time) > 1 else 10.0
        )
        
        # Classify
        if self.classifier is None:
            self._load_classifier()
        
        decisions, probabilities = self.classifier.classify(
            features, 
            threshold=self.cd_threshold
        )
        
        # Interpolate to match time vector if needed
        if len(decisions) != len(time):
            # Interpolate decisions
            dec_inter = np.interp(time, feat_time, decisions.astype(np.float64))
            decisions = (dec_inter > 0.5).astype(np.int32)
            
            # Interpolate probabilities
            if return_probabilities:
                probabilities = np.interp(time, feat_time, probabilities)
            
            # Interpolate features (each column separately)
            if return_features:
                features_inter = np.zeros((len(time), features.shape[1]))
                for col in range(features.shape[1]):
                    features_inter[:, col] = np.interp(time, feat_time, features[:, col])
                features = features_inter
        
        # Return appropriate values
        prob_out = probabilities if return_probabilities else None
        feat_out = features if return_features else None
        
        return decisions, prob_out, feat_out
    
    def process(self, audio: np.ndarray, sample_rate: int = 16000,
                use_am: bool = True, use_cd: bool = True,
                frame_shift_ms: float = 10.0) -> Dict[str, np.ndarray]:
        """
        Simple process method (backward compatibility).
        Returns only decisions and antimode.
        """
        return self.process_extended(
            audio, sample_rate, use_am, use_cd, frame_shift_ms,
            return_features=False, return_probabilities=False
        )


# Backward compatibility alias
CreakDetector = CreakDetectorExtended
