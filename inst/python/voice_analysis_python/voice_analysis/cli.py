"""
Command-Line Interface for Voice Analysis Toolbox
"""

import argparse
import sys
import numpy as np
import json
from pathlib import Path

from .core import analyze_voice_file


def main():
    parser = argparse.ArgumentParser(
        description='Voice Analysis Toolbox - Compute dysphonia measures from audio',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  voice-analysis audio.wav
  voice-analysis audio.wav --f0-min 75 --f0-max 300
  voice-analysis audio.wav --f0-algorithm PRAAT --output results.json
  
Citation:
  Tsanas, A., Little, M., McSharry, P., & Ramig, L. (2011).
  Nonlinear speech analysis algorithms mapped to a standard metric achieve
  clinically useful quantification of average Parkinson's disease symptom severity.
  Journal of the Royal Society Interface, 8(59), 842-855.
        """
    )
    
    parser.add_argument(
        'audio_file',
        type=str,
        help='Path to audio file (WAV format)'
    )
    
    parser.add_argument(
        '--f0-min',
        type=float,
        default=50,
        help='Minimum F0 in Hz (default: 50)'
    )
    
    parser.add_argument(
        '--f0-max',
        type=float,
        default=500,
        help='Maximum F0 in Hz (default: 500)'
    )
    
    parser.add_argument(
        '--f0-algorithm',
        type=str,
        choices=['SWIPE', 'PRAAT'],
        default='SWIPE',
        help='F0 estimation algorithm (default: SWIPE)'
    )
    
    parser.add_argument(
        '--output', '-o',
        type=str,
        default=None,
        help='Output file for results (JSON format)'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Verbose output'
    )
    
    args = parser.parse_args()
    
    # Check if file exists
    if not Path(args.audio_file).exists():
        print(f"Error: File not found: {args.audio_file}", file=sys.stderr)
        sys.exit(1)
    
    # Run analysis
    print(f"Analyzing: {args.audio_file}")
    print(f"F0 range: {args.f0_min}-{args.f0_max} Hz")
    print(f"F0 algorithm: {args.f0_algorithm}")
    print("-" * 60)
    
    try:
        measures, F0 = analyze_voice_file(
            args.audio_file,
            f0_min=args.f0_min,
            f0_max=args.f0_max,
            f0_algorithm=args.f0_algorithm
        )
        
        print("-" * 60)
        print(f"\n✓ Computed {len(measures)} measures")
        
        # F0 statistics
        F0_valid = F0[F0 > 0]
        if len(F0_valid) > 0:
            print(f"\nF0 Statistics:")
            print(f"  Mean: {np.mean(F0_valid):.2f} Hz")
            print(f"  Std:  {np.std(F0_valid):.2f} Hz")
            print(f"  Min:  {np.min(F0_valid):.2f} Hz")
            print(f"  Max:  {np.max(F0_valid):.2f} Hz")
        
        # Display measures
        if args.verbose:
            print(f"\nAll Measures:")
            print("-" * 60)
            for key in sorted(measures.keys()):
                print(f"  {key:40s}: {measures[key]:12.6f}")
        
        # Save to file if requested
        if args.output:
            output_data = {
                'file': args.audio_file,
                'parameters': {
                    'f0_min': args.f0_min,
                    'f0_max': args.f0_max,
                    'f0_algorithm': args.f0_algorithm
                },
                'measures': {k: float(v) if not np.isnan(v) else None 
                            for k, v in measures.items()},
                'F0_statistics': {
                    'mean': float(np.mean(F0_valid)) if len(F0_valid) > 0 else None,
                    'std': float(np.std(F0_valid)) if len(F0_valid) > 0 else None,
                    'min': float(np.min(F0_valid)) if len(F0_valid) > 0 else None,
                    'max': float(np.max(F0_valid)) if len(F0_valid) > 0 else None,
                }
            }
            
            with open(args.output, 'w') as f:
                json.dump(output_data, f, indent=2)
            
            print(f"\n✓ Results saved to: {args.output}")
        
    except Exception as e:
        print(f"\nError during analysis: {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
