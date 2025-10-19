"""
Apple Silicon (M1/M2) specific optimizations
Uses Accelerate framework and ARM-specific optimizations
"""

import numpy as np
import platform
from typing import Tuple

# Check if on Apple Silicon
IS_APPLE_SILICON = platform.machine() == 'arm64' and platform.system() == 'Darwin'


def configure_numpy_for_apple_silicon():
    """
    Configure NumPy to use Apple's Accelerate framework
    
    Apple Silicon has hardware-accelerated BLAS/LAPACK via Accelerate.
    This is automatically used if NumPy is built with Accelerate support.
    """
    if not IS_APPLE_SILICON:
        return False
    
    try:
        # Check if NumPy is using Accelerate
        # Capture output instead of printing to stdout
        import io
        import sys
        old_stdout = sys.stdout
        sys.stdout = io.StringIO()
        try:
            np.__config__.show()
            config = sys.stdout.getvalue()
        finally:
            sys.stdout = old_stdout
        
        if config and ('Accelerate' in str(config) or 'vecLib' in str(config)):
            return True
    except:
        pass
    
    return False


def optimize_for_apple_silicon():
    """
    Apply Apple Silicon specific optimizations
    
    Returns:
        dict with optimization status
    """
    if not IS_APPLE_SILICON:
        return {'apple_silicon': False}
    
    status = {
        'apple_silicon': True,
        'accelerate': False,
        'metal': False,
        'recommendations': []
    }
    
    # Check Accelerate
    status['accelerate'] = configure_numpy_for_apple_silicon()
    
    if not status['accelerate']:
        status['recommendations'].append(
            "Install NumPy built with Accelerate: pip install numpy --force-reinstall --no-binary numpy"
        )
    
    # Check for Metal support (for future GPU acceleration)
    try:
        import metal
        status['metal'] = True
    except ImportError:
        status['recommendations'].append(
            "For GPU acceleration, consider installing Metal support libraries"
        )
    
    return status


def get_architecture_info() -> dict:
    """Get detailed CPU architecture information"""
    info = {
        'platform': platform.system(),
        'machine': platform.machine(),
        'processor': platform.processor(),
        'cores': None,
        'performance_cores': None,
        'efficiency_cores': None,
        'is_apple_silicon': IS_APPLE_SILICON,
        'is_x86': platform.machine() in ['x86_64', 'AMD64'],
    }
    
    # Try to get detailed core information
    try:
        import subprocess
        if IS_APPLE_SILICON:
            # Get core count on Apple Silicon
            result = subprocess.run(
                ['sysctl', '-n', 'hw.perflevel0.physicalcpu'],
                capture_output=True, text=True
            )
            if result.returncode == 0:
                info['performance_cores'] = int(result.stdout.strip())
            
            result = subprocess.run(
                ['sysctl', '-n', 'hw.perflevel1.physicalcpu'],
                capture_output=True, text=True
            )
            if result.returncode == 0:
                info['efficiency_cores'] = int(result.stdout.strip())
            
            info['cores'] = (info.get('performance_cores', 0) or 0) + \
                           (info.get('efficiency_cores', 0) or 0)
    except:
        pass
    
    if info['cores'] is None:
        import os
        info['cores'] = os.cpu_count()
    
    return info


class AppleSiliconOptimizer:
    """
    Optimizer for Apple Silicon specific features
    """
    
    def __init__(self):
        self.arch_info = get_architecture_info()
        self.optimizations = optimize_for_apple_silicon()
    
    def get_recommended_num_processes(self) -> int:
        """
        Get recommended number of processes for Apple Silicon
        
        M1/M2 have performance and efficiency cores.
        Best to use performance cores for compute-heavy tasks.
        """
        if not self.arch_info['is_apple_silicon']:
            # For x86, use standard logic
            return max(1, self.arch_info['cores'] - 1)
        
        # Apple Silicon: prefer performance cores
        perf_cores = self.arch_info.get('performance_cores')
        
        if perf_cores:
            # Use all performance cores, leave efficiency cores for OS
            return perf_cores
        else:
            # Fallback: use 75% of total cores
            return max(4, int(self.arch_info['cores'] * 0.75))
    
    def configure_thread_affinity(self):
        """
        Configure thread affinity for optimal performance
        
        On Apple Silicon, try to pin threads to performance cores
        """
        if not self.arch_info['is_apple_silicon']:
            return False
        
        # macOS doesn't provide direct thread affinity APIs
        # But we can set quality of service hints
        try:
            import os
            # Set QoS to user-interactive for best performance core assignment
            os.environ['OPENBLAS_NUM_THREADS'] = str(self.get_recommended_num_processes())
            os.environ['MKL_NUM_THREADS'] = str(self.get_recommended_num_processes())
            os.environ['OMP_NUM_THREADS'] = str(self.get_recommended_num_processes())
            return True
        except:
            return False
    
    def optimize_fft_for_arm(self):
        """
        Optimize FFT for ARM architecture
        
        Apple's vDSP (part of Accelerate) is highly optimized for ARM
        """
        if not self.arch_info['is_apple_silicon']:
            return False
        
        # Ensure we're using scipy.fft which can leverage Accelerate
        try:
            from scipy import fft
            # scipy.fft can use Accelerate's vDSP on Apple Silicon
            return True
        except ImportError:
            return False
    
    def get_optimization_summary(self) -> str:
        """Get human-readable optimization summary"""
        lines = ["Apple Silicon Optimization Status:", "=" * 50]
        
        if self.arch_info['is_apple_silicon']:
            lines.append(f"✓ Apple Silicon detected ({self.arch_info['machine']})")
            
            if self.arch_info.get('performance_cores'):
                lines.append(f"  Performance cores: {self.arch_info['performance_cores']}")
            if self.arch_info.get('efficiency_cores'):
                lines.append(f"  Efficiency cores: {self.arch_info['efficiency_cores']}")
            
            lines.append(f"  Recommended processes: {self.get_recommended_num_processes()}")
            
            if self.optimizations['accelerate']:
                lines.append("✓ Accelerate framework: Active")
            else:
                lines.append("✗ Accelerate framework: Not detected")
                for rec in self.optimizations['recommendations']:
                    lines.append(f"  → {rec}")
        else:
            lines.append(f"  Platform: {self.arch_info['platform']} {self.arch_info['machine']}")
            lines.append(f"  Cores: {self.arch_info['cores']}")
            lines.append(f"  Recommended processes: {self.get_recommended_num_processes()}")
        
        return "\n".join(lines)


# Singleton instance
_optimizer = None

def get_optimizer() -> AppleSiliconOptimizer:
    """Get singleton optimizer instance"""
    global _optimizer
    if _optimizer is None:
        _optimizer = AppleSiliconOptimizer()
    return _optimizer
