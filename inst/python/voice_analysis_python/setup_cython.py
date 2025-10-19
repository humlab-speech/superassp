"""
Setup script for Voice Analysis Toolbox with Cython optimizations

Supports:
- Pure Python installation (fallback)
- Cython compilation for performance
- Platform-specific optimizations (ARM NEON, AVX-512)
- R reticulate compatibility
"""

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import sys
import os
import platform
import numpy

# Try to import Cython
try:
    from Cython.Build import cythonize
    CYTHON_AVAILABLE = True
except ImportError:
    CYTHON_AVAILABLE = False
    print("Warning: Cython not available. Installing without compiled extensions.")
    print("For better performance, install Cython: pip install cython")


class BuildExtWithFallback(build_ext):
    """
    Custom build_ext that gracefully handles compilation failures
    """
    def run(self):
        try:
            build_ext.run(self)
        except Exception as e:
            print(f"Warning: Cython compilation failed: {e}")
            print("Falling back to pure Python/Numba implementation")
            print("Performance will be slightly reduced")
    
    def build_extension(self, ext):
        try:
            build_ext.build_extension(self, ext)
        except Exception as e:
            print(f"Warning: Failed to compile {ext.name}: {e}")
            print(f"Skipping this extension (will use fallback)")


def get_compile_args():
    """
    Get platform-specific compilation arguments
    """
    compile_args = ["-O3"]
    
    system = platform.system()
    machine = platform.machine()
    
    # Platform-specific optimizations
    if system == "Darwin":  # macOS
        if machine == "arm64":  # M1/M2/M3
            print("Detected Apple Silicon - enabling ARM NEON optimizations")
            compile_args.extend([
                "-march=armv8-a+simd",
                "-mtune=apple-m1",
                "-ftree-vectorize",
            ])
        else:  # Intel Mac
            compile_args.extend([
                "-march=native",
                "-mtune=native",
            ])
    
    elif system == "Linux":
        # Check for AVX-512 support (common on AMD EPYC and modern Intel Xeon)
        try:
            import subprocess
            cpu_info = subprocess.check_output(['lscpu'], text=True)
            
            if 'avx512' in cpu_info.lower():
                print("Detected AVX-512 support - enabling AVX-512 optimizations")
                compile_args.extend([
                    "-march=native",
                    "-mavx512f",
                    "-mavx512dq",
                    "-mfma",
                ])
            elif 'avx2' in cpu_info.lower():
                print("Detected AVX2 support - enabling AVX2 optimizations")
                compile_args.extend([
                    "-march=native",
                    "-mavx2",
                ])
            else:
                compile_args.append("-march=native")
        except:
            # Fallback: use native without specific flags
            compile_args.append("-march=native")
    
    elif system == "Windows":
        # Windows MSVC flags
        compile_args = ["/O2", "/favor:INTEL64"]
    
    return compile_args


def get_extensions():
    """
    Define Cython extensions
    """
    if not CYTHON_AVAILABLE:
        return []
    
    compile_args = get_compile_args()
    
    extensions = [
        Extension(
            "voice_analysis.features.rpde_cython",
            ["voice_analysis/features/rpde_cython.pyx"],
            include_dirs=[numpy.get_include()],
            extra_compile_args=compile_args,
            extra_link_args=[],
        ),
        Extension(
            "voice_analysis.utils.perturbation_cython",
            ["voice_analysis/utils/perturbation_cython.pyx"],
            include_dirs=[numpy.get_include()],
            extra_compile_args=compile_args,
            extra_link_args=[],
        ),
    ]
    
    return extensions


# Read version from __init__.py
def get_version():
    """Extract version from package __init__.py"""
    init_file = os.path.join('voice_analysis', '__init__.py')
    if os.path.exists(init_file):
        with open(init_file, 'r') as f:
            for line in f:
                if line.startswith('__version__'):
                    return line.split('=')[1].strip().strip('"').strip("'")
    return "0.1.0"


# Read long description from README
long_description = ""
if os.path.exists('README.md'):
    with open('README.md', 'r', encoding='utf-8') as f:
        long_description = f.read()


# Core dependencies (required)
install_requires = [
    'numpy>=1.20.0',
    'scipy>=1.7.0',
    'soundfile>=0.10.0',
    'pysptk>=0.2.0',
    'librosa>=0.9.0',
    'nolds>=0.5.0',
]

# Performance dependencies (recommended)
performance_extras = [
    'numba>=0.56.0',
    'joblib>=1.0.0',
]

# Optional feature dependencies
optional_extras = [
    'PyWavelets>=1.3.0',
    'EMD-signal>=1.6.0',  # Note: may need folder rename pyemd → PyEMD
]

# R integration dependencies
r_extras = [
    'pandas>=1.3.0',  # For R data.frame conversion
]

# Development dependencies
dev_extras = [
    'pytest>=6.0.0',
    'pytest-cov>=2.0.0',
    'cython>=0.29.0',
]

# All extras
all_extras = performance_extras + optional_extras + r_extras + dev_extras


# Setup configuration
setup_config = {
    'name': 'voice-analysis',
    'version': get_version(),
    'description': 'Comprehensive voice analysis toolkit - MATLAB reimplementation in Python',
    'long_description': long_description,
    'long_description_content_type': 'text/markdown',
    'author': 'Voice Analysis Team',
    'author_email': '',
    'url': 'https://github.com/your-repo/voice-analysis',
    'packages': [
        'voice_analysis',
        'voice_analysis.features',
        'voice_analysis.f0_estimation',
        'voice_analysis.utils',
        'voice_analysis.dypsa',
    ],
    'install_requires': install_requires,
    'extras_require': {
        'performance': performance_extras,
        'optional': optional_extras,
        'r': r_extras,
        'dev': dev_extras,
        'all': all_extras,
    },
    'python_requires': '>=3.8',
    'classifiers': [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Medical Science Apps.',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Cython',
    ],
    'keywords': 'voice analysis speech jitter shimmer RPDE DFA acoustic features',
    'entry_points': {
        'console_scripts': [
            'voice-analysis=voice_analysis.cli:main',
        ],
    },
}

# Add Cython extensions if available
if CYTHON_AVAILABLE:
    extensions = get_extensions()
    
    setup_config['ext_modules'] = cythonize(
        extensions,
        compiler_directives={
            'language_level': 3,
            'boundscheck': False,
            'wraparound': False,
            'cdivision': True,
            'initializedcheck': False,
            'nonecheck': False,
        },
        annotate=False,  # Set to True to generate HTML annotation files
    )
    setup_config['cmdclass'] = {'build_ext': BuildExtWithFallback}
    
    print("\n" + "="*70)
    print("Building with Cython optimizations")
    print("Target platform:", platform.system(), platform.machine())
    print("Compiler flags:", ' '.join(get_compile_args()))
    print("="*70 + "\n")
else:
    print("\n" + "="*70)
    print("Building without Cython (pure Python/Numba mode)")
    print("For better performance: pip install cython && pip install --no-binary :all: voice-analysis")
    print("="*70 + "\n")

setup(**setup_config)
