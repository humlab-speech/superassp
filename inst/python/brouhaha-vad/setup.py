from pathlib import Path
import sys
import numpy as np

from setuptools import find_packages, setup, Extension

here = Path(__file__).parent

with open(here / "requirements.txt") as file:
    requirements = file.read().splitlines()

with open(here / "README.md") as file:
    readme = file.read()

# ============================================================================
# Cython Extension Configuration
# ============================================================================

# Try to import Cython
USE_CYTHON = False
try:
    from Cython.Build import cythonize
    USE_CYTHON = True
    print("Cython found - building optimized extensions")
except ImportError:
    print("Cython not found - skipping optimized extensions")
    print("Install with: pip install cython")

# Define extensions if Cython is available
extensions = []
if USE_CYTHON:
    # Compiler flags for maximum performance
    extra_compile_args = []
    extra_link_args = []

    if sys.platform == 'darwin':  # macOS
        extra_compile_args = ['-O3', '-march=native', '-ffast-math']
    elif sys.platform.startswith('linux'):  # Linux
        extra_compile_args = ['-O3', '-march=native', '-ffast-math', '-fopenmp']
        extra_link_args = ['-fopenmp']
    elif sys.platform == 'win32':  # Windows
        extra_compile_args = ['/O2', '/openmp']

    # Define Cython extensions
    extensions = [
        Extension(
            "brouhaha.utils.collate_fast",
            ["brouhaha/utils/collate_fast.pyx"],
            include_dirs=[np.get_include()],
            extra_compile_args=extra_compile_args,
            extra_link_args=extra_link_args,
            define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],
        ),
        Extension(
            "brouhaha.utils.metrics_fast",
            ["brouhaha/utils/metrics_fast.pyx"],
            include_dirs=[np.get_include()],
            extra_compile_args=extra_compile_args,
            extra_link_args=extra_link_args,
            define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],
        ),
    ]

    # Cythonize with compiler directives
    extensions = cythonize(
        extensions,
        compiler_directives={
            'language_level': 3,
            'boundscheck': False,
            'wraparound': False,
            'cdivision': True,
            'initializedcheck': False,
            'nonecheck': False,
        },
        annotate=False,  # Set to True to generate HTML annotation files for optimization analysis
    )

# ============================================================================
# Setup Configuration
# ============================================================================

setup(
    name="brouhaha",
    version='0.9.1',  # Incremented for optimization release
    packages=find_packages(),
    description="Pyannote extension for SNR and C50 predictions - Optimized Edition",
    long_description=readme,
    long_description_content_type="text/markdown",
    author="Marianne Métais",
    url="https://github.com/marianne-m/brouhaha-vad",
    install_requires=requirements,
    extras_require={
        'optimization': [
            'cython>=0.29.0',
            'numba>=0.56.0',
        ],
        'dev': [
            'pytest>=6.0',
            'pytest-benchmark>=3.4.1',
            'cython>=0.29.0',
            'numba>=0.56.0',
        ],
    },
    ext_modules=extensions,
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering",
    ],
    zip_safe=False,  # Required for Cython extensions
)
