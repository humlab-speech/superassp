"""
Setup script for formant tracking with optional Cython acceleration.

Install regular version:
    pip install -e .

Build with Cython acceleration:
    pip install cython
    python setup.py build_ext --inplace
"""

from setuptools import setup, find_packages, Extension
import numpy as np

# Try to import Cython
try:
    from Cython.Build import cythonize
    USE_CYTHON = True
except ImportError:
    USE_CYTHON = False

# Define Cython extensions (only if Cython is available)
ext_modules = []
if USE_CYTHON:
    extensions = [
        Extension(
            "ftrack.gloat.pitch_cython",
            ["ftrack/gloat/pitch_cython.pyx"],
            include_dirs=[np.get_include()],
            extra_compile_args=["-O3", "-ffast-math"],
            define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
        ),
        Extension(
            "ftrack.gloat.gci_cython",
            ["ftrack/gloat/gci_cython.pyx"],
            include_dirs=[np.get_include()],
            extra_compile_args=["-O3", "-ffast-math"],
            define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
        ),
        Extension(
            "ftrack.gloat.utils_cython",
            ["ftrack/gloat/utils_cython.pyx"],
            include_dirs=[np.get_include()],
            extra_compile_args=["-O3", "-ffast-math"],
            define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
        ),
    ]
    ext_modules = cythonize(
        extensions,
        compiler_directives={
            'language_level': '3',
            'boundscheck': False,
            'wraparound': False,
            'cdivision': True,
            'initializedcheck': False,
        }
    )

setup(
    name="ftrack",
    version="0.1.0",
    description="Python implementation of TVWLP formant tracking (with optional Cython)",
    author="Python port from MATLAB implementation by D.Gowda",
    packages=find_packages(),
    ext_modules=ext_modules,
    install_requires=[
        "numpy>=1.20.0",
        "scipy>=1.7.0",
        "matplotlib>=3.3.0",
        "librosa>=0.9.0",
        "numba>=0.57.0",  # For JIT acceleration
    ],
    extras_require={
        "cython": ["cython>=0.29.0"],
    },
    python_requires=">=3.8",
    zip_safe=False,
)
