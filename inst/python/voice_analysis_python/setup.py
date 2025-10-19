"""
Voice Analysis Toolbox - Python Implementation
A faithful reimplementation of the MATLAB Voice Analysis Toolbox
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="voice-analysis-toolbox",
    version="0.1.0",
    author="Voice Analysis Team",
    author_email="",
    description="Comprehensive voice analysis computing 132 dysphonia measures",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.21.0",
        "scipy>=1.7.0",
        "soundfile>=0.10.0",
        "librosa>=0.9.0",
        "pywt>=1.1.1",
        "pysptk>=0.1.0",
        "nolds>=0.5.0",
        "PyEMD>=1.3.0",
        "numba>=0.54.0",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.0",
            "black>=21.0",
            "flake8>=3.9",
            "sphinx>=4.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "voice-analysis=voice_analysis.cli:main",
        ],
    },
)
