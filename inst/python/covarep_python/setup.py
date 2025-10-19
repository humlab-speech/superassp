"""
Setup script for COVAREP Python implementation

Python reimplementation of COVAREP - Cooperative Voice Analysis Repository
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read long description from README
long_description = (Path(__file__).parent / "README.md").read_text()

setup(
    name='covarep',
    version='0.1.0',
    description='COVAREP - Cooperative Voice Analysis Repository (Python Implementation)',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='COVAREP Python Team',
    author_email='',
    url='https://github.com/covarep/covarep-python',
    packages=find_packages(),
    install_requires=[
        'numpy>=1.20.0',
        'scipy>=1.7.0',
        'soundfile>=0.10.0',
        'librosa>=0.9.0',
        'pysptk>=0.2.0',
        'numba>=0.56.0',
    ],
    extras_require={
        'dev': [
            'pytest>=6.0.0',
            'pytest-cov>=2.0.0',
            'jupyter>=1.0.0',
            'matplotlib>=3.3.0',
        ],
        'docs': [
            'sphinx>=4.0.0',
            'sphinx-rtd-theme>=1.0.0',
        ],
        'all': [
            'pyworld>=0.3.0',
            'resampy>=0.2.0',
            'cython>=0.29.0',
        ],
    },
    python_requires='>=3.8',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
    ],
    keywords='speech audio voice analysis glottal f0 pitch vocoder',
)
