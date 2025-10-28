"""
Setup script for Voxit Python module.

Installation:
    pip install -e .
    
With optimizations:
    pip install -e .[numba]
    pip install -e .[cython]
    pip install -e .[all]
"""

from setuptools import setup, find_packages
import os

# Read README
readme_path = os.path.join(os.path.dirname(__file__), 'README.md')
if os.path.exists(readme_path):
    with open(readme_path, 'r', encoding='utf-8') as f:
        long_description = f.read()
else:
    long_description = "Voxit: Voice and articulation complexity measures"

setup(
    name='voxit',
    version='1.0.0',
    description='Voice and articulation complexity measures for speech analysis',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Voxit Contributors',
    author_email='',
    url='https://github.com/[voxit-url]',
    packages=find_packages(),
    python_requires='>=3.7',
    
    install_requires=[
        'numpy>=1.19.0',
        'scipy>=1.5.0',
        'lempel_ziv_complexity>=0.2.0',
    ],
    
    extras_require={
        'numba': ['numba>=0.50.0'],
        'cython': ['Cython>=0.29.0'],
        'all': ['numba>=0.50.0', 'Cython>=0.29.0'],
        'dev': [
            'pytest>=6.0.0',
            'pytest-cov>=2.10.0',
            'flake8>=3.8.0',
        ],
    },
    
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Topic :: Multimedia :: Sound/Audio :: Speech',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'License :: OSI Approved :: MIT License',
    ],
    
    keywords='speech prosody complexity voice articulation acoustics',
    
    project_urls={
        'Documentation': 'https://[docs-url]',
        'Source': 'https://github.com/[voxit-url]',
        'Bug Reports': 'https://github.com/[voxit-url]/issues',
    },
    
    include_package_data=True,
    zip_safe=False,
)
