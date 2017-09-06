#!/usr/bin/env bash
r"""
Installation script for the Macaw module
It depends on distutils
"""

try:
    from sage.env import SAGE_SRC
except ImportError:
    raise ValueError("this package currently installs only inside SageMath (http://www.sagemath.org)")

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import sys, os


setup(name='Macaw',
      version='0.23',
      description="Macaw",
      author='Balazs Strenner',
      author_email='strenner@math.gatech.edu',
      url='http://www.math.gatech.edu/~bstrenner7',
      license="GPL v3",
      packages=['macaw',
                'macaw/train_tracks',
                'macaw/train_tracks/dehn_thurston'],
      classifiers=[
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)',
          'Operating System :: OS Independent',
          'Programming Language :: Python',
          'Topic :: Scientific/Engineering :: Mathematics',
      ],
      keywords='surfaces, mapping class groups, train tracks, measured laminations',
)
