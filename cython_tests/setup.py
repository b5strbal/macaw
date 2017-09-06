from distutils.core import setup
from Cython.Build import cythonize

setup(ext_modules = cythonize(
           "arrays.pyx",                 # our Cython source
    sources=["arrays.cpp"],  # additional source file(s)
           language="c++",             # generate C++ code
      ))
