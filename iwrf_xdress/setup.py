import os
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy as np

incdirs = [os.path.join(os.getcwd(), '..','include'),
           os.path.join(os.getcwd(), '..','rsm'),
           ]

ext_modules = [
 Extension("iwrf.iwrf_functions", ["iwrf/iwrf_functions.pyx"],
           include_dirs=incdirs, language="c"),
 Extension("iwrf.c_iwrf_functions", ["iwrf/c_iwrf_functions.pyx"],
           include_dirs=incdirs, language="c"),
]

setup(name="iwrf",
      cmdclass={'build_ext': build_ext},
    ext_modules = ext_modules,
    packages = ['iwrf']
)