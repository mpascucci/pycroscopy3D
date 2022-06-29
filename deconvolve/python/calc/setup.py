
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#  
#   Copyright (C) 2018-2020
#    Laboratory of Systems Biology, Department of Cybernetics,
#    School of Science, Tallinn University of Technology
#   This file is part of project: IOCBIO Deconvolve


import numpy
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

# setup(
#     ext_modules = [cythonize('c_model_single_exp.pyx')],
#     include_dirs = [numpy.get_include()],
#     libraries=["m"],
# )

ext_modules=[
    Extension("cpp_calc",
              sources=["cpp_calc.pyx"],
              libraries=["deconvolve"], # Unix-like specific # "gfortran",
              language="c++",
              library_dirs = ["../../cpp"],
              include_dirs = [numpy.get_include(), "../../cpp/include"],
              extra_compile_args = ["-std=c++11", "-funroll-loops"],
          )
]

setup(
    ext_modules = cythonize(ext_modules)
)

# run as:
# python3 setup.py build_ext --inplace
# creating html version of generated c++ file
# cython --cplus -a c_model_single_exp.pyx
