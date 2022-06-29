# import setuptools
from setuptools import Extension, setup, find_packages
from Cython.Build import cythonize
import numpy
from glob import glob


VERSION_MAJOR = 0
VERSION_MINOR = 1
VERSION_MICRO = 6

with open("README.txt", "r") as fh:
    long_description = fh.read()

# just build the module with a previously compiled library
iocbio_deconvolve = [
    Extension("iocbio_deconvolve",
              language="c++",
              sources=[
                  "deconvolve/python/calc/cpp_calc.pyx",
                  *glob("deconvolve/cpp/src/*.cpp")],
              # include directories needed to compile the lib
              include_dirs=[
                  numpy.get_include(), "deconvolve/cpp/include"],
              # the shared librairies linked at import time
              libraries=['pthread', 'fftw3', 'fftw3f'],
              extra_compile_args=[
                  "-std=c++11", "-funroll-loops", "-Ofast", "-march=native", "-fPIC"],
              )
]

setup(name='Pycroscopy3D',
      version='{}.{}.{}'.format(
          VERSION_MAJOR, VERSION_MINOR, VERSION_MICRO),
      description='Python tools for microscopy',
      author='Marco Pascucci',
      long_description=long_description,
      long_description_content_type="text/markdown",
      author_email='marpas.paris@gmail.com',
      url='https://github.com/mpascucci/pycroscopy',
      ext_modules=cythonize(iocbio_deconvolve),
      packages=find_packages(),
      package_data={'': ['pycroscopy3D/Matlab/*']},
      include_package_data=True,
      install_requires=['numpy', 'matplotlib', 'tqdm', 'psutil',
                        'connected-components-3d', 'antspyx', 'multipagetiff'],
      entry_points={'console_scripts': [
          'pycro_register=pycroscopy3D.cli.registration:main',
          'pycro_unpad=pycroscopy3D.cli.unpad:main',
          'pycro_sum_stacks=pycroscopy3D.cli.sum_stacks:main',
          'pycro_deconvolve=pycroscopy3D.cli.deconvolution:main',
          'pycro_convert=pycroscopy3D.cli.ants_to_tif:main',
          'pycro_skew_correct_one=pycroscopy3D.cli.skew_correct_one:main',
          'pycro_skew_correct_many=pycroscopy3D.cli.skew_correct_many:main',
          'pycro_psf=pycroscopy3D.cli.psf_mean:main']
      },
      classifiers=[
          "Programming Language :: Python",
          "License :: OSI Approved :: GNU General Public License v2 (GPLv3)",
          "Operating System :: OS Independent",
      ]
      )
