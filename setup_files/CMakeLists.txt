# This file is needed to compile IOCBIO/Deconvolve
# and should be copied in the folder "deconvolve/cpp"

# Copyright: Marco Pascucci, 2022 (marpas.paris@gmail.com)
# distribued under GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007

project(iocbio_deconvolve)
cmake_minimum_required(VERSION 3.18)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include_directories(include)
add_compile_options(-Ofast -march=native -fPIC)

file(GLOB src "src/*")
add_library(iocbio_deconvolve ${src})

#target_link_libraries(iocbio_deconvolve -lpthread -lfftw3 -lfftw3f)


