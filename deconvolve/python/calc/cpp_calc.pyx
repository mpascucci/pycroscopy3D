#!/usr/bin/env python3

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



import numpy as np
cimport numpy as np
cimport cython

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

FTYPE = np.float32
ctypedef np.float32_t FTYPE_t

from libcpp.vector cimport vector


cdef extern from "deconvolve.hpp" namespace "deconvolve":
    ctypedef int (*callbackfunc)(void *user_data, size_t iteration_number,
                                 double cmin, double cmax, double csum,
                                 double nrm2_prev, double nrm2_prevprev,
                                 double lmbda, double lambda_factor, double snr)
    cdef cppclass Deconvolve[T]:
        Deconvolve() except +
        void set_psf(const vector[T] &data, size_t n1, size_t n2, size_t n3, T v1, T v2, T v3)
        void set_callback(callbackfunc cb, void *user_data)
        void clear_callback()
        void enable_regularization();
        void disable_regularization();
        void set_snr(T snr)
        void clear_snr()
        void set_max_iterations(size_t iters)
        void clear_max_iterations()
        int regularized()
        vector[T] convolve(const vector[T] &data, size_t n1, size_t n2, size_t n3, T v1, T v2, T v3)
        vector[T] deconvolve(const vector[T] &data, size_t n1, size_t n2, size_t n3, T v1, T v2, T v3)

cdef int callback_for_deconvolution(void *f, size_t iteration_number, 
                                    double cmin, double cmax, double csum, 
                                    double nrm2_prev, double nrm2_prevprev, 
                                    double lmbda, double lambda_factor, double snr):
     return (<object>f)(iteration_number=iteration_number, cmin=cmin, cmax=cmax, csum=csum,
                        nrm2_prev=nrm2_prev, nrm2_prevprev=nrm2_prevprev, lmbda=lmbda, lmbda_factor=lambda_factor, snr=snr)

cdef class PyDeconvolve:
    cdef Deconvolve[double] *thisptr # hold a C++ instance which we're wrapping
    
    def __cinit__(self):
        '''
        Parameters:
        -----------

        '''
        self.thisptr = new Deconvolve[double]()

    def __dealloc__(self):
        del self.thisptr

    def enable_regularization(self):
        self.thisptr.enable_regularization()

    def disable_regularization(self):
        self.thisptr.disable_regularization()

    def regularized(self):
        return self.thisptr.regularized()

    def set_snr(self, v):
        self.thisptr.set_snr(v)

    def clear_snr(self):
        self.thisptr.clear_snr()

    def set_max_iterations(self, v):
        self.thisptr.set_max_iterations(v)

    def clear_max_iterations(self):
        self.thisptr.clear_max_iterations()

    cpdef void set_psf(self, np.ndarray[DTYPE_t, ndim=1, mode="c"] data, size_t n1, size_t n2, size_t n3, double v1, double v2, double v3):
        '''
        Parameters:
        -----------

        '''
        self.thisptr.set_psf(data, n1, n2, n3, v1, v2, v3)

    cpdef vector[double] convolve(self, np.ndarray[DTYPE_t, ndim=1, mode="c"] data, size_t n1, size_t n2, size_t n3, double v1, double v2, double v3):
        '''
        Parameters:
        -----------

        '''
        return self.thisptr.convolve(data, n1, n2, n3, v1, v2, v3)

    cpdef vector[double] deconvolve(self, np.ndarray[DTYPE_t, ndim=1, mode="c"] data, size_t n1, size_t n2, size_t n3, double v1, double v2, double v3, callback=None):
        '''
        Parameters:
        -----------

        '''

        if callback is None:
            self.thisptr.clear_callback()
        else:
            self.thisptr.set_callback(callback_for_deconvolution, <void*>callback)
            
        return self.thisptr.deconvolve(data, n1, n2, n3, v1, v2, v3)


# float
cdef class PyDeconvolveFloat:
    cdef Deconvolve[float] *thisptr # hold a C++ instance which we're wrapping
    
    def __cinit__(self):
        '''
        Parameters:
        -----------

        '''
        self.thisptr = new Deconvolve[float]()

    def __dealloc__(self):
        del self.thisptr

    def enable_regularization(self):
        self.thisptr.enable_regularization()

    def disable_regularization(self):
        self.thisptr.disable_regularization()

    def regularized(self):
        return self.thisptr.regularized()

    def set_snr(self, v):
        self.thisptr.set_snr(v)

    def clear_snr(self):
        self.thisptr.clear_snr()

    def set_max_iterations(self, v):
        self.thisptr.set_max_iterations(v)

    def clear_max_iterations(self):
        self.thisptr.clear_max_iterations()

    cpdef void set_psf(self, np.ndarray[FTYPE_t, ndim=1, mode="c"] data, size_t n1, size_t n2, size_t n3, double v1, double v2, double v3):
        '''
        Parameters:
        -----------

        '''
        self.thisptr.set_psf(data, n1, n2, n3, v1, v2, v3)

    cpdef vector[float] convolve(self, np.ndarray[FTYPE_t, ndim=1, mode="c"] data, size_t n1, size_t n2, size_t n3, double v1, double v2, double v3):
        '''
        Parameters:
        -----------

        '''
        return self.thisptr.convolve(data, n1, n2, n3, v1, v2, v3)

    cpdef vector[float] deconvolve(self, np.ndarray[FTYPE_t, ndim=1, mode="c"] data, size_t n1, size_t n2, size_t n3, double v1, double v2, double v3, callback=None):
        '''
        Parameters:
        -----------

        '''
        if callback is None:
            self.thisptr.clear_callback()
        else:
            self.thisptr.set_callback(callback_for_deconvolution, <void*>callback)
            
        return self.thisptr.deconvolve(data, n1, n2, n3, v1, v2, v3)
