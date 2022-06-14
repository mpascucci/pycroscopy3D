
/*
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *  
 *   Copyright (C) 2018-2020
 *    Laboratory of Systems Biology, Department of Cybernetics,
 *    School of Science, Tallinn University of Technology
 *   This file is part of project: IOCBIO Deconvolve
 */


#ifndef IOCBIO_FFTWINTERFACE_HPP
#define IOCBIO_FFTWINTERFACE_HPP

#include <fftw3.h>

#include <functional>

// Types used in FFTW interface 

namespace deconvolve {

  /// \brief Generalization of FFTW types and plan functions, empty structure
  ///
  /// See fftw_implementation<double> and fftw_implementation<float> for
  /// used types.
  ///
  template <typename T>
  struct fftw_implementation {
  };

  /// \brief Specialization of FFTW types and plan functions for `double` precision
  ///
  template<>
  struct fftw_implementation<double> {    
    typedef fftw_plan plan_type;       ///< Type of FFTW plan
    typedef fftw_complex complex_type; ///< Type for complex numbers

    /// \brief Function type used to create forward and inverse FFTW plans
    ///
    typedef std::function<plan_type(double *data, int n0, int n1, int n2)> plan_function;

    /// \brief Function type for clearing (destroying) FFTW plan
    ///
    typedef std::function<void(plan_type)> clear_function;
  };

  /// \brief Specialization of FFTW types and plan functions for `float` precision
  ///
  template<>
  struct fftw_implementation<float> {
    typedef fftwf_plan plan_type;       ///< Type of FFTW plan
    typedef fftwf_complex complex_type; ///< Type for complex numbers

    /// \brief Function type used to create forward and inverse FFTW plans
    ///
    typedef std::function<plan_type(float *data, int n0, int n1, int n2)> plan_function;

    /// \brief Function type for clearing (destroying) FFTW plan
    ///
    typedef std::function<void(plan_type)> clear_function;
  };
}

#endif
