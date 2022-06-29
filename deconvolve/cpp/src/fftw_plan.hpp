
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


#ifndef IOCBIO_FFTWPLAN_HPP
#define IOCBIO_FFTWPLAN_HPP

#include "fftw_interface.hpp"
#include "image_settings.hpp"

#include <fftw3.h>

#include <functional>
#include <memory>

namespace deconvolve {

  /// \brief Generalization of FFTW types and plan functions
  ///
  /// See fftw_implementation_detail<double> and fftw_implementation_detail<float> for
  /// used members.
  ///
  template <typename T>
  struct fftw_implementation_detail {
  };

  /// \brief Default FFTW handlers for `double` precision
  ///
  /// Note that these handlers use directly functions from FFTW. This
  /// is in contrast to the handlers that can be specified by
  /// Deconvolve::set_fftw_handlers . Such difference is made to
  /// simplify the library implementation.
  template<>  
  struct fftw_implementation_detail<double> {
    /// \brief Default FFTW handler for the forward plan generation.
    const std::function<fftw_implementation<double>::plan_type(int n0, int n1, int n2,
                                                              double *in, fftw_implementation<double>::complex_type *out,
                                                              unsigned flags)> forward{fftw_plan_dft_r2c_3d};
    /// \brief Default FFTW handler for the inverse plan generation.
    const std::function<fftw_implementation<double>::plan_type(int n0, int n1, int n2,
                                                              fftw_implementation<double>::complex_type *in, double *out,
                                                              unsigned flags)> inverse{fftw_plan_dft_c2r_3d};
    /// \brief Default function for FFTW plan destruction.
    const std::function<void(fftw_implementation<double>::plan_type)> clear{fftw_destroy_plan};
  };

  template<>
  struct fftw_implementation_detail<float> {
    /// \brief Default FFTW handler for the forward plan generation.
    const std::function<fftw_implementation<float>::plan_type(int n0, int n1, int n2,
                                                             float *in, fftw_implementation<float>::complex_type *out,
                                                             unsigned flags)> forward{fftwf_plan_dft_r2c_3d};
    /// \brief Default FFTW handler for the inverse plan generation.
    const std::function<fftw_implementation<float>::plan_type(int n0, int n1, int n2,
                                                             fftw_implementation<float>::complex_type *in, float *out,
                                                             unsigned flags)> inverse{fftwf_plan_dft_c2r_3d};
    /// \brief Default function for FFTW plan destruction.
    const std::function<void(fftw_implementation<float>::plan_type)> clear{fftwf_destroy_plan};
  };

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  
  /// \brief Wrapper around FFTW plan for real-to-comlpex and complex-to-real in-place 3D transformations
  ///
  /// This wrapper can keep either forward or inverse FFTW plan. The
  /// plan is allocated using \ref forward or \ref inverse
  /// methods. After allocation, the plan can be used using \ref
  /// execute method. On destruction or calling \ref clean, the
  /// allocated plan is destroyed.
  template <typename T>
  class FFTWPlan {
  public:
    /// \brief Construct a wrapper around FFTW plans using FFTW handlers as specified in ImageSettings
    FFTWPlan(std::shared_ptr< ImageSettings<T> > settings);
    FFTWPlan() = delete;
    
    ~FFTWPlan();

    /// \brief Create FFTW plan for forward FFT
    void forward(T *data, int n0, int n1, int n2);

    /// \brief Create FFTW plan for inverse FFT
    void inverse(T *data, int n0, int n1, int n2);

    /// \brief Destroy FFTW plan
    void clear();

    /// \brief Execute FFTW plan
    ///
    /// Note that the plan has to be created first using \ref forward
    /// or \ref inverse methods.
    void execute();

    /// \brief Check if the plan is allocated
    ///
    /// \return `true` if the plan is allocated
    operator bool() const { return m_plan != NULL; }

  protected:
    typename fftw_implementation<T>::plan_type m_plan = NULL; ///< FFTW plan
    std::shared_ptr< ImageSettings<T> > m_settings;          ///< ImageSettings specifying FFTW plan handler functions 
  };
}


#endif

