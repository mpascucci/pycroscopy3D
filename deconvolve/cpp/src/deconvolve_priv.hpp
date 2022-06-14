
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


#ifndef IOCBIO_DECONVOLVE_PRIV_HPP
#define IOCBIO_DECONVOLVE_PRIV_HPP

#include "deconvolve.hpp"
#include "image.hpp"
#include "image_settings.hpp"
#include "psf.hpp"

#include <deque>
#include <memory>
#include <vector>


namespace deconvolve {

  ////////////////////////////////////////////////////////////////////////
  /// \brief Implementation of deconvolution interface API
  ///
  /// See class \ref Deconvolve for documentation of public methods.
  template <typename T>
  class DeconvolvePrivate {
  public:
    DeconvolvePrivate();
    ~DeconvolvePrivate();

    /// \brief Set point spread function for convolution and deconvolution operations
    void set_psf(const std::vector<T> &data, size_t n1, size_t n2, size_t n3, T v1, T v2, T v3);

    void set_callback(callback_cpp_type callback); ///< Set callback function.
    void clear_callback(); ///< Drop the specified callback.

    void enable_regularization();  ///< Use deconvolution with regularization
    void disable_regularization(); ///< Use deconvolution without regularization
    bool regularized() const { return m_regularize; } ///< Current state of regularization

    void set_snr(T snr); ///< Set SNR for the image
    void clear_snr();    ///< Estimate SNR by the default algorithm

    void set_max_iterations(size_t iters) { m_max_iterations = iters; } ///< Set maximal number of iterations.
    void clear_max_iterations() { m_max_iterations = const_max_iterations; } ///< Use default maximal number of iterations.
    size_t max_iterations() const { return m_max_iterations; } ///< Current number of maximal number of iterations.

    /// \brief Set FFTW plan handlers
    void set_fftw_handlers( const typename fftw_implementation<T>::plan_function &forward,
                            const typename fftw_implementation<T>::plan_function &inverse,
                            const typename fftw_implementation<T>::clear_function &clear );
    /// \brief Use default FFTW plan handlers
    void clear_fftw_handlers();

    /// \brief Convolve image
    void convolve(std::vector<T> &data, size_t n1, size_t n2, size_t n3, T v1, T v2, T v3);

    /// \brief Deconvolve image
    void deconvolve(std::vector<T> &data, size_t n1, size_t n2, size_t n3, T v1, T v2, T v3);

  protected:

    /// \brief Default callback for deconvolution
    ///
    /// This callback prints out iteration statistics on stdout and
    /// stops iterations either on reaching maximal iteration count
    /// (\sa set_max_iterations and \sa m_max_iterations) or if
    /// regularization factor lambda is decreasing for
    /// const_lambda_stack_size iterations
    ///
    int callback_default(size_t iteration_number,
                         double min, double max, double sum,
                         double nrm2_prev, double nrm2_prevprev,
                         double lambda, double lambda_factor, double snr);

  protected:

    const size_t const_lambda_stack_size{3}; ///< Default number of last lambda values to compare to
    const size_t const_max_iterations{100};  ///< Default maximal number of iterations

    std::shared_ptr< ImageSettings<T> > m_settings; ///< Current image settings

    PSF<T> m_psf; ///< Current PSF
    
    std::deque<T> m_lambda_evolution; ///< Used to track lambda changes during deconvolution by default callback

    callback_cpp_type m_callback; ///< Callback function specified by the user

    bool m_regularize{true}; ///< Whether to use regularization or not

    size_t m_max_iterations{const_max_iterations}; ///< Maximal number of iterations allowed by the default callback

    T m_snr{-1}; ///< Positive when specified by the user

  };
}

#endif
