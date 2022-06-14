
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


#ifndef IOCBIO_DECONVOLVE_HPP
#define IOCBIO_DECONVOLVE_HPP

#include "fftw_interface.hpp"

#include <vector>
#include <memory>
#include <functional>
#include <stddef.h>

namespace deconvolve {

  template <typename T> class DeconvolvePrivate;

  /// \brief Callback function type used to communicate during deconvolution. For use from C
  ///
  /// Callback function is called during each iteration and should
  /// return zero if the iteration should be stopped. To continue
  /// iteration, return non-zero value, for example 1.
  ///
  /// Note that the callback function is called first before the first
  /// iteration.
  ///
  /// Function arguments:
  ///
  /// \param iteration_number Current iteration (from zero onwards)
  /// \param min minimal value of deconvolved image. Set for iterations with `iteration_number` 1 and larger
  /// \param max maximal value of deconvolved image. Set for iterations with `iteration_number` 1 and larger
  /// \param sum sum of all values of deconvolved image. Set for iterations with `iteration_number` 1 and larger
  /// \param nrm2_prev Euclidean norm between current and previous estimate of deconvolved image. Set for iterations with `iteration_number` 1 and larger
  /// \param nrm2_prevprev Euclidean norm between current and estimate from two steps before of deconvolved image. Set for iterations with `iteration_number` 2 and larger
  /// \param lambda current lambda value. Set for iterations with `iteration_number` 1 and larger
  /// \param lambda_factor factor used in the calculations of lambda. Corresponds to "C" in equation 5 of https://doi.org/10.1111/j.1365-2818.2011.03486.x .Set for iterations with `iteration_number` 1 and larger 
  /// \param snr used SNR estimate for original image. If not set by user, SNR is calculated before iterations. Set for all calls.
  /// \return Zero to stop iterations, continue iterations for non-zero values
  /// 
  typedef int (*callback_type)(size_t iteration_number,
                               double min, double max, double sum,
                               double nrm2_prev, double nrm2_prevprev,
                               double lambda, double lambda_factor, double snr);
  
  /// \brief Extended callback function type used to communicate during deconvolution
  ///
  /// This is an extension of \ref callback_type to include user-supplied data pointer
  /// that is passed on every call. All other parameters and return values are the
  /// same as for \ref callback_type
  ///
  /// \sa callback_type
  typedef int (*callback_extended_type)(void *user_data, size_t iteration_number,
                                        double min, double max, double sum,
                                        double nrm2_prev, double nrm2_prevprev,
                                        double lambda, double lambda_factor, double snr);
  
  /// \brief Callback function type used to communicate during deconvolution. For use from C++
  ///
  typedef std::function<int(size_t iteration_number,
                            double min, double max, double sum,
                            double nrm2_prev, double nrm2_prevprev,
                            double lambda, double lambda_factor, double snr)> callback_cpp_type;
  
  ///////////////////////////////////////////////////////////////////////////////////////////////
  
  /// \brief Deconvolution and convolution of 3D images
  ///
  /// This is an externally accessible class that encapsulates
  /// functionality of the package. It is provided as a template class
  /// which is instantiated in `double` or `float` precision. Select
  /// the desired precision and use accordingly.
  ///
  /// For deconvolving or convolving images, you have to specify point
  /// spread function first using \ref set_psf method. Next, to
  /// control deconvolution, it is advisable to specify the callback
  /// function using one of the \ref set_callback overloaded
  /// versions. After that, images can be deconvolved using \ref
  /// deconvolve. For convolution, only point spread function is
  /// needed and convolution cal be performed using \ref convolve
  /// method.
  ///
  /// It is possible to specify FFTW plan handlers to incorporate FFTW
  /// wisdom management in the application. See \ref set_fftw_handlers for
  /// the description of the corresponding API.
  ///
  /// Note that all voxel sizes are expected in / meters. Image is
  /// expected to have Poisson noise.
  ///
  /// This template class is available either for double or float types.
  ///
  /// \sa deconvolve::callback_type
  /// \sa set_psf
  /// \sa set_fftw_handlers
  ///
  template <typename T>
  class Deconvolve {
  public:

    Deconvolve();
    ~Deconvolve();

    /// \brief Set point spread function for convolution and deconvolution operations
    ///
    /// Sets point spread function that will be used when performing
    /// convolution or deconvolution. Point spread function can be set
    /// many times with the last value used in the following
    /// operations. Note that point spread function has to be
    /// specified before \ref convolve or \ref deconvolve is called.
    ///
    /// \param data vector of size n1*n2*n3 with the last dimension (n3) changing the fastest
    /// \param n1 the slowest changing dimension
    /// \param n2 the medium changing dimension
    /// \param n3 the fastest changing dimension
    /// \param v1 voxel size along dimension 1, in meters
    /// \param v2 voxel size along dimension 2, in meters
    /// \param v3 voxel size along dimension 3, in meters
    void set_psf(const std::vector<T> &data, size_t n1, size_t n2, size_t n3, T v1, T v2, T v3);

    /// \brief Set callback function to communicate during deconvolution
    ///
    /// Callback function can be specified to provide user interface
    /// for deconvolution and determine whether to stop deconvolution
    /// process. When the callback function is specified, the default
    /// console messages will not be printed. To restore the default
    /// callback, call \ref clear_callback
    ///
    /// See \ref callback_type for description of the callback
    /// function.
    ///
    /// \sa clear_callback
    ///
    void set_callback(callback_type callback);

    /// \brief Set callback function to communicate during deconvolution
    ///
    /// Overloaded version. Here, a version of callback function more
    /// suitable for C++ (\ref callback_cpp_type) is used.
    ///
    /// \sa set_callback
    ///
    void set_callback(callback_cpp_type callback);

    /// \brief Set callback function to communicate during deconvolution
    ///
    /// Overloaded version. Here, a version with additional user data passed
    /// as an argument is used.
    ///
    /// \sa set_callback
    /// \sa callback_extended_type
    /// 
    void set_callback(callback_extended_type callback, void *user_data);

    /// \brief Drop the specified callback and use the default one
    ///
    /// Drops the callback previously specified by one of the
    /// \ref set_callback functions and will use the default internal
    /// callback with console messages and stopping criteria. It is
    /// safe to call this method even if set_callback was not called
    /// before.
    ///
    /// \sa set_callback
    ///
    void clear_callback();

    /// \brief Use deconvolution with regularization
    ///
    /// Use deconvolution algorithm with regularization
    /// component. This is the default setting.
    void enable_regularization();

    /// \brief Use deconvolution without regularization
    ///
    /// Use deconvolution algorithm without regularization
    /// component.
    void disable_regularization();

    /// \brief Current state of regularization
    ///
    /// \return `true` if regularized deconvolution algorithm is set to be used
    ///
    bool regularized() const;

    /// \brief Set SNR for the image
    ///
    /// Sets SNR that will be used to estimate lambda factor at the
    /// beginning of regularized deconvolution process. If not set,
    /// SNR will be estimated assuming that the image data represents
    /// counts which are not scaled by gain. See "Image requirements"
    /// in the main [README](index.html) of the library.
    ///
    /// \param snr SNR for the image. Should be positive
    ///
    void set_snr(T snr);
    
    /// \brief Clear previously set SNR
    ///
    /// Clears SNR estimate set by \ref set_snr and will use the
    /// default estimation for SNR assuming that the image represents
    /// photon counts
    void clear_snr();

    /// \brief Set maximal number of iterations
    ///
    /// Sets maximal number of iterations for the following
    /// deconvolution. Note that this iteration limit is used only if
    /// the callback is not specified. For regularized deconvolution,
    /// iterations will be stopped by the default callback either if
    /// this number of iterations is exceeded or regularization factor
    /// will start to decline.
    ///
    /// \param iters Maximal number of iterations
    ///
    void set_max_iterations(size_t iters);

    /// \brief Clear previously set maximal number of iterations
    ///
    /// Clears iteration count limit set by \ref set_max_iterations
    /// and use the compiled in default
    ///
    /// \sa set_max_iterations
    ///
    void clear_max_iterations();

    /// \brief Current maximal number of iterations
    ///
    /// Maximal number of iterations allowed by the default callback
    /// function.
    ///
    /// \return maximal number of iterations
    ///
    size_t max_iterations() const;
    
    /// \brief Set FFTW plan handling functions
    ///
    /// Set functions that allocate FFTW plans and destroy them. This
    /// facility allows to implement a central handling of FFTW plan
    /// construction and destruction through user-provided
    /// functions. While the library has thread-safe plan handling
    /// functions, it maybe feasible to provide plan handling by
    /// user's functions allowing users to handle saving and loading
    /// of the wisdom. To revert to the default FFTW plan handling,
    /// call \ref clear_fftw_handlers.
    ///
    /// There are three functions that need to be provided as
    /// handlers: two to create forward and inverse plans and one to
    /// destroy them. When it is needed to provide plan handling
    /// through object methods, use `std::bind` to create
    /// corresponding `std::function` object.
    ///
    /// All plans are designed to be used for 3D DFTs of
    /// real data with data input and output stored in the same memory
    /// array. Note that the memory allocation is handled by the
    /// library itself, so there is no need to be concerned about
    /// it. The plans will be constructed and destroyed in accordance
    /// with the memory allocation and release.
    /// 
    /// Forward FFTW plan should be created by a function that takes
    /// the pointer to the data structure and its dimensions as an
    /// integer. The created plan should use the data for input and
    /// output. An example implementation for double precision
    /// assuming single threaded application would be
    ///
    ///     fftw_plan forward(double *data, int n0, int n1, int n2) {
    ///         return fftw_plan_dft_r2c_3d(n0, n1, n2, data, (fftw_complex*)data, FFTW_ESTIMATE);
    ///     }
    ///
    /// In this example, no synchronization primitives were used. Note
    /// that, in multi-threaded application, you would want to protect
    /// plan handlers with some synchronization primitive, such as a
    /// `std::mutex`.
    ///
    /// Inverse FFTW plan uses the same approach as the forward plan:
    /// data pointer is used for input and output. An example
    /// implementation for double assuming single threaded application
    /// would be
    ///
    ///     fftw_plan inverse(double *data, int n0, int n1, int n2) {
    ///         return fftw_plan_dft_c2r_3d(n0, n1, n2, (fftw_complex*)data, data, FFTW_ESTIMATE);
    ///     }
    ///
    /// Finally, the function for destruction of the plan has to be
    /// provided. This function takes the plan pointer as an argument
    /// and destroys it. For example,
    ///
    ///     void clear(fftw_plan plan) {
    ///         fftw_destroy_plan(plan);
    ///     }
    ///
    /// \param forward Function to create FFTW plan for 3D forward DFT
    /// of real data. By default, it is handled by
    /// `fftw_plan_dft_r2c_3d` or `fftwf_plan_dft_r2c_3d` through
    /// wrappers.
    ///
    /// \param inverse Function to create FFTW plan for 3D inverse DFT
    /// to the real data. By default, it is handled by
    /// `fftw_plan_dft_c2r_3d` or `fftwf_plan_dft_c2r_3d` through
    /// wrappers.
    ///
    /// \param clear Function to destroy FFTW plan. By default, it is
    /// handled by `fftw_destroy_plan` or `fftwf_destroy_plan` through
    /// wrappers.
    ///
    /// \sa clear_fftw_handlers
    ///
    void set_fftw_handlers( const typename fftw_implementation<T>::plan_function &forward,
                            const typename fftw_implementation<T>::plan_function &inverse,
                            const typename fftw_implementation<T>::clear_function &clear );

    /// \brief Reset FFTW handlers to the library-provided defaults
    ///
    /// Use to restore the library-provided defaults for FFTW
    /// handlers.
    ///
    /// \sa set_fftw_handlers
    ///
    void clear_fftw_handlers();
    
    /// \brief Convolve image with the point spread function
    ///
    /// Convolves given image with the point spread function specified
    /// earlier by \ref set_psf.
    ///
    /// \param data vector of size n1*n2*n3 with the last dimension (n3) changing the fastest
    /// \param n1 the slowest changing dimension
    /// \param n2 the medium changing dimension
    /// \param n3 the fastest changing dimension
    /// \param v1 voxel size along dimension 1, in meters
    /// \param v2 voxel size along dimension 2, in meters
    /// \param v3 voxel size along dimension 3, in meters
    ///
    /// \return convolved image in the same format as data
    ///
    std::vector<T> convolve(const std::vector<T> &data, size_t n1, size_t n2, size_t n3, T v1, T v2, T v3);
    
    /// \brief Deconvolve image
    ///
    /// Deconvolves given image by taking into account the point
    /// spread function specified earlier by \ref set_psf.
    ///
    /// \param data vector of size n1*n2*n3 with the last dimension (n3) changing the fastest
    /// \param n1 the slowest changing dimension
    /// \param n2 the medium changing dimension
    /// \param n3 the fastest changing dimension
    /// \param v1 voxel size along dimension 1, in meters
    /// \param v2 voxel size along dimension 2, in meters
    /// \param v3 voxel size along dimension 3, in meters
    ///
    /// \return deconvolved image in the same format as data
    ///
    std::vector<T> deconvolve(const std::vector<T> &data, size_t n1, size_t n2, size_t n3, T v1, T v2, T v3);

  private:

    std::unique_ptr<DeconvolvePrivate<T> > m_dec; ///< Pointer to the private implementation of deconvolution class
  };
}

#endif
