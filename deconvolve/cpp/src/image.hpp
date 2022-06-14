
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


#ifndef IOCBIO_IMAGE_HPP
#define IOCBIO_IMAGE_HPP

#include "constants.hpp"
#include "fftw_plan.hpp"
#include "image_settings.hpp"

#include <array>
#include <complex>
#include <memory>
#include <vector>

#include <stddef.h>

namespace deconvolve {

  /// \brief Image class defining memory allocation and all mathematical operations required for deconvolution
  ///
  /// This is a core class that defines memory allocation, operations
  /// required for deconvolution and interaction with the
  /// user-provided ImageSettings. This class is responsible for
  /// taking the user data in the form of continuous data arrays,
  /// storing it in the format that is suitable for the used internal
  /// functions, perform mathematical operations and, on request,
  /// return the data to the user as a continuous data array.
  ///
  /// The current implementation is based on FFTW library and uses
  /// FFTW representation for real and complex data (FFTW is referred
  /// to as a backend below). In future, if there will be more
  /// approaches to implement deconvolution operations then this class
  /// should be separated into an abstract interface, FFTW-based Image
  /// class and the classes defining other approaches. For example,
  /// deconvolution using GPUs is expected to be implemented in such a
  /// way.
  ///
  template <typename T>
  class Image {
  public:

    /// \brief Constructor that copies user-provided data
    ///
    /// Use this constructor to allocate memory in backend-supported
    /// format and copy the data from user-provided vector.
    ///
    /// \param settings image settings describing memory storage options and operations
    /// \param data image data vector of size n1*n2*n3 to be copied
    /// \param n1 the slowest changing dimension
    /// \param n2 the medium changing dimension
    /// \param n3 the fastest changing dimension
    /// \param v1 voxel size along dimension 1, in meters
    /// \param v2 voxel size along dimension 2, in meters
    /// \param v3 voxel size along dimension 3, in meters
    ///
    Image(std::shared_ptr< ImageSettings<T> > settings, const std::vector<T> &data, size_t n1, size_t n2, size_t n3, T v1, T v2, T v3);

    /// \brief Constructor that only allocates memory
    ///
    /// Use this constructor to allocate memory in backend-supported
    /// format. The constructed image data is not initialized and can
    /// be used later to store results of mathematical operations on
    /// other images.
    ///
    /// \param settings image settings describing memory storage options and operations
    /// \param n1 the slowest changing dimension
    /// \param n2 the medium changing dimension
    /// \param n3 the fastest changing dimension
    /// \param v1 voxel size along dimension 1, in meters
    /// \param v2 voxel size along dimension 2, in meters
    /// \param v3 voxel size along dimension 3, in meters
    ///
    Image(std::shared_ptr< ImageSettings<T> > settings, size_t n1, size_t n2, size_t n3, T v1, T v2, T v3);

    ~Image();

    // delete default and copy constructors and operator
    Image() = delete;
    Image(const Image&) = delete;
    Image& operator=(const Image&) = delete;

    /// \brief Test if image has allocated memory for data
    ///
    /// \return `true` if memory has been allocated
    ///
    operator bool() const { return m_data != nullptr;}

    /// \brief Set image data
    ///
    /// Makes a copy of the data in backend-supported internal format
    ///
    /// \param data image data vector of size n1*n2*n3 to be copied
    /// \param n1 the slowest changing dimension
    /// \param n2 the medium changing dimension
    /// \param n3 the fastest changing dimension
    /// \param v1 voxel size along dimension 1, in meters
    /// \param v2 voxel size along dimension 2, in meters
    /// \param v3 voxel size along dimension 3, in meters
    ///
    void set(const std::vector<T> &data, size_t n1, size_t n2, size_t n3, T v1, T v2, T v3);

    /// \brief Copy data from image
    ///
    /// Makes a copy of the data into `this` from the provided
    /// image. Note that the images have to be compatible (same
    /// dimensions and voxel sizes).
    ///
    /// \param im image to copy the data from
    ///
    void copy_data(const Image& im);

    /// \brief Get image as a continuous data vector 
    ///
    /// Use to get resulting image assuming that stored `this->m_data`
    /// is real. The given vector is resized to fill the data and then
    /// filled.
    ///
    /// \param data vector to fill the data to
    ///
    void get_image(std::vector<T> &data);

    /// \brief Swaps images between this and provided image
    ///
    /// Swap image data and corresponding structures between this and
    /// image given as an argument
    ///
    /// \param image Image to be swapped with
    void swap(Image &image);

    /// \brief Compare given dimensions with the dimensions of `this` image
    ///
    /// Use to compare `this` image with the dimensions given as arguments.
    ///
    /// \return `true` if the image has the same number of voxels in all dimensions
    ///
    bool same_dims(size_t n1, size_t n2, size_t n3) const;
    
    /// \brief Compare dimensions of the other image with the dimensions of `this` image
    ///
    /// Overloaded version of `same_dims`. Use to compare `this` image
    /// with the dimensions of the image given as an argument.
    ///
    /// \return `true` if the image has the same number of voxels in all dimensions as the other image
    ///
    bool same_dims(const Image &other) const;

    
    /// \brief Compare given voxel sizes with the sizes of `this` image
    ///
    /// Use to compare `this` image voxel sizes with the sizes given as arguments.
    ///
    /// \return `true` if the image has the same voxel sizes in all dimensions
    ///
    bool same_voxel(T v1, T v2, T v3) const;
    
    /// \brief Compare voxel sizes of the other image with the sizes of `this` image
    ///
    /// Overloaded version of `same_voxel`. Use to compare `this`
    /// image voxel sizes with the sizes of the image given as an
    /// argument.
    ///
    /// \return `true` if the image has the same voxel sizes in all dimensions as the other image
    ///
    bool same_voxel(const Image &other) const;

    /// \brief Compare if the image settings are the same as the settings given by an argument
    ///
    /// \param settings Settings to be compared to
    /// \return `true` if the settings are identical
    ///
    bool same_settings(const std::shared_ptr< ImageSettings<T> > settings) { return m_settings && settings && m_settings->same(settings); }

    /// \brief Tests if other image is compatible to `this` image
    ///
    /// Checks if `this` and the other image have the same dimensions,
    /// voxels, and have both allocated memory
    ///
    /// \return `true` if the memory was allocated for the both images and voxels and dimensions are the same.
    ///
    bool compatible(const Image &other) const;
    
    /// \brief Applies FFT to this image and keeps it in place
    ///
    /// It it assumed that `this` holds real data. As a result of this
    /// operation, `this` will be in Fourier space
    void fft();

    /// \brief Applies IFFT to this image
    ///
    /// It is assumed that `this` is in Fourier space. After `ifft`, real
    /// data is recovered into `this`.
    void ifft();

    /// \brief Convolve this image by the given kernel
    ///
    /// It is assumed that `this` holds real image and the image given
    /// as a kernel is in Fourier space. In addition, the kernel and
    /// this image are expected to be \ref compatible. After
    /// operation, `this` holds the real data corresponding to the
    /// convolved image.
    void convolve(const Image &kernel);

    /// \brief Convolve this image by conjugate of the given kernel
    ///
    /// It is assumed that `this` holds real image and the image given
    /// as a kernel is in Fourier space. In addition, the kernel and
    /// this image are expected to be \ref compatible. After
    /// operation, `this` holds the real data corresponding to the
    /// convolved image.
    void convolve_conj(const Image &kernel);

    /// \brief Operation used in deconvolution: this=image/this
    ///
    /// The operation takes value of each voxel in `image`, divides it
    /// by the corresponding voxel value in `this` and stores the
    /// result as a new value for voxel in `this`. Images (`this` and
    /// `image`) are expected to be real and \ref compatible.
    void invdivide_image(const Image &image);

    /// \brief Operation used in deconvolution: this=image*this
    ///
    /// The operation takes value of each voxel in `image`, multiplies it
    /// by the corresponding voxel value in `this` and stores the
    /// result as a new value for voxel in `this`. Images (`this` and
    /// `image`) are expected to be real and \ref compatible.
    void prod_image(const Image &image);

    /// \brief Operation used in deconvolution: this=div(grad(image)/mod(grad(image)))
    ///
    /// Calculation of divergence of the normalized gradient of the
    /// image given as an argument. The result is stored in `this`. It
    /// is expected that `image` is \ref compatible with `this` and
    /// `image` stores real data.
    void div_unit_grad(const Image &image);

    /// \brief Operation used in deconvolution: this=this*(image/(1-lambda*div))
    ///
    /// Calculation of operation used in regularized
    /// deconvolution. Stores result in `this`. It is expected that
    /// `image` is \ref compatible with `this` and `image` stores real data.
    void prod_regularized(const Image &image, T lambda, const Image &div);

    /// \brief Peak signal-to-noise of this image
    ///
    /// Find the peak signal-to-noise ratio assuming that the image
    /// corresponds to the recordings of Poisson process. This
    /// implementation finds the maximal average intensity for `this`
    /// image in a small box with the given kernel size. From the
    /// average value, signal-to-noise can be estimated assuming that
    /// the recordings are in accordance with the Poisson process.
    ///
    /// \param convolution_kernel_size single edge of the box used to find the local average
    /// \return estimated signal-to-noise for `this` image
    ///
    T snr(size_t convolution_kernel_size) const;

    /// \brief Image characteristics in terms of extreme values and sum
    ///
    /// Find minimum, maximum, and the sum of the voxel values in `this` image
    ///
    /// \param `cmin` on return, minimum value
    /// \param `cmax` on return, maximum value
    /// \param `csum` on return, sum of all voxel values
    ///
    void get_stats(T &cmin, T &cmax, T &csum) const;

    /// \brief Euclidean norm between `this` and provided image
    ///
    /// Estimate the difference between `this` and the image given as
    /// an argument using Euclidean norm.
    ///
    /// \return Euclidean norm between this and provided image
    T nrm2(const Image &image) const;

    /// \brief Lambda LSQ calculation
    ///
    /// Calculate Lambda LSQ (Equation 5 in [Laasmaa et
    /// al](https://doi.org/10.1111/j.1365-2818.2011.03486.x))
    ///
    /// \param cconv Image containing `i/(o[s] x h) x conj(h)` from the formula
    /// \param div Image containing divergence of the normalized gradient of `o[s]`
    /// \return Lambda LSQ
    ///
    static T lambda_lsq(const Image &cconv, const Image &div);

  protected:
    void allocate_data(); ///< Allocate data in the format suitable for the backend
    void release_data();  ///< Release the allocated data

    /// \brief Implementation of convolution
    ///
    /// Implementation of generalized convolution allowing to convolve
    /// `this` with the kernel or its conjugated form.
    void convolve_implementation( const Image<T> &kernel,
                                  void (*p)(std::complex<T> *im, std::complex<T> *ker, T scale, size_t sz) );
    
    size_t last_dim() const;  ///< Returns the size of the last dimension in the format used by FFTW

    size_t data_size() const; ///< Returns data size assuming that its all stored as real numbers

    T& operator()(size_t i, size_t j, size_t k);      ///< Access to the voxel value assuming that `this` holds real data
    T operator()(size_t i, size_t j, size_t k) const; ///< Access to the voxel value assuming that `this` holds real data

  protected:
    std::shared_ptr< ImageSettings<T> > m_settings;   ///< Settings used to allocate memory and keeping configuration of relevant parameters
    
    T *m_data = nullptr;                              ///< Pointer to the allocated data
    std::array<size_t, NDIMS> m_n{0,0,0};             ///< Image dimenstions
    std::array<T, NDIMS> m_voxel;                     ///< Image voxel sizes

    FFTWPlan<T> m_plan_forward;                       ///< FFTW forward FFT plan
    FFTWPlan<T> m_plan_inverse;                       ///< FFTW inverse FFT plan
  };
}

#endif
