
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


#ifndef IOCBIO_PSF_HPP
#define IOCBIO_PSF_HPP

#include "image.hpp"
#include "constants.hpp"

#include <array>
#include <vector>
#include <memory>

namespace deconvolve {

  /// \brief Point spread function
  ///
  /// Object of this class keeps point spread function (PSF) and
  /// provides the facility to calculate OTF after scaling PSF into
  /// the image-based dimensions and voxel sizes.
  template <typename T>
  class PSF {
  public:

    PSF();

    /// \brief Construct PSF from the user-provided data
    ///
    /// See \ref set for description
    ///
    PSF(const std::vector<T> &data, size_t n1, size_t n2, size_t n3, T v1, T v2, T v3);
    
    ~PSF();

    /// \brief Set PSF data
    ///
    /// Set PSF on the basis of provided data, its dimensions, and voxel sizes.
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

    /// \brief Calculate and return OTF of the PSF
    ///
    /// Calculate OTF for the PSF after scaling it to fill the given
    /// dimensions taking into account the voxel sizes. The calculated
    /// OTF is kept in the memory by PSF object and is returned from
    /// this cache if it is requested again with the same parameters.
    ///
    /// OTF is recalculated if a new PSF is given by \ref set, image
    /// settings were changed or some other image parameters
    /// (dimension or voxel size) have changed.
    ///
    /// \param settings image settings describing memory storage options and operations
    /// \param n1 the slowest changing dimension for OTF
    /// \param n2 the medium changing dimension for OTF
    /// \param n3 the fastest changing dimension for OTF
    /// \param v1 voxel size along dimension 1, in meters (for OTF)
    /// \param v2 voxel size along dimension 2, in meters (for OTF)
    /// \param v3 voxel size along dimension 3, in meters (for OTF)
    ///    
    Image<T>& otf(std::shared_ptr< ImageSettings<T> > settings, size_t n1, size_t n2, size_t n3, T v1, T v2, T v3);

    /// \brief Test if PSF has been allocated
    ///
    /// \return `true` if `this` contains PSF data
    operator bool() const { return m_data.size() > 0; }
    
  protected:
    std::vector<T> m_data;                ///< PSF voxel values
    std::array<size_t, NDIMS> m_n{0,0,0}; ///< PSF dimensions
    std::array<T, NDIMS> m_voxel;    ///< PSF voxel sizes

    std::unique_ptr< Image<T> > m_otf;    ///< OTF that was calculated earlier
  };
  
}

#endif
