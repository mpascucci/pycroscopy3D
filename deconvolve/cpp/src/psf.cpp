
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


#include "psf.hpp"
#include "boost/multi_array.hpp"

#include <algorithm>
#include <exception>
#include <numeric>

using namespace deconvolve;
using namespace std;

template <typename T>
PSF<T>::PSF()
{
}

template <typename T>
PSF<T>::PSF(const std::vector<T> &data, size_t n1, size_t n2, size_t n3, T v1, T v2, T v3)
{
  set(data, n1, n2, n3, v1, v2, v3);
}

template <typename T>
PSF<T>::~PSF()
{
}


template <typename T>
void PSF<T>::set(const std::vector<T> &data, size_t n1, size_t n2, size_t n3, T v1, T v2, T v3)
{
  if ( data.size()!=0 && data.size() != n1*n2*n3 )
    throw std::runtime_error(EXCPT_USER "Size of PSF data as represented by vector inconsistent with the given dimensions");

  m_n[0] = n1;
  m_n[1] = n2;
  m_n[2] = n3;

  m_voxel[0] = v1;
  m_voxel[1] = v2;
  m_voxel[2] = v3;

  m_data = data;
  m_otf.reset();
}

template <typename T>
static T getind(T distance, T voxel, size_t elem)
{
  return distance / voxel + elem*0.5 - 0.5;
}

template <typename T>
Image<T>& PSF<T>::otf(std::shared_ptr< ImageSettings<T> > settings,  size_t n1, size_t n2, size_t n3, T v1, T v2, T v3)
{
  if ( m_data.size() < 1 )
    throw std::runtime_error(EXCPT_INTERNAL "Requesting OTF from empty PSF");
  
  if (m_otf &&
      m_otf->same_settings(settings) && 
      m_otf->same_dims(n1, n2, n3) &&
      m_otf->same_voxel(v1, v2, v3) )
    return *m_otf;

  std::vector<T> psf_interp_data(n1*n2*n3, 0);
  boost::multi_array_ref<T, NDIMS> psf_interp(psf_interp_data.data(), boost::extents[n1][n2][n3]);
  boost::multi_array_ref<T, NDIMS> psf(m_data.data(), m_n);

  // tri-linear interpolation of PSF data into OTF_DATA (before OTF is calculated by FFT)
#pragma omp parallel for
  for (int i1=0; i1 < n1; ++i1)
    {
      T d1 = v1 * ((i1+0.5) - n1*0.5);
      T ind1 = getind<>(d1, m_voxel[0], m_n[0]);
      int j1 = (int)ind1;
      T x1 = ind1-j1;

      for (int i2=0; i2 < n2; ++i2)
        {
          T d2 = v2 * ((i2+0.5) - n2*0.5);
          T ind2 = getind<>(d2, m_voxel[1], m_n[1]);
          int j2 = (int)ind2;
          T x2 = ind2-j2;

          for (int i3=0; i3 < n3; ++i3)
            {
              T d3 = v3 * ((i3+0.5) - n3*0.5);
              T ind3 = getind<>(d3, m_voxel[2], m_n[2]);
              int j3 = (int)ind3;
              T x3 = ind3-j3;

              T interp_value = 0;

              // check for out of bounds
              if ( j1 >= 0 && j1 < m_n[0]-1 &&
                   j2 >= 0 && j2 < m_n[1]-1 &&
                   j3 >= 0 && j3 < m_n[2]-1 )
                {
                  // https://en.wikipedia.org/wiki/Trilinear_interpolation
                  T c000 = psf[j1][j2][j3];
                  T c001 = psf[j1][j2][j3+1];
                  T c010 = psf[j1][j2+1][j3];
                  T c011 = psf[j1][j2+1][j3+1];
                  T c100 = psf[j1+1][j2][j3];
                  T c101 = psf[j1+1][j2][j3+1];
                  T c110 = psf[j1+1][j2+1][j3];
                  T c111 = psf[j1+1][j2+1][j3+1];

                  T c00 = c000*(1-x1) + c100*x1;
                  T c01 = c001*(1-x1) + c101*x1;
                  T c10 = c010*(1-x1) + c110*x1;
                  T c11 = c011*(1-x1) + c111*x1;

                  T c0 = c00*(1-x2) + c10*x2;
                  T c1 = c01*(1-x2) + c11*x2;

                  T c = c0*(1-x3) + c1*x3;

                  interp_value = c;
                }

              // shift PSF to enshure that the convolved image will not
              // move around
              size_t s1 = (i1 + n1/2 + 1) % n1;
              size_t s2 = (i2 + n2/2 + 1) % n2;
              size_t s3 = (i3 + n3/2 + 1) % n3;

              psf_interp[s1][s2][s3] = interp_value;
            }
        }
    }

  // normalize psf sum to 1.0
  T s = std::accumulate(psf_interp_data.begin(), psf_interp_data.end(), 0.0);
  for (T &v: psf_interp_data)
    v /= s;

  m_otf.reset(new Image<T>(settings, psf_interp_data, n1, n2, n3, v1, v2, v3));
  m_otf->fft();
  return *m_otf;
}



// instantiate
namespace deconvolve {
  template class PSF<float>;
  template class PSF<double>;
}
