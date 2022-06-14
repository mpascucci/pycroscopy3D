
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


#include "image.hpp"

#include <complex.h>
#include <complex>
#include <fftw3.h>

#include <iostream>
#include <cassert>
#include <utility>
#include <exception>

#include <string.h> // memcpy
#include <cmath>    // sqrt

using namespace deconvolve;

template <typename T>
Image<T>::Image(std::shared_ptr< ImageSettings<T> > settings, const std::vector<T> &data, size_t n1, size_t n2, size_t n3, T v1, T v2, T v3):
  m_settings(settings),
  m_plan_forward(settings),
  m_plan_inverse(settings)
{
  set(data, n1, n2, n3, v1, v2, v3);
}

template <typename T>
Image<T>::Image(std::shared_ptr< ImageSettings<T> > settings, size_t n1, size_t n2, size_t n3, T v1, T v2, T v3):
  m_settings(settings),
  m_plan_forward(settings),
  m_plan_inverse(settings)
{
  set(std::vector<T>(), n1, n2, n3, v1, v2, v3);
}

template <typename T>
Image<T>::~Image()
{
  release_data();
}

namespace deconvolve {
  /// \brief Used to indicate failure to allocate FFTW data object.
  class fftw_bad_alloc: public std::bad_alloc {
  public:
    virtual const char* what() const throw() { return EXCPT_MEMORY "Failed to allocate memory for FFTW"; }
  };
}
      
template <typename T>
void Image<T>::allocate_data()
{
  assert(m_data == NULL);
  assert(m_n[0]>0 && m_n[1]>0 && m_n[2]>0);

  m_data = (T*)fftw_malloc(sizeof(T) * data_size());

  if (m_data == NULL)
    throw fftw_bad_alloc();
}

template <typename T>
void Image<T>::release_data()
{
  m_plan_forward.clear();
  m_plan_inverse.clear();

  fftw_free(m_data);
  m_data = nullptr;
  
  m_n.fill(0);
  m_voxel.fill(0);
}


template <typename T>
size_t Image<T>::last_dim() const
{
  return 2*(m_n[2]/2+1);
}

template <typename T>
size_t Image<T>::data_size() const
{
  return m_n[0]*m_n[1]*last_dim();
}

template <typename T>
bool Image<T>::same_dims(size_t n1, size_t n2, size_t n3) const
{
  return (n1==m_n[0] && n2==m_n[1] && n3==m_n[2]);
}

template <typename T>
bool Image<T>::same_voxel(T v1, T v2, T v3) const
{
  const T tol = 1e-13;
  return ( abs(v1-m_voxel[0])<tol && abs(v2-m_voxel[1])<tol && abs(v3-m_voxel[2])<tol );
}

template <typename T>
bool Image<T>::same_dims(const Image &im) const
{
  return same_dims(im.m_n[0], im.m_n[1], im.m_n[2]);
}

template <typename T>
bool Image<T>::same_voxel(const Image &im) const
{
  return same_voxel(im.m_voxel[0], im.m_voxel[1], im.m_voxel[2]);
}

template <typename T>
bool Image<T>::compatible(const Image &im) const
{
  return ( *(this) && im && 
           same_dims(im) && same_voxel(im) );
}

template <typename T>
void Image<T>::set(const std::vector<T> &data, size_t n1, size_t n2, size_t n3, T v1, T v2, T v3)
{
  if ( data.size()!=0 && data.size() != n1*n2*n3 )
    throw std::runtime_error(EXCPT_USER "Size of image data as represented by vector inconsistent with the given dimensions");
  
  release_data();

  m_n[0] = n1;
  m_n[1] = n2;
  m_n[2] = n3;

  m_voxel[0] = v1;
  m_voxel[1] = v2;
  m_voxel[2] = v3;

  allocate_data();

  if (data.size()!=0)
    {      
      // copy data over into FFTW format
      for (size_t i = 0; i < n1*n2; ++i)
        {
          T *d = m_data + i*last_dim();
          const T *s = data.data() + i*n3;
          for (size_t j=0; j < n3; ++j, ++d, ++s)
            *d = *s;
        }
    }
}

template <typename T>
void Image<T>::copy_data(const Image &im)
{
  if (!compatible(im) )
    throw std::runtime_error(EXCPT_INTERNAL "Trying to copy data between incompatible images");
  
  memcpy( (void*)m_data, (void*)im.m_data, data_size()*sizeof(T) );
}

template <typename T>
void Image<T>::get_image(std::vector<T> &data)
{
  if (!(*this))
    throw std::runtime_error(EXCPT_INTERNAL "Trying to get data from empty Image object");
  
  data.resize(m_n[0]*m_n[1]*m_n[2]);
  
  size_t n12 = m_n[0]*m_n[1];
  size_t n3_real = m_n[2];
  size_t n3 = last_dim();
  T *tgt = data.data();
  for (size_t i = 0; i < n12; ++i)
    {
      const T *d = m_data + i*n3;
      for (size_t j=0; j < n3_real; ++j, ++d, ++tgt)
        (*tgt) = (*d);
    }
}

template <typename T>
void Image<T>::swap(Image &image)
{
  using std::swap;

  swap(this->m_data, image.m_data);
  swap(this->m_n, image.m_n);
  swap(this->m_voxel, image.m_voxel);
  swap(this->m_plan_forward, image.m_plan_forward);
  swap(this->m_plan_inverse, image.m_plan_inverse);
}


template <typename T>
void Image<T>::fft()
{
  if (!(*this))
    throw std::runtime_error(EXCPT_INTERNAL "Trying to perform FFT on an empty Image object");

  if (!m_plan_forward)
    m_plan_forward.forward( m_data, m_n[0], m_n[1], m_n[2] );

  m_plan_forward.execute();
}

template <typename T>
void Image<T>::ifft()
{
  if (!(*this))
    throw std::runtime_error(EXCPT_INTERNAL "Trying to perform inverse FFT on an empty Image object");

  if (!m_plan_inverse)
      m_plan_inverse.inverse( m_data, m_n[0], m_n[1], m_n[2] );

  m_plan_inverse.execute();
}

template <typename T>
static void prod(std::complex<T> *im, std::complex<T> *ker, T scale, size_t sz)
{
  for (size_t i=0; i < sz; ++i, ++im, ++ker)
    (*im) *= (*ker) / scale;
}

template <typename T>
static void prod_conj(std::complex<T> *im, std::complex<T> *ker, T scale, size_t sz)
{
  for (size_t i=0; i < sz; ++i, ++im, ++ker)
    (*im) *= std::conj((*ker)) / scale;
}

template <typename T>
void Image<T>::convolve_implementation( const Image<T> &kernel,
                                        void (*p)(std::complex<T> *im, std::complex<T> *ker, T scale, size_t sz) )
{
  if (!compatible(kernel))
    throw std::runtime_error(EXCPT_INTERNAL "Convolution attempted with incompatible kernel");

  fft();

  T scale = m_n[0]*m_n[1]*m_n[2];
  std::complex<T> *im = (std::complex<T> *)m_data;
  std::complex<T> *ker = (std::complex<T> *)kernel.m_data;
  p(im, ker, scale, m_n[0]*m_n[1]*(m_n[2]/2+1));

  ifft();
}

template <typename T>
void Image<T>::convolve(const Image<T> &kernel)
{
  convolve_implementation(kernel, prod);
}

template <typename T>
void Image<T>::convolve_conj(const Image<T> &kernel)
{
  convolve_implementation(kernel, prod_conj);
}

template <typename T>
void Image<T>::invdivide_image(const Image<T> &image)
{
  if ( !compatible(image) )
    throw std::runtime_error(EXCPT_INTERNAL "invDivide attempted between incompatible images");
  
  size_t n12 = m_n[0]*m_n[1];
  size_t n3_real = m_n[2];
  size_t n3 = last_dim();
  
#pragma omp parallel for
  for (size_t i = 0; i < n12; ++i)
    {
      T *d = m_data + i*n3;
      const T *im = image.m_data + i*n3;
      for (size_t j=0; j < n3_real; ++j, ++d, ++im)
        {
          const T v = *d;
          if (v <= 0) *d = 0;
          else (*d) = (*im) / v;
        }
    }
}

template <typename T>
void Image<T>::prod_image(const Image<T> &image)
{
  if ( !compatible(image) )
    throw std::runtime_error(EXCPT_INTERNAL "prod_image attempted between incompatible images");

  size_t n12 = m_n[0]*m_n[1];
  size_t n3_real = m_n[2];
  size_t n3 = last_dim();
  
#pragma omp parallel for
  for (size_t i = 0; i < n12; ++i)
    {
      T *d = m_data + i*n3;
      const T *im = image.m_data + i*n3;
      for (size_t j=0; j < n3_real; ++j, ++d, ++im)
        (*d) = (*im) * (*d);
    }
}


template <typename T>
void Image<T>::prod_regularized(const Image<T> &image, T lambda, const Image<T> &div)
{
  if ( !compatible(image) )
    throw std::runtime_error(EXCPT_INTERNAL "prod_regularized attempted between incompatible images");

  size_t n12 = m_n[0]*m_n[1];
  size_t n3_real = m_n[2];
  size_t n3 = last_dim();
  
#pragma omp parallel for
  for (size_t i = 0; i < n12; ++i)
    {
      T *result = m_data + i*n3;
      const T *im = image.m_data + i*n3;
      const T *di = div.m_data + i*n3;
      for (size_t j=0; j < n3_real; ++j, ++result, ++im, ++di)
        (*result) = (*result) * (*im) / (1.0 - lambda*(*di));
    }
}


template <typename T>
T& Image<T>::operator()(size_t i, size_t j, size_t k)
{
  return *(m_data + i*m_n[1]*last_dim() + j*last_dim() + k);
}

template <typename T>
T Image<T>::operator()(size_t i, size_t j, size_t k) const
{
  return *(m_data + i*m_n[1]*last_dim() + j*last_dim() + k);
}

template <typename T>
static T m(T a, T b)
{
  if (a<0 && b<0)
    {
      if (a >= b) return a;
      return b;
    }
  if (a>0 && b>0)
    {
      if (a < b) return a;
      return b;
    }
  return 0.0;
}

template <typename T>
static T hypot3(T a, T b, T c)
{
  return std::sqrt(a*a + b*b + c*c);
}

template <typename T>
void Image<T>::div_unit_grad(const Image<T> &image)
{
  if ( !compatible(image) )
    throw std::runtime_error(EXCPT_INTERNAL "div_unit_grad attempted between incompatible images");

  const T h0 = image.m_voxel[0];
  const T h1 = image.m_voxel[1];
  const T h2 = image.m_voxel[2];
  const T eps = 0.0;
  
  size_t n1 = m_n[0];
  size_t n2 = m_n[1];
  size_t n3 = m_n[2];

  #pragma omp parallel for
  for (size_t i=0; i<n1; ++i)
    {
      size_t im1 = (i?i-1:0);
      size_t im2 = (im1?im1-1:0);
      size_t ip1 = (i+1==m_n[0]?i:i+1);

      for (size_t j=0; j<n2; ++j)
        {
          size_t jm1 = (j?j-1:0);
          size_t jm2 = (jm1?jm1-1:0);
          size_t jp1 = (j+1==m_n[1]?j:j+1);
          
          for (size_t k=0; k<n3; ++k)
            {
              size_t km1 = (k?k-1:0);
              size_t km2 = (km1?km1-1:0);
              size_t kp1 = (k+1==m_n[2]?k:k+1);
              
              T fimjm = image(im1, jm1, k);
              T fim = image(im1, j, k);
              T fimkm = image(im1, j, km1);
              T fimkp = image(im1, j, kp1);
              T fimjp = image(im1, jp1, k);

              T fjmkm = image(i, jm1, km1);
              T fjm = image(i, jm1, k);
              T fjmkp = image(i, jm1, kp1);

              T fkm = image(i, j, km1);
              T fijk = image(i, j, k);
              T fkp = image(i, j, kp1);

              T fjpkm = image(i, jp1, km1);
              T fjp = image(i, jp1, k);

              T fipjm = image(ip1, jm1, k);
              T fipkm = image(ip1, j, km1);
              T fip = image(ip1, j, k);

              T Dxpf = (fip - fijk) / h0;
              T Dxmf = (fijk - fim) / h0;
              T Dypf = (fjp - fijk) / h1;
              T Dymf = (fijk - fjm) / h1;
              T Dzpf = (fkp - fijk) / h2;
              T Dzmf = (fijk - fkm) / h2;
              T aijk = hypot3(Dxpf, m(Dypf, Dymf), m(Dzpf, Dzmf));
              T bijk = hypot3(Dypf, m(Dxpf, Dxmf), m(Dzpf, Dzmf));
              T cijk = hypot3(Dzpf, m(Dypf, Dymf), m(Dxpf, Dxmf));

              aijk = (aijk>eps?Dxpf / aijk:0.0);
              bijk = (bijk>eps?Dypf / bijk:0.0);
              cijk = (cijk>eps?Dzpf / cijk:0.0); 
		  

              Dxpf = (fijk - fim) / h0;
              Dypf = (fimjp - fim) / h1;
              Dymf = (fim - fimjm) / h1;
              Dzpf = (fimkp - fim) / h2;
              Dzmf = (fim - fimkm) / h2;
              T aim = hypot3(Dxpf, m(Dypf, Dymf), m(Dzpf, Dzmf));

              aim = (aim>eps?Dxpf/aim:0.0); 


              Dxpf = (fipjm - fjm) / h0;
              Dxmf = (fjm - fimjm) / h0;
              Dypf = (fijk - fjm) / h1;
              Dzpf = (fjmkp - fjm) / h2;
              Dzmf = (fjm - fjmkm) / h2;
              T bjm = hypot3(Dypf, m(Dxpf, Dxmf), m(Dzpf, Dzmf));

              bjm = (bjm>eps?Dypf/bjm:0.0);
		  

              Dxpf = (fipkm - fkm) / h0;
              Dxmf = (fjm - fimkm) / h0;
              Dypf = (fjpkm - fkm) / h1;
              Dymf = (fkm - fjmkm) / h1;
              Dzpf = (fijk - fkm) / h2;
              T ckm = hypot3(Dzpf, m(Dypf, Dymf), m(Dxpf, Dxmf));

              ckm = (ckm>eps?Dzpf/ckm:0.0); 

              T Dxma = (aijk - aim) / h0;
              T Dymb = (bijk - bjm) / h1;
              T Dzmc = (cijk - ckm) / h2;

              this->operator()(i, j, k) = Dxma + Dymb + Dzmc;
            }
        }
    }
}


template <typename T>
T Image<T>::snr(size_t convolution_kernel_size) const
{
  if (!(*this))
    throw std::runtime_error(EXCPT_INTERNAL "Cannot determine SNR of an empty image");
  
  T snr = 0.0;
  
  size_t n1 = m_n[0];
  size_t n2 = m_n[1];
  size_t n3 = m_n[2];

  int ker = convolution_kernel_size;

#pragma omp parallel for reduction(max:snr) 
  for (size_t i1 = convolution_kernel_size; i1 < n1-convolution_kernel_size; ++i1)
    for (size_t i2 = convolution_kernel_size; i2 < n2-convolution_kernel_size; ++i2)
      for (size_t i3 = convolution_kernel_size; i3 < n3-convolution_kernel_size; ++i3)
        {
          T s = 0.0; //this->operator()( i1, i2, i3 );
          for (int j1=-ker; j1 <= ker; ++j1)
            for (int j2=-ker; j2 <= ker; ++j2)
              for (int j3=-ker; j3 <= ker; ++j3)
                s += this->operator()( i1+j1, i2+j2, i3+j3 );

          snr = std::max(s, snr);
        }

  return std::sqrt( snr / std::pow(convolution_kernel_size*2+1, 3) );
}

template <typename T>
void Image<T>::get_stats(T &cmin, T &cmax, T &csum) const
{
  if (!(*this))
    throw std::runtime_error(EXCPT_INTERNAL "Cannot determine image statistics of an empty image");
  
  csum = 0;
  cmax = cmin = m_data[0];

  size_t n12 = m_n[0]*m_n[1];
  size_t n3_real = m_n[2];
  size_t n3 = last_dim();

#pragma omp parallel for reduction(+:csum) reduction(min: cmin) reduction(max:cmax) 
  for (size_t i = 0; i < n12; ++i)
    {
      const T *d = m_data + i*n3;
      for (size_t j=0; j < n3_real; ++j, ++d)
        {
          const T v = *d;
          cmin = std::min(v, cmin);
          cmax = std::max(v, cmax);
          csum += v;
        }
    }
}

template <typename T>
T Image<T>::nrm2(const Image &image) const
{
  if (!(*this))
    throw std::runtime_error(EXCPT_INTERNAL "Cannot determine image norm of an empty image");
  
  T nrm = 0.0;
  
  size_t n12 = m_n[0]*m_n[1];
  size_t n3_real = m_n[2];
  size_t n3 = last_dim();
  
#pragma omp parallel for reduction(+:nrm) 
  for (size_t i = 0; i < n12; ++i)
    {
      const T *d = m_data + i*n3;
      const T *d2 = image.m_data + i*n3;
      for (size_t j=0; j < n3_real; ++j, ++d, ++d2)
        {
          const T t = *d - *d2;
          nrm += t*t;
        }
    }

  return nrm;
}


template <typename T>
T Image<T>::lambda_lsq(const Image &cconv, const Image &div)
{
  if (!cconv.compatible(div))
    throw std::runtime_error(EXCPT_INTERNAL "Cannot determine lambda for incompatible or empty images");
  
  T divsqrsum = 0.0;
  T lambda = 0.0;

  size_t n12 = cconv.m_n[0]*cconv.m_n[1];
  size_t n3_real = cconv.m_n[2];
  size_t n3 = cconv.last_dim();

#pragma omp parallel for reduction(+:lambda,divsqrsum) 
  for (size_t i = 0; i < n12; ++i)
    {
      const T *cc = cconv.m_data + i*n3;
      const T *d  = div.m_data   + i*n3;
      for (size_t j=0; j < n3_real; ++j, ++d, ++cc)
        {
          T cval = (*cc);
          T dval = (*d);

          lambda += (1.0 - cval)*dval;
          divsqrsum += dval*dval;
        }
    }
  
  return lambda / divsqrsum;
}


///////////////////////////////
// instantiate
namespace deconvolve {
  template class Image<float>;
  template class Image<double>;
}
