
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


#include "deconvolve.hpp"
#include "deconvolve_priv.hpp"
#include "constants.hpp"

#include <iostream>
#include <cassert>
#include <exception>
#include <functional>

namespace deconvolve {

  /// \brief Used to indicate failure to allocate private deconvolution implementation object
  class dec_bad_alloc: public std::bad_alloc {
  public:
    virtual const char* what() const throw() { return EXCPT_MEMORY "Failed to allocate memory for DeconvolvePrivate class"; }
  };
  
  ////////////////////////////////////////////////////////////////////////
  // Globally used API
  template <typename T> 
  Deconvolve<T>::Deconvolve()
  {
    m_dec.reset(new DeconvolvePrivate<T>());

    if (!m_dec)
      throw dec_bad_alloc();
  }

  template <typename T> 
  Deconvolve<T>::~Deconvolve()
  {
  }

  // NB! Internally, voxel sizes are expected in nanometers. Hence converting them from meters to nanometers while calling
  // internal methods
  
  template <typename T> 
  void Deconvolve<T>::set_psf(const std::vector<T> &data, size_t n1, size_t n2, size_t n3, T v1, T v2, T v3)
  {
    m_dec->set_psf(data, n1, n2, n3, v1*1e9, v2*1e9, v3*1e9);
  }

  template <typename T> 
  void Deconvolve<T>::set_callback(callback_type callback)
  {
    m_dec->set_callback(callback);
  }
  
  template <typename T> 
  void Deconvolve<T>::set_callback(callback_cpp_type callback)
  {
    m_dec->set_callback(callback);
  }
  
  template <typename T> 
  void Deconvolve<T>::set_callback(callback_extended_type callback, void *user_data)
  {
    using namespace std::placeholders;
    m_dec->set_callback( std::bind(callback, user_data, _1, _2, _3, _4, _5, _6, _7, _8, _9) );
  }
  
  template <typename T> 
  void Deconvolve<T>::clear_callback()
  {
    m_dec->clear_callback();
  }

  template <typename T> 
  void Deconvolve<T>::enable_regularization()
  {
    m_dec->enable_regularization();
  }
  
  template <typename T> 
  void Deconvolve<T>::disable_regularization()
  {
    m_dec->disable_regularization();
  }

  template <typename T> 
  bool Deconvolve<T>::regularized() const
  {
    return m_dec->regularized();
  }
  
  template <typename T> 
  void Deconvolve<T>::set_snr(T snr)
  {
    m_dec->set_snr(snr);
  }
  
  template <typename T> 
  void Deconvolve<T>::clear_snr()
  {
    m_dec->clear_snr();
  }

  template <typename T> 
  void Deconvolve<T>::set_max_iterations(size_t iters)
  {
    m_dec->set_max_iterations(iters);
  }
  
  template <typename T> 
  void Deconvolve<T>::clear_max_iterations()
  {
    m_dec->clear_max_iterations();
  }
  
  template <typename T> 
  size_t Deconvolve<T>::max_iterations() const
  {
    return m_dec->max_iterations();
  }

  template <typename T>
  void Deconvolve<T>::set_fftw_handlers( const typename fftw_implementation<T>::plan_function &forward,
                                         const typename fftw_implementation<T>::plan_function &inverse,
                                         const typename fftw_implementation<T>::clear_function &clear )
  {
    m_dec->set_fftw_handlers(forward, inverse, clear);
  }
  
  template <typename T>
  void Deconvolve<T>::clear_fftw_handlers()
  {
    m_dec->clear_fftw_handlers();
  }
  
  template <typename T> 
  std::vector<T> Deconvolve<T>::convolve(const std::vector<T> &data, size_t n1, size_t n2, size_t n3, T v1, T v2, T v3)
  {
    std::vector<T> result = data;
    m_dec->convolve(result, n1, n2, n3, v1*1e9, v2*1e9, v3*1e9);
    return result;
  }
  
  template <typename T> 
  std::vector<T> Deconvolve<T>::deconvolve(const std::vector<T> &data, size_t n1, size_t n2, size_t n3, T v1, T v2, T v3)
  {
    std::vector<T> result = data;
    m_dec->deconvolve(result, n1, n2, n3, v1*1e9, v2*1e9, v3*1e9);
    return result;
  }

  
  ////////////////////////////////////////
  // instantiate
  template class Deconvolve<float>;
  template class Deconvolve<double>;
}
