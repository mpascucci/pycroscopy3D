
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


#include "image_settings.hpp"

#include "constants.hpp"
#include "fftw_plan.hpp"

#include <iostream>
#include <mutex>

#ifdef USE_FFTW_THREADS
#include <omp.h>
#endif

#define FFTWFLAG FFTW_ESTIMATE


using namespace deconvolve;

////////////////////////////////////////////////////////////////////////////
// FFTW handling mutex, initialization state, and init function
namespace deconvolve {
  static std::mutex fftw_mutex;
  static bool fftw_initialized{false};

  /// \brief Initialize FFTW if using internal FFTW handling
  ///
  static void fftw_init()
  {
    if (fftw_initialized) return;

#ifdef USE_FFTW_THREADS
    if (fftw_init_threads() == 0)
      throw std::runtime_error(EXCPT_MEMORY "Failed to init FFTW threads");
#endif

    fftw_initialized = true;
  }
}


////////////////////////////////////////////////////////////////////////////
// ImageSettings implementation

// Constructors

template <typename T>
ImageSettings<T>::ImageSettings()
{
}

template <typename T>
ImageSettings<T>::ImageSettings(const ImageSettings<T> &old,
                                const typename fftw_implementation<T>::plan_function &forward,
                                const typename fftw_implementation<T>::plan_function &inverse,
                                const typename fftw_implementation<T>::clear_function &clear):
  ImageSettings(old, true)
{
  m_fftw_forward_plan = forward;
  m_fftw_inverse_plan = inverse;
  m_fftw_clear_plan = clear;
}


template <typename T>
ImageSettings<T>::ImageSettings(const ImageSettings<T> &old, bool increment_id)
{
  m_fftw_forward_plan = old.m_fftw_forward_plan;
  m_fftw_inverse_plan = old.m_fftw_inverse_plan;
  m_fftw_clear_plan = old.m_fftw_clear_plan;

  if (increment_id) m_id = old.m_id + 1;
  else m_id = old.m_id;
}

// template <typename T>
// void ImageSettings<T>::fftw_reset()
// {
//   m_fftw_forward_plan = fftw_plan_function();
//   m_fftw_inverse_plan = fftw_plan_function();
//   m_fftw_clear_plan = fftw_clear_function();
//   ++m_id;
// }

// FFTW handling
template <typename T>
typename fftw_implementation<T>::plan_type ImageSettings<T>::fftw_forward_plan(T *data, int n0, int n1, int n2)
{
  if (m_fftw_forward_plan) return m_fftw_forward_plan(data, n0, n1, n2);
  else
    {
      std::lock_guard<std::mutex> _lk(fftw_mutex);
      
      fftw_init();
      
#ifdef USE_FFTW_THREADS
      fftw_plan_with_nthreads(omp_get_max_threads());
#endif
      
      fftw_implementation_detail<T> fi;
      typename fftw_implementation<T>::plan_type plan = fi.forward(n0, n1, n2,
                                                                  data, (typename fftw_implementation<T>::complex_type*) data,
                                                                  FFTWFLAG);
      
      if (plan == NULL)
        throw std::runtime_error(EXCPT_MEMORY "Couldn't allocate forward FFT plan");
      
      return plan;
    }
}

template <typename T>
typename fftw_implementation<T>::plan_type ImageSettings<T>::fftw_inverse_plan(T *data, int n0, int n1, int n2)
{
  if (m_fftw_inverse_plan) return m_fftw_inverse_plan(data, n0, n1, n2);
  else
    {
      std::lock_guard<std::mutex> _lk(fftw_mutex);
      
      fftw_init();
      
#ifdef USE_FFTW_THREADS
      fftw_plan_with_nthreads(omp_get_max_threads());
#endif
      
      fftw_implementation_detail<T> fi;
      typename fftw_implementation<T>::plan_type plan = fi.inverse(n0, n1, n2,
                                                                  (typename fftw_implementation<T>::complex_type*)data, data,
                                                                  FFTWFLAG);
      if (plan == NULL)
        throw std::runtime_error(EXCPT_MEMORY "Couldn't allocate forward FFT plan");
      
      return plan;
    }
}

template <typename T>
void ImageSettings<T>::fftw_clear_plan(typename fftw_implementation<T>::plan_type plan)
{
  if (m_fftw_clear_plan) return m_fftw_clear_plan(plan);
  else
    {
      std::lock_guard<std::mutex> _lk(fftw_mutex);
      
      fftw_init();
      
      if (plan != NULL)
        {
          fftw_implementation_detail<T> fi;
          fi.clear(plan);
        }
    }
}
  

///////////////////////////////
// instantiate
namespace deconvolve {
  template class ImageSettings<double>;
  template class ImageSettings<float>;
}
