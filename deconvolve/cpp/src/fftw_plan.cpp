
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


#include "fftw_plan.hpp"
#include "constants.hpp"

#include <complex.h>
#include <complex>

#include <iostream>

#ifdef USE_FFTW_THREADS
 #include <omp.h>
#endif

#define FFTWFLAG FFTW_ESTIMATE

using namespace deconvolve;

template <typename T>
FFTWPlan<T>::FFTWPlan(std::shared_ptr< ImageSettings<T> > settings):
  m_settings(settings)
{
  if (!m_settings)
    throw std::runtime_error(EXCPT_INTERNAL "Cannot make FFT plans without any settings");
}

template <typename T>
FFTWPlan<T>::~FFTWPlan()
{
  clear();
}

template <typename T>
void FFTWPlan<T>::clear()
{
  m_settings->fftw_clear_plan(m_plan);
  m_plan = NULL;
}

template <typename T>
void FFTWPlan<T>::forward(T *data, int n0, int n1, int n2)
{
  clear();
  m_plan = m_settings->fftw_forward_plan(data, n0, n1, n2);
  if (!(*this))
    throw std::runtime_error(EXCPT_MEMORY "Couldn't allocate forward FFT plan");
}

template <typename T>
void FFTWPlan<T>::inverse(T *data, int n0, int n1, int n2)
{
  clear();
  m_plan = m_settings->fftw_inverse_plan(data, n0, n1, n2);
  if (!(*this))
    throw std::runtime_error(EXCPT_MEMORY "Couldn't allocate inverse FFT plan");
}


// specialization of members
namespace deconvolve {

  //////////
  // double

  template<>
  void FFTWPlan<double>::execute()
  {
    if (!(*this))
      throw std::runtime_error(EXCPT_INTERNAL "Cannot execute unallocated FFT plan");

    fftw_execute(m_plan);
  }

  //////////
  // float

  template<>
  void FFTWPlan<float>::execute()
  {
    if (!(*this))
      throw std::runtime_error(EXCPT_INTERNAL "Cannot execute unallocated FFT plan");

    fftwf_execute(m_plan);
  }

  
  // instantiate
  template class FFTWPlan<float>;
  template class FFTWPlan<double>;
}
