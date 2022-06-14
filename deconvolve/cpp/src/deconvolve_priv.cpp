
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


#include "deconvolve_priv.hpp"
#include "constants.hpp"

#include <exception>
#include <iostream>

using namespace deconvolve;

// private class implementation
template <typename T>
DeconvolvePrivate<T>::DeconvolvePrivate():
    m_settings(new ImageSettings<T>())
{
}

template <typename T>
DeconvolvePrivate<T>::~DeconvolvePrivate()
{
}

template <typename T>
void DeconvolvePrivate<T>::set_psf(const std::vector<T> &data, size_t n1, size_t n2, size_t n3, T v1, T v2, T v3)
{
  m_psf.set(data, n1, n2, n3, v1, v2, v3);  
}

template <typename T>
void DeconvolvePrivate<T>::set_callback(callback_cpp_type callback)
{
  m_callback = callback;
}

template <typename T>
void DeconvolvePrivate<T>::clear_callback()
{
  m_callback = callback_type();
}

template <typename T>
void DeconvolvePrivate<T>::enable_regularization()
{
  m_regularize = true;
}

template <typename T>
void DeconvolvePrivate<T>::disable_regularization()
{
  m_regularize = false;
}

template <typename T>
void DeconvolvePrivate<T>::set_snr(T snr)
{
  m_snr = snr;
}

template <typename T>
void DeconvolvePrivate<T>::clear_snr()
{
  m_snr = -1;
}

template <typename T>
void DeconvolvePrivate<T>::set_fftw_handlers( const typename fftw_implementation<T>::plan_function &forward,
                                              const typename fftw_implementation<T>::plan_function &inverse,
                                              const typename fftw_implementation<T>::clear_function &clear )
{
  m_settings.reset(new ImageSettings<T>(*m_settings, forward, inverse, clear));
}
  
template <typename T>
void DeconvolvePrivate<T>::clear_fftw_handlers()
{
  m_settings.reset(new ImageSettings<T>(*m_settings,
                                        typename fftw_implementation<T>::plan_function(),
                                        typename fftw_implementation<T>::plan_function(),
                                        typename fftw_implementation<T>::clear_function() ));
}

template <typename T>
void DeconvolvePrivate<T>::convolve(std::vector<T> &data, size_t n1, size_t n2, size_t n3, T v1, T v2, T v3)
{
  if (!m_psf)
    throw std::runtime_error(EXCPT_USER "Cannot convolve without PSF. Please set PSF before calling convolve");

  Image<T> image(m_settings, data, n1, n2, n3, v1, v2, v3);
  Image<T> &otf = m_psf.otf(m_settings, n1, n2, n3, v1, v2, v3);

  image.convolve(otf);
  image.get_image(data);
}


template <typename T>
void DeconvolvePrivate<T>::deconvolve(std::vector<T> &data, size_t n1, size_t n2, size_t n3, T v1, T v2, T v3)
{
  if (!m_psf)
    throw std::runtime_error(EXCPT_USER "Cannot deconvolve without PSF. Please set PSF before calling deconvolve");

  Image<T> image(m_settings, data, n1, n2, n3, v1, v2, v3);
  Image<T> &otf = m_psf.otf(m_settings, n1, n2, n3, v1, v2, v3);

  // temporary images
  Image<T> oC(m_settings, data, n1, n2, n3, v1, v2, v3); // current iteration
  Image<T> o0(m_settings, data, n1, n2, n3, v1, v2, v3); // previous iteration
  Image<T> om1(m_settings, n1, n2, n3, v1, v2, v3);      // 2 iterations ago
  Image<T> div(m_settings, n1, n2, n3, v1, v2, v3);      // divergence

  // clear lambda stack
  m_lambda_evolution.clear();

  // estimate snr
  T snr = m_snr;
  if (snr < 0) // not specified, have to calcuate
    snr = oC.snr(1);

  // first estimation is convolved original image
  oC.convolve(otf);

  // iteration
  T lambda_factor = -1;
  T lambda = 0;
  T cmin = 0, cmax = 0, csum = 0, nrm2_prev = 0, nrm2_prevprev = 0;                   
  for (size_t iter = 0;
       (m_callback &&
        m_callback(iter, cmin, cmax, csum, nrm2_prev, nrm2_prevprev,
                   lambda, lambda_factor, snr)) ||
         (!m_callback &&
          callback_default(iter, cmin, cmax, csum, nrm2_prev, nrm2_prevprev,
                           lambda, lambda_factor, snr));
       ++iter)
    {
      oC.convolve(otf);
      oC.invdivide_image(image);
      oC.convolve_conj(otf);

      if ( !m_regularize )
        oC.prod_image(o0);

      else
        {
          div.div_unit_grad(o0);
          
          lambda = Image<T>::lambda_lsq(oC, div);

          if (lambda < 0 && iter == 0)
            throw std::runtime_error(EXCPT_NOBODYS_FAULT " First estimate of regularization factor is negative, cannot continue "
                                     "(lambda = " +
                                     std::to_string(lambda) + ")");
          
          if (iter == 0)
            lambda_factor = 50 / snr / lambda;

          if (lambda < 0)
            lambda = 0.0;
          
          lambda *= lambda_factor;
      
          oC.prod_regularized(o0, lambda, div);
        }

      oC.get_stats(cmin, cmax, csum);

      nrm2_prev = oC.nrm2(o0);
      if (iter > 1)
        nrm2_prevprev = oC.nrm2(om1);

      om1.swap(o0);
      o0.copy_data(oC);
    }

  oC.get_image(data);
}


template <typename T>
int DeconvolvePrivate<T>::callback_default(size_t iter,
                                           double min, double max, double sum,
                                           double nrm2_prev, double nrm2_prevprev,
                                           double lambda, double lambda_factor, double snr)
{
  using namespace std;

  // check if lambda is decaying if the stack is full
  bool done = false;
  if (m_regularize && m_lambda_evolution.size() >= const_lambda_stack_size)
    {
      done = true;
      for (const auto i: m_lambda_evolution)
        done = (done && i > lambda);
    }

  cout << "Iter: " << iter << " " << "Min/Max/Sum: " << min << " " << max << " " << sum
       << "  Nrm2 (i)-(i-1)/(i)-(i-2): " << nrm2_prev << " " << nrm2_prevprev
       << "  Lambda: " << lambda
       << "  LFactor: " << lambda_factor << "  SNR: " << snr << endl;

  m_lambda_evolution.push_back(lambda);
  if (m_lambda_evolution.size() > const_lambda_stack_size)
    m_lambda_evolution.pop_front();
  
  return !done && (iter < m_max_iterations);
}


////////////////////////////////////////
// instantiate
namespace deconvolve {
  template class DeconvolvePrivate<float>;
  template class DeconvolvePrivate<double>;
}
