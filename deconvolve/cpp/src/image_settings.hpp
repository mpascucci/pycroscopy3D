
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


#ifndef IOCBIO_IMAGE_SETTINGS_HPP
#define IOCBIO_IMAGE_SETTINGS_HPP

#include "fftw_interface.hpp"

#include <functional>
#include <memory>

namespace deconvolve {

  template <typename T> struct fftw_implementation;

  /// \brief Image settings
  ///
  /// This class keeps relevant image settings that are expected to
  /// describe memory storage model, FFT operations, and such. At
  /// present, it is holding FFTW plan handlers. In future, if more
  /// backends are going to be implemented, this class can be used to
  /// specify the backend as well as carry the backend-specific
  /// settings, as done for FFTW now.
  ///
  /// To avoid static variables, ImageSettings can be set only through
  /// constructors by supplying the base settings (old settings
  /// object) and the new settings. Internally, in the library,
  /// ImageSettings objects are created so that the comparison between
  /// the settings can be performed just by comparing the internal
  /// counter `m_id` which is incremented every time the new settings
  /// are created on the basis of the old ones.
  ///
  template <typename T>
  class ImageSettings {

  public:
    
  public:
    ImageSettings(); ///< Default constructor used for the first initialization

    /// \brief Constructor with new FFTW handlers
    ///
    /// Constructs the new settings with supplied FFTW handling
    /// routines. The FFTW plan handlers can be empty functions. In
    /// that case, the default thread-safe implementation of the
    /// handlers is used. See \ref Deconvolve::set_fftw_handlers for
    /// the description of used API.
    ///
    /// \param old Settings to base the new settings on
    /// \param forward function returning forward FFTW plan
    /// \param inverse function returning inverse FFTW plan
    /// \param clear function to destroy FFTW plan
    ///
    /// \sa Deconvolve::set_fftw_handlers
    ///
    ImageSettings(const ImageSettings &old,
                  const typename fftw_implementation<T>::plan_function &forward,
                  const typename fftw_implementation<T>::plan_function &inverse,
                  const typename fftw_implementation<T>::clear_function &clear ); 

    /// \brief Check if the settings are the same as the ones in the argument
    bool same(const ImageSettings &other) { return other.m_id == m_id; }

    /// \brief Check if the settings are the same as the ones in the argument
    bool same(const std::shared_ptr<ImageSettings> other) { return other && this->same(*other); }

    /// \brief Create FFTW plan for FFT
    ///
    /// \sa Deconvolve::set_fftw_handlers
    typename fftw_implementation<T>::plan_type fftw_forward_plan(T *data, int n0, int n1, int n2);

    /// \brief Create FFTW plan for iFFT
    ///
    /// \sa Deconvolve::set_fftw_handlers
    typename fftw_implementation<T>::plan_type fftw_inverse_plan(T *data, int n0, int n1, int n2);

    /// \brief Clear (destroy) allocated FFTW plan
    ///
    /// \sa Deconvolve::set_fftw_handlers
    void fftw_clear_plan(typename fftw_implementation<T>::plan_type plan);

  protected:    
    ImageSettings(const ImageSettings &old, bool increment_id); ///< Copy old settings before changing them (used internally)

  protected:
    /// \brief id of the current settings
    ///
    /// `m_id` is incremented every time new settings are generated on
    /// the basis of the older settings.
    size_t m_id{0};

    typename fftw_implementation<T>::plan_function m_fftw_forward_plan; ///< Current handler for FFTW plan creation. If not specified, a default handler is used
    typename fftw_implementation<T>::plan_function m_fftw_inverse_plan; ///< Current handler for FFTW plan creation. If not specified, a default handler is used
    typename fftw_implementation<T>::clear_function m_fftw_clear_plan; ///< Current handler for FFTW plan destruction. If not specified, a default handler is used
  };

}

#endif

