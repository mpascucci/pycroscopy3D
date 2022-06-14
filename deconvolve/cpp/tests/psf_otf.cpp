
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


#include <iostream>
#include <fstream>
#include <vector>

#include <math.h>

#include "psf.hpp"

int main()
{
  deconvolve::PSF<double> psf;

  std::vector<double> psf_orig;
  size_t n12=41;
  size_t n3=81;
  double v12 = 0.05;
  double v3 = 0.15;

  for (int i1=-20; i1 < 21; ++i1)
    for (int i2=-20; i2 < 21; ++i2)
      for (int i3=-40; i3 < 41; ++i3)
        psf_orig.push_back( exp( -pow(i1*v12/0.2,2) ) *
                            exp( -pow(i2*v12/0.2,2) ) *
                            exp( -pow(i3*v3/0.8,2) ) );
    
  psf.set( psf_orig.data(), n12, n12, n3, v12, v12, v3 );

  size_t s12 = 128, s3=120;
  deconvolve::Image<double> &image = psf.otf(s12, s12, s3, 0.01, 0.01, 0.05);

  std::vector<double> interp;
  if ( !image.get_image(interp) )
    {
      std::cerr << "Failed to get image" << std::endl;
      return -1;
    }

  std::ofstream fout("psf_otf.data");
  fout << s12 << " " << s12 << " " << s3 << " ";
  for (size_t i=0; i < interp.size(); ++i)
    fout << interp[i] << " ";
  
  return 0;
}
