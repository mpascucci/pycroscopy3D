#!/usr/bin/env python3

#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#  
#   Copyright (C) 2018-2020
#    Laboratory of Systems Biology, Department of Cybernetics,
#    School of Science, Tallinn University of Technology
#   This file is part of project: IOCBIO Deconvolve



import numpy as np
import pylab as plt
from scipy import ndimage
from scipy.interpolate import RegularGridInterpolator
from simple_plot import read_data


def gaussian_3d(shape, sigma, voxel, a=1.):
    # a : amplitude
    Z, Y, X = np.meshgrid(np.arange(shape[0]),
                          np.arange(shape[1]),
                          np.arange(shape[2]),
                          indexing='ij')
    sz, sy, sx = np.array(sigma) / np.array(voxel)
    z0, y0, x0 = [int(np.floor(i/2)) for i in shape]

    return a*np.exp( -( ((X-x0)/sx)**2 + ((Y-y0)/sy)**2 + ((Z-z0)/sz)**2) )


if __name__ == '__main__':

    # shape = (41, 21, 21)
    # voxel = (0.3, 0.1, 0.1)
    shape = (81, 41, 41)
    voxel = (0.15, 0.05, 0.05)

    sigma = (0.8, 0.2, 0.2)
    z0, y0, x0 = [int(np.floor(i/2)) for i in shape]

    psf = gaussian_3d(shape, sigma, voxel)
    print('psf', psf.shape, psf.min(), psf.max())

    # nshape = (41, 128, 128)
    # nvoxel = (0.3, 0.15, 0.15)
    nshape = (120, 128, 128)
    nvoxel = (0.05, 0.01, 0.01)

    nz0, ny0, nx0 = [int(np.floor(i/2)) for i in nshape]

    # scipy.ndimage.interpolation.zoom
    if 0:
        zpsf = ndimage.interpolation.zoom(psf, np.array(voxel)/np.array(nvoxel), order=1)

    # scipy.interpolate.RegularGridInterpolator
    if 1:
        trilin_func = RegularGridInterpolator((np.arange(shape[0])-z0, np.arange(shape[1])-y0, np.arange(shape[2])-x0), psf)
        zshape = [int(i) for i in np.array(voxel) / np.array(nvoxel) * np.array(shape)]
        z = np.linspace(0, shape[0]-1, zshape[0])-z0
        y = np.linspace(0, shape[1]-1, zshape[1])-y0
        x = np.linspace(0, shape[2]-1, zshape[2])-x0
        pts = []
        for k in z:
            for l in y:
                for m in x:
                    pts.append([k,l,m])

        zpsf = trilin_func(pts)
        zpsf.resize(zshape)

    zz0, zy0, zx0 = [int(np.floor(i/2)) for i in zpsf.shape]
    print('zpsf', zpsf.shape, zpsf.min(), zpsf.max(), zpsf.sum())

    if 0:
        npsf = np.zeros(nshape)
        npsf[nz0-zz0:nz0+zpsf.shape[0]-zz0, ny0-zy0:ny0+zpsf.shape[1]-zy0, nx0-zx0:nx0+zpsf.shape[2]-zx0] = zpsf
    if 1:
        npsf = zpsf[zz0-60:zz0+60, zy0+1-64:zy0+1+64, zx0+1-64:zx0+1+64]

    fig = plt.figure()
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)

    ax1.imshow(psf[z0,:,:], interpolation='nearest')
    ax1.set_title('XY orig')

    ax2.imshow(psf[:,y0,:], interpolation='nearest')
    ax2.set_title('ZX orig')

    fname = '/export/kybi/home/markov/code/deconvolve/cpp/psf_otf.data'
    data = read_data(fname)
    dif = npsf-data

    print('npsf', npsf.shape, npsf.min(), npsf.max(), npsf.sum(), ndimage.measurements.maximum_position(npsf))
    print('data', data.shape, data.min(), data.max(), data.sum(), ndimage.measurements.maximum_position(data))
    print('dif', dif.shape, dif.min(), dif.max(), dif.sum())


    ax3.imshow(dif[nz0,:,:], interpolation='nearest')
    ax3.set_title('XY new')

    ax4.imshow(dif[:,ny0,:], interpolation='nearest')
    ax4.set_title('ZX new')

    #
    # ax3.imshow(npsf[nz0,:,:], interpolation='nearest')
    # ax3.set_title('XY new')
    #
    # ax4.imshow(npsf[:,ny0,:], interpolation='nearest')
    # ax4.set_title('ZX new')

    plt.show()
